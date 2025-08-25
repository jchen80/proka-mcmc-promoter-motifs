"""
Author: Julie Chen
calculate_kmer_stats.py

For a given species directory containing the output of mcmc/cds_randomizer.py, 
which should be a *kmer_frequencies_k{6/7}.csv file, with the format column 0 (kmer sequence)
column 1 (original genome frequencies), columns 2-onward (randomized genome frequencies)

then, will (1) generate a file with statistics per kmer of enrichment score (ratio) and 
significance of enrichment relative to null, assuming a normal distribution approximation

(2) make plots with histogram showing distribution of enrichment scores of real vs random_mean, and random_i vs random_mean
overlaid with boxplots showing distribution scores of select kmers of interest
"""

import argparse
import pandas as pd
import math
import os
import numpy as np
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

from itertools import product

# IUPAC ambiguity codes
IUPAC_BASES = {
    'A': ['A'],
    'T': ['T'],
    'C': ['C'],
    'G': ['G'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'M': ['A', 'C'],
    'K': ['G', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}
VALID_BASES = ['A', 'T', 'C', 'G']
BACTERIAL_KMER_SUBSET = ['TTGACA', 'TATAAT', 'AGGAGG']
ARCAHEA_KMER_SUBSET = ['TTTATA', 'TTATAT', 'TATATA']
DAM_DCM_MOTIFS = ['GATC', 'CCAGG', 'CCTGG']

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, base) for base in reversed(seq.upper()))

def test_kmer_enrichment(kmer_df, rs_sites, dam_dcm_motifs):
    """
    given dataframe of kmers and their frequencies 
    column 1 (kmer) sequence
    column 2 (original) true genomic frequencies of kmer 
    columns 3-n (random_i) frequencies of kmer in the randomized genomes
    """
    result_df = pd.DataFrame(kmer_df['kmer'])

    mean_random = kmer_df.iloc[:, 2:].mean(axis=1)
    std_random = kmer_df.iloc[:, 2:].std(axis=1)
    result_df['es'] = (kmer_df['original'] / mean_random).apply(lambda x: math.log2(x))

    # Calculate Z-scores and p-values for kmer frequencies
    z_scores = (kmer_df['original'] - mean_random) / std_random
    p_values = 2 * norm.sf(abs(z_scores))  # sf = 1 - cdf

    result_df['z_score'] = z_scores
    result_df['p_value'] = p_values

    # p-value correction using FDR method (Benjamini-Hochberg)
    # Drop NaNs and correct
    valid = result_df['p_value'].notna()
    rejected, corrected_pvals, _, _ = multipletests(result_df.loc[valid, 'p_value'], method='fdr_bh')

    # Assign back
    result_df.loc[valid, 'p_value_adj'] = corrected_pvals
    result_df.loc[valid, 'significant'] = rejected
    result_df = result_df.sort_values('p_value_adj', ascending=True)

    # annotate restriction modification sites
    result_df['is_rsite'] = result_df['kmer'].isin(with_reverse_complements(rs_sites))
    result_df['is_mt_motif'] = result_df['kmer'].apply(lambda x: contains_any_motif(x, dam_dcm_motifs))

    return result_df

def with_reverse_complements(seq_list):
    """Return a set containing the original sequences and their reverse complements."""
    rc_list = [reverse_complement(seq) for seq in seq_list]
    return list(set(seq_list) | set(rc_list))

def contains_any_motif(kmer, motif_list):
    """Return True if any motif or its reverse complement is a substring of the kmer."""
    all_motifs = with_reverse_complements(motif_list)
    return any(m in kmer for m in all_motifs)

def find_kmers_containing_motif(kmer_list, motifs):
    """Return all kmers from kmer_list that contain any of the motifs or their RCs."""
    motifs_rc = with_reverse_complements(motifs)
    return [k for k in kmer_list if any(m in k for m in motifs_rc)]

def plot_kmer_histogram(kmer_df, kmers_of_interest, plot_filename, monomeric_repeats=None, 
                        restriction_sites=None, dam_dcm_motifs=None):
    """
    kmer_df: Pandas df output of the mcmc/cds_randomizer.py pipeline which should contain kmer frequencies for real/randomized genomes
    kmers_of_interest: list of kmers of interest (TATA-like, etc)
    plot_filename: output name for the plot
    monomeric_repeats, restriction_sites, dam_dcm_motifs: optional args to specify additional categories of kmers of interest

    generates histogram + boxplot figure 
    histogram: distribution of real vs random_mean and random_i vs random_mean enrichment scores
    * idea is to visualize the variance of true kmer distribution from random frequencies 
    * relative to the overall noisiness of the genome
    boxplots: distribution of real vs random_i enrichment scores for kmers of interest
    """
    # Compute mean randomized counts per kmer
    # assumes that the kmer_df has column 0 (kmer sequence), 1 (original) 2- (random)
    n_random_mean = kmer_df.iloc[:, 2:].mean(axis=1)

    # Compute enrichment scores for real vs each random genome
    enrichment_real_vs_random = pd.DataFrame(kmer_df['kmer'])
    for i in range(2, kmer_df.shape[1]):
        enrichment_real_vs_random[f'random_{i-1}'] = np.log2(kmer_df['original'] / kmer_df.iloc[:, i])  # shape: (kmers, n_random_genomes)

    # Compute enrichment scores for random mean vs each random genome
    enrichment_mean_vs_random = pd.DataFrame(kmer_df['kmer'])
    for i in range(2, kmer_df.shape[1]):
        enrichment_mean_vs_random[f'random_{i-1}'] = np.log2(n_random_mean / kmer_df.iloc[:, i])  # shape: (kmers, n_random_genomes)

    # Flatten for histogram plotting
    real_vs_random_vals = enrichment_real_vs_random.iloc[:, 1:].values.flatten()
    mean_vs_random_vals = enrichment_mean_vs_random.iloc[:, 1:].values.flatten()

    # organize data into format suitable for boxplots
    boxplot_data = []
    boxplot_labels = []

    # Collect kmer groups with RCs in deterministic order
    seen = set()
    grouped_kmers = []
    for kmer in kmers_of_interest:
        rc = reverse_complement(kmer)
        canonical = min(kmer, rc)
        pair = (canonical, max(kmer, rc))
        if pair not in seen:
            seen.add(pair)
            grouped_kmers.append(pair)

    # Sort the grouped kmers alphabetically
    grouped_kmers.sort()
    for kmer1, kmer2 in grouped_kmers:
        group = [kmer1, kmer2]
        subset = enrichment_real_vs_random[enrichment_real_vs_random['kmer'].isin(group)]
        if not subset.empty:
            boxplot_data.append(subset.iloc[:, 1:].values.flatten())
            boxplot_labels.append(f"{kmer1}/{kmer2}")

    # DAM/DCM motif-matching kmers
    if dam_dcm_motifs != None:
        methylation_kmers = find_kmers_containing_motif(
            enrichment_real_vs_random['kmer'].tolist(), dam_dcm_motifs
        )
        subset_methyl = enrichment_real_vs_random[enrichment_real_vs_random['kmer'].isin(methylation_kmers)]
        if not subset_methyl.empty:
            boxplot_data.append(subset_methyl.iloc[:, 1:].values.flatten())
            boxplot_labels.append("DAM/DCM motifs")

    # Restriction sites (with RCs) as one aggregate boxplot
    if restriction_sites != None:
        restriction_all = with_reverse_complements(restriction_sites)
        subset_restrict = enrichment_real_vs_random[enrichment_real_vs_random['kmer'].isin(restriction_all)]
        if not subset_restrict.empty:
            boxplot_data.append(subset_restrict.iloc[:, 1:].values.flatten())
            boxplot_labels.append('Restriction sites')

    # Monomeric repeats (with RCs) as one aggregate boxplot
    if monomeric_repeats != None:
        repeats_all = monomeric_repeats
        subset_repeats = enrichment_real_vs_random[enrichment_real_vs_random['kmer'].isin(repeats_all)]
        if not subset_repeats.empty:
            boxplot_data.append(subset_repeats.iloc[:, 1:].values.flatten())
            boxplot_labels.append('Monomeric repeats')

    # Generate boxplots + histogram
    plt.rcParams['font.family'] = 'Helvetica'
    fig = plt.figure(figsize=(8, 6))
    gs = fig.add_gridspec(2, height_ratios=[len(boxplot_labels)*0.5, 3], hspace=0.05)

    # Boxplots (top)
    ax_box = fig.add_subplot(gs[0, 0])
    bp = ax_box.boxplot(boxplot_data, vert=False, patch_artist=True,
                        labels=boxplot_labels, showfliers=False, boxprops=dict(facecolor="#E99A55", color='black'), 
                        medianprops=dict(color="#D35815", linewidth=1.5))

    ax_box.set_xlabel('')
    ax_box.set_yticklabels(boxplot_labels, fontsize=12)
    ax_box.invert_yaxis()
    #ax_box.grid(axis='x')
    ax_box.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

    for whisker in bp['whiskers']:
        whisker.set_visible(False)
    for cap in bp['caps']:
        cap.set_visible(False)

    species_name = plot_filename.split('_kmer_')[0].split('_')[0][0] + ". " + plot_filename.split('_kmer_')[0].split('_')[1]

    # Histogram (bottom)
    ax_hist = fig.add_subplot(gs[1, 0], sharex=ax_box)
    ax_hist.hist(real_vs_random_vals, bins=300, color="#5473CF", label=species_name, alpha=0.6)
    ax_hist.hist(mean_vs_random_vals, bins=300, alpha=0.6, label='Random', color="#545E52")
    ax_hist.set_xlabel(r'$\log_{2}$' + '(real / randomized)', fontsize=14, fontweight='bold')
    ax_hist.set_ylabel('Frequency', fontsize=14, fontweight='bold')
    ax_hist.tick_params(axis='x', labelsize=12)
    ax_hist.tick_params(axis='y', labelsize=12)
    ax_hist.grid(False)
    ax_hist.set_xlim(-1.25, 1.25)
    #ax_hist.grid(axis='y', visible=False)
    ax_hist.legend(fontsize=12)

    # add a dotted vertical line at x=0
    ax_box.axvline(0, color="#7A7E7A", linestyle="--", linewidth=2)
    ax_hist.axvline(0, color="#7A7E7A", linestyle="--", linewidth=2)

    plt.savefig(plot_filename, dpi=1200, bbox_inches='tight')
   
def get_restriction_sites(rs_filename, k, org_filter="Escherichia coli"):
    """
    given a list of restriction sites in REBASE format, collect a list of restriction enzyme sites
    of length `k` which correspond to the species given by `org_filter`
    """
    def expand_ambiguous_sequence(seq):
        """Expand IUPAC ambiguous sequence into all unambiguous combinations."""
        try:
            bases = [IUPAC_BASES[base] for base in seq]
        except KeyError as e:
            print(f"Warning: Invalid base {e} in sequence {seq}")
            return []
        return [''.join(p) for p in product(*bases)]
    
    with open(rs_filename, 'r') as f:
        lines = f.readlines()

    valid_sites = set()
    current_organism = ""

    for i, line in enumerate(lines):
        if line.startswith('<3>'):
            current_organism = line.strip()[3:].strip()

        elif line.startswith('<5>') and org_filter in current_organism:
            raw_seq = line.strip()[3:].strip()

            # Trim anything from '(' onward
            raw_seq = raw_seq.split('(')[0].strip()
            if not raw_seq:
                continue

            clean_seq = raw_seq.replace('^', '')
            if len(clean_seq) == k and all(base in IUPAC_BASES for base in clean_seq):
                expanded = expand_ambiguous_sequence(clean_seq)
                valid_sites.update(expanded)

    return list(valid_sites)

def main():
    parser = argparse.ArgumentParser(description="Extract N-bp restriction sites from a REBASE file and match against k-mer frequency data.")

    parser.add_argument('kmer_file', type=str, help='Path to kmer frequency csv file')
    parser.add_argument('-k', '--length', type=int, required=True, help='Length of kmers (e.g., 6)')
    parser.add_argument('-d',  '--domain', type=str, required=True, choices=['bacteria', 'archaea'], help='domain: either bacteria or archaea')
    parser.add_argument('-r', '--rs_filename', type=str, required=False, help='txt file containing restriction enzyme target sites')

    args = parser.parse_args()

    # Extract species name prefix
    base_name = os.path.basename(args.kmer_file)
    species_name = '_'.join(base_name.split('_kmer_')[0].split('_')[:2])
    k = args.length

    stats_output = f"{species_name}_kmer_stats_k{k}.csv"
    plot_output = f"{species_name}_kmer_plot_k{k}.png"

    # Load k-mer frequencies
    kmer_df = pd.read_csv(args.kmer_file)
    kmer_df.replace(0, 1, inplace=True)  # avoid divisions by 0

    # Restriction enzyme targets
    species_name_formatted = species_name.replace('_', ' ')
    rs_sites = get_restriction_sites(args.rs_filename, args.length, org_filter=species_name_formatted)
    print(len(rs_sites), "restriction enzyme target sites")

    # Calculate enrichment
    result_df = test_kmer_enrichment(kmer_df, rs_sites, DAM_DCM_MOTIFS)
    result_df.to_csv(stats_output, index=False)

    monomer_repeats = [args.length * base for base in VALID_BASES]

    # Plotting
    if args.domain.lower() == "archaea":
        plot_kmer_histogram(kmer_df, ARCAHEA_KMER_SUBSET, plot_output, monomer_repeats)
    else:
        plot_kmer_histogram(kmer_df, BACTERIAL_KMER_SUBSET, plot_output, monomer_repeats, restriction_sites=rs_sites)

if __name__ == "__main__":
    main()
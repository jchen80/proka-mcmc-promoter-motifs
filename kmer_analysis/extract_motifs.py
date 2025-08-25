"""
Author: Julie Chen
extract_motifs.py

Gets the enrichment scores of a variety of promoter motifs (as defined by sliding window over consensus regexes
or by best learned consensus motif from tRNA promoter scan)
"""

import argparse
import pandas as pd
import numpy as np
import os
import itertools

ARCHAEAL_FIXED = ['TTTATATA'] 
ARCHAEAL_REGEX = ["TTTTAAA", "TTTWWW", "TTTAWATA"] 

BACTERIA_FIXED = ['TTGACA', 'TATAAT']

def expand_w(seq):
    """Expand ambiguous 'W' bases in a motif to all possible combinations (A/T)."""
    positions = [( ['A','T'] if ch == 'W' else [ch]) for ch in seq]
    return [''.join(p) for p in itertools.product(*positions)]

def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    complement = str.maketrans("ATGC", "TACG")
    return seq.translate(complement)[::-1]

def extract_mean_enrichment_score(kmer_file, kmer, include_rc=True):
    """Compute mean log2 enrichment score for a kmer in a species kmer file."""
    df = pd.read_csv(kmer_file)
    df.replace(0, 1, inplace=True) # avoid division by 0

    random_cols = [col for col in df.columns if col.startswith("random_")]
    subset = df[df['kmer'] == kmer]

    if subset.empty:
        return np.nan

    original = subset['original'].values[0]
    enrichment = np.mean(np.log2(original / subset[random_cols].values.flatten()))

    if include_rc:
        rc_subset = df[df['kmer'] == reverse_complement(kmer)]
        if not rc_subset.empty:
            rc_original = rc_subset['original'].values[0]
            rc_enrichment = np.mean(np.log2(rc_original / rc_subset[random_cols].values.flatten()))
            enrichment = np.mean([enrichment, rc_enrichment])

    return enrichment

def sliding_window(seq, k):
    """Return all length-k substrings of seq."""
    return [seq[i:i+k] for i in range(len(seq)-k+1)]

def is_monorepeat(seq):
        """Return True if all characters in seq are the same."""
        return len(set(seq)) == 1

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('parent_dir', help='Parent directory with species subdirs')
    parser.add_argument('learned_kmers', help='CSV file with learned tRNA scan kmers')
    parser.add_argument('--nap_file', type=str, required=True, help='CSV with species and NAP metadata')
    parser.add_argument('-d', '--domain', type=str, required=True, choices=['bacteria', 'archaea'])
    parser.add_argument('--k', type=int, required=True, help='kmer length for sliding window motifs')
    parser.add_argument('--include_rc', action='store_true', default=True, help='Include reverse complement kmer')
    args = parser.parse_args()

    # Load tRNA scan motifs
    learned_df = pd.read_csv(args.learned_kmers)
    learned_df.rename(columns={'kmer': 'tRNA_scan_motif'}, inplace=True)

    # Prepare dynamic kmer lists, i.e. kmers of interest
    if args.domain == "archaea":
        fixed_kmers = list(itertools.chain.from_iterable([sliding_window(seq, args.k) for seq in ARCHAEAL_FIXED]))
        regex_windows = list(itertools.chain.from_iterable([sliding_window(seq, args.k) for seq in ARCHAEAL_REGEX]))
        regex_kmers = list(itertools.chain.from_iterable([expand_w(kmer) for kmer in regex_windows]))

        # Filter fixed kmers
        fixed_kmers = [kmer for kmer in fixed_kmers if not is_monorepeat(kmer)]

        # Filter regex kmers
        regex_kmers = [kmer for kmer in regex_kmers if not is_monorepeat(kmer)]
        regex_kmers = [kmer for kmer in regex_kmers if kmer not in fixed_kmers]
        all_dynamic_kmers = fixed_kmers + regex_kmers
    
    else:
        fixed_kmers = list(itertools.chain.from_iterable([sliding_window(seq, args.k) for seq in BACTERIA_FIXED]))
        fixed_kmers = [kmer for kmer in fixed_kmers if not is_monorepeat(kmer)]
        all_dynamic_kmers = fixed_kmers

    results = []
    for species in os.listdir(args.parent_dir):
        species_dir = os.path.join(args.parent_dir, species)
        if not os.path.isdir(species_dir):
            continue

        kmer_file = os.path.join(species_dir, f"{species}_kmer_frequencies_k{args.k}.csv")
        if not os.path.exists(kmer_file):
            continue

        row = {'species': species}

        # Compute enrichment for tRNA motif
        species_motif_row = learned_df[learned_df['species'] == species]
        if not species_motif_row.empty:
            tRNA_motif = species_motif_row['tRNA_scan_motif'].values[0]
            row['tRNA_scan_sequence'] = tRNA_motif
            row['tRNA_scan_motif'] = extract_mean_enrichment_score(kmer_file, tRNA_motif, args.include_rc)
        else:
            row['tRNA_scan_sequence'] = np.nan
            row['tRNA_scan_motif'] = np.nan

        # Compute enrichment for each dynamic kmer
        for kmer in all_dynamic_kmers:
            row[kmer] = extract_mean_enrichment_score(kmer_file, kmer, args.include_rc)

        if args.domain == "archaea":
            # fixed_kmers
            fixed_vals = {kmer: row[kmer] for kmer in fixed_kmers if not np.isnan(row[kmer])}
            if fixed_vals:
                min_fixed_kmer = min(fixed_vals, key=fixed_vals.get)  # kmer with smallest value
                row['ARCHAEAL_FIXED_min'] = fixed_vals[min_fixed_kmer]
                row['ARCHAEAL_FIXED_min_kmer'] = min_fixed_kmer
            else:
                row['ARCHAEAL_FIXED_min'] = np.nan
                row['ARCHAEAL_FIXED_min_kmer'] = np.nan

            # regex_kmers
            regex_vals = {kmer: row[kmer] for kmer in regex_kmers if not np.isnan(row[kmer])}
            if regex_vals:
                min_regex_kmer = min(regex_vals, key=regex_vals.get)
                row['ARCHAEAL_REGEX_min'] = regex_vals[min_regex_kmer]
                row['ARCHAEAL_REGEX_min_kmer'] = min_regex_kmer
            else:
                row['ARCHAEAL_REGEX_min'] = np.nan
                row['ARCHAEAL_REGEX_min_kmer'] = np.nan

        # fixed (+ regex) kmers
            tata_vals = {kmer: row[kmer] for kmer in all_dynamic_kmers if not np.isnan(row[kmer])}
            if tata_vals:
                min_tata_kmer = min(tata_vals, key=tata_vals.get)
                row['TATA_min'] = tata_vals[min_tata_kmer]
                row['TATA_min_kmer'] = min_tata_kmer
            else:
                row['TATA_min'] = np.nan
                row['TATA_min_kmer'] = np.nan
            
            # Compute top 3 most negative enrichments across TATA motifs
            row['TATA_avg_top3'] = np.mean(sorted(tata_vals.values())[:3]) if tata_vals else np.nan

        # overall min across all dynamic kmers + tRNA motif
        all_vals = {kmer: row[kmer] for kmer in all_dynamic_kmers if not np.isnan(row[kmer])}
        if not np.isnan(row['tRNA_scan_motif']):
            all_vals[row['tRNA_scan_sequence']] = row['tRNA_scan_motif']

        if all_vals:
            min_overall_kmer = min(all_vals, key=all_vals.get)
            row['overall_min'] = all_vals[min_overall_kmer]
            row['overall_min_kmer'] = min_overall_kmer
        else:
            row['overall_min'] = np.nan
            row['overall_min_kmer'] = np.nan

        results.append(row)

    # Merge with nap_file metadata
    nap_df = pd.read_csv(args.nap_file)
    result_df = pd.DataFrame(results)
    merged = nap_df.merge(result_df, on='species', how='left')

    merged.to_csv(f"motif_enrichment_summary_{args.domain}_k{args.k}.csv", index=False)
    print("Done! Output saved to:", f"motif_enrichment_summary_{args.domain}_k{args.k}.csv")

if __name__ == '__main__':
    main()


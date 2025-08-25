"""
Author: Julie Chen
plot_kmer_compare_boxplot.py

Given an input csv file that contains URL links to base directory for download, 
download the cds_from_genomic, genomic, and translated_cds files
"""

import os
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def reverse_complement(seq):
    complement = str.maketrans("ATGC", "TACG")
    return seq.translate(complement)[::-1]

def shorten_species_name(name):
    """Shorten species name for plotting, e.g. Archaeoglobus_fulgidus -> A. fulgidus"""
    parts = name.split('_')
    if len(parts) >= 2:
        return f"{parts[0][0]}. {parts[1]}"
    else:
        return name

def extract_enrichment_scores(kmer_file, kmer, include_rc=False):
    """
    gets enrichment scores (log2(real/randomized)) for a specific kmer from a CSV file 
    containing kmer frequencies in real vs randomized genomes (kmer_file)
    """
    df = pd.read_csv(kmer_file)
    df.replace(0, 1, inplace=True)
    random_cols = [col for col in df.columns if col.startswith("random_")]
    subset = df[df['kmer'].isin([kmer])]

    if subset.empty:
        raise ValueError(f"Kmer '{kmer}' not found in file '{kmer_file}'")
    
    original = subset['original'].values[0]
    enrichment = np.log2(original / subset[random_cols].values.flatten())
    
    if include_rc:
        rc_subset = df[df['kmer'].isin([reverse_complement(kmer)])]
        if rc_subset.empty:
            raise ValueError(f"Reverse complement of {kmer} not found in {kmer_file}")
        rc_original = rc_subset['original'].values[0]
        rc_enrichment = np.log2(rc_original / rc_subset[random_cols].values.flatten())
        enrichment = np.concatenate((enrichment, rc_enrichment))

    return enrichment

def generate_boxplots(merged, kmer_name, compare_name=None):
    merged['species_reformatted'] = merged["species"].apply(shorten_species_name)
    # Order species by median enrichment of the main kmer
    species_order = (
        merged[merged['source'] == kmer_name]
        .groupby('species_reformatted')['enrichment']
        .median()
        .sort_values()
        .index
        .tolist()
    )
    merged['species_reformatted'] = pd.Categorical(merged['species_reformatted'], categories=species_order, ordered=True)

    # Figure
    plt.figure(figsize=(8, max(6, len(species_order) * 0.6)))

    palette = ["#F2923E", "#7D9AF1"]

    ax = sns.boxplot(
        data=merged, y='species_reformatted', x='enrichment',
        hue='source', showfliers=False, whis=0, linewidth=1, palette=palette, saturation=1
    )
    ax.axvline(0, color='lightgrey', linestyle='--', linewidth=2, zorder=1)
    ax.set_ylabel("")
    ax.tick_params(axis='x', labelsize=14, labelfontfamily="Helvetica")
    ax.tick_params(axis='y', labelsize=14, labelfontfamily="Helvetica")
    ax.set_xlabel(r'$\log_{2}$' + '(real / randomized)', fontsize=14, fontname="Helvetica")
    #ax.set_title(f"Kmer enrichment for {kmer_name}")
    plt.legend()

    if compare_name:
        plt.savefig(f"{kmer_name}_vs_{compare_name}.png", dpi=1200, bbox_inches="tight")
    else:
        plt.savefig(f"{kmer_name}.png", dpi=1200, bbox_inches="tight")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Plot enrichment of kmers across archaeal species.")
    parser.add_argument("parent_dir", help="Parent directory containing species subfolders with kmer CSVs")
    parser.add_argument("kmer", help="Kmer of interest, or CSV of motifs per species")
    parser.add_argument("--compare", required=True, help="A single kmer sequence to compare against")
    parser.add_argument("--include_rc", default=True, action="store_true")

    args = parser.parse_args()

    kmers_as_file = args.kmer.endswith('.csv')
    if kmers_as_file:
        query_kmers = pd.read_csv(args.kmer)
        k = len(query_kmers['kmer'].values[0])
    else:
        k = len(args.kmer)

    kmer_name = args.kmer if not kmers_as_file else f"tRNA promoter motifs k={k}"

    data = []
    for species in os.listdir(args.parent_dir):
        species_dir = os.path.join(args.parent_dir, species)
        if not os.path.isdir(species_dir):
            continue

        kmer_file = os.path.join(species_dir, f"{species}_kmer_frequencies_k{k}.csv")
        if not os.path.exists(kmer_file):
            continue

        if kmers_as_file:
            if species not in query_kmers['species'].values:
                continue
            query_kmer = query_kmers.loc[query_kmers['species'] == species, 'kmer'].values[0]
        else:
            query_kmer = args.kmer

        enrichment = extract_enrichment_scores(kmer_file, query_kmer, args.include_rc)
        for score in enrichment:
            data.append({'species': species, 'enrichment': score, 'source': kmer_name})

        compare_enrichment = extract_enrichment_scores(kmer_file, args.compare, args.include_rc)
        for score in compare_enrichment:
            data.append({'species': species, 'enrichment': score, 'source': args.compare})
    
    merged = pd.DataFrame(data)

    generate_boxplots(merged, kmer_name, compare_name=args.compare)

if __name__ == "__main__":
    main()

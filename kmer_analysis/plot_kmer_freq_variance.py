"""
Author: Julie Chen
plot_kmer_freq_variance.py

plots the mean vs normalized variance of kmer frequencies across the multiple randomized genomes 
key idea is to address the question: are kmers which have low overall frequency also more variable across the ranodomizations?
"""

import pandas as pd
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Plot mean vs variance of kmer frequencies across randomized genomes.")
    parser.add_argument('--parent_dir', type=str, required=True, help='Path to top-level directory containing species subfolders')
    parser.add_argument('--kmer_length', type=int, required=True, help='K-mer length (e.g., 6)')
    args = parser.parse_args()

    parent = Path(args.parent_dir)
    pattern = f"*kmer_frequencies_k{args.kmer_length}.csv"
    
    all_kmer_dfs = []

    for filepath in parent.rglob(pattern):
        df = pd.read_csv(filepath)
        all_kmer_dfs.append(df)

    if not all_kmer_dfs:
        print("No matching kmer frequency files found.")
        return

    combined_df = pd.concat(all_kmer_dfs, ignore_index=True)

    # Extract columns corresponding to randomized frequencies
    random_cols = [col for col in combined_df.columns if col.startswith('random_')]

    # Compute mean and variance for each k-mer across randomized genomes
    combined_df['random_mean'] = combined_df[random_cols].mean(axis=1)
    combined_df['random_var'] = combined_df[random_cols].var(axis=1)
    
    # Normalize variance by mean² (CV²)
    combined_df = combined_df[combined_df['random_mean'] > 0]  # Avoid divide-by-zero
    combined_df['normalized_var'] = combined_df['random_var'] / (combined_df['random_mean'] ** 2)

    # Plot mean vs normalized variance
    plt.figure(figsize=(6, 5))
    plt.scatter(combined_df['random_mean'], combined_df['normalized_var'], alpha=0.4, s=10)
    plt.xlabel('Mean k-mer frequency across random genomes')
    plt.ylabel('Normalized Variance (CV²)')
    plt.title(f'K={args.kmer_length} mean vs normalized variance across species')
    plt.grid(True)
    plt.tight_layout()
    plt.xscale('log')

    plt.show()

if __name__ == "__main__":
    main()

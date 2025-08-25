"""
Author: Julie Chen
compute_distance_motif.py

Visualization of Levenshtein distance between the found promoter motifs of different species
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import Levenshtein
from tqdm import tqdm
from itertools import product
import sys

def reverse_complement(seq):
    complement = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement)[::-1]

def avg_mean_levenshtein(set1, set2):
    """
    Compute the average of the n smallest Levenshtein distances
    between all motif pairs from set1 and set2.
    """
    distances = sorted(Levenshtein.distance(a, b) for a, b in product(set1, set2))
    return np.mean(distances)

def process_species_motifs(filename):
    df = pd.read_csv(filename).dropna(subset=["kmer"])
    
    # Group motifs by species
    species_to_kmers = {}
    for species, motifs in df.groupby("species")["kmer"]:
        expanded_set = set()
        for motif in motifs:
            expanded_set.add(motif)
            expanded_set.add(reverse_complement(motif))
        species_to_kmers[species] = expanded_set

    species_list = list(species_to_kmers.keys())
    n = len(species_list)
    dist_matrix = np.zeros((n, n))

    for i in tqdm(range(n), desc="Computing species distances"):
        for j in range(i, n):
            d = avg_mean_levenshtein(
                species_to_kmers[species_list[i]],
                species_to_kmers[species_list[j]])
            dist_matrix[i, j] = d
            dist_matrix[j, i] = d

    dist_df = pd.DataFrame(dist_matrix, index=species_list, columns=species_list)
    return dist_df


def main():
    if len(sys.argv) != 2:
        print("Usage: python3 plot_motif.py <kingdom>")
        sys.exit(1)

    kingdom = sys.argv[1]

    # Example usage for k=6
    k6_filename = f"tRNA_promoter/{kingdom}/consensus_combined_top_motifs_k6.csv"
    k7_filename = f"tRNA_promoter/{kingdom}/consensus_combined_top_motifs_k7.csv"

    def generate_distance_plot(filename):
        dist_df = process_species_motifs(filename)

        k = str.split(filename, ".")[0][-1]

        # Save and plot
        sns.clustermap(dist_df, cmap="viridis", figsize=(10, 10))
        plt.title(f"Species Levenshtein Distance Based on Expanded Motif Sets (k={k})")
        plt.savefig(f"tRNA_promoter/{kingdom}/k{k}_promoter_levenshtein_distance.png")
        plt.close()

    generate_distance_plot(k6_filename)
    generate_distance_plot(k7_filename)

if __name__ == "__main__":
    main()



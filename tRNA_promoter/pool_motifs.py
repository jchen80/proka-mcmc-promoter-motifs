"""
Author: Julie Chen
pool_motifs.py 

Loops through species directories in a main parent directory and gathers the 
best consensus sequence per species to collate into a single file
"""

import os
import pandas as pd
import sys

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 pool_motifs.py <base_dir> <output_dir>")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    out_dir = sys.argv[2]
    k_values = [6, 7]

    for k in k_values:
        species_to_motifs = {}

        # loop over all species directories
        for species in os.listdir(base_dir):
            species_dir = os.path.join(base_dir, species)
            if not os.path.isdir(species_dir):
                continue

            fname = f"{species}_consensus_meme_k{k}_scores.csv"
            fpath = os.path.join(species_dir, fname)
            if not os.path.exists(fpath):
                continue

            df = pd.read_csv(fpath)
            df['species'] = species
            species_to_motifs[species] = df

        if not species_to_motifs:
            print(f"No motif score files found for k={k}.")
            continue

        combined_df = pd.concat(species_to_motifs).reset_index(drop=True)
        out_path = os.path.join(out_dir, f"consensus_combined_top_motifs_k{k}.csv")
        combined_df.to_csv(out_path, index=False)
        print(f"Saved balanced combined top k={k} motifs to {out_path}")


if __name__ == "__main__":
    main()
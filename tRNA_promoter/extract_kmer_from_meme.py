"""
Author: Julie Chen
extract_kmer_from_meme.py

This script takes the output from MEME motif searching algorithm and finds the *best* motif of length k={6,7}
that is not monomeric

*best* motif is defined as the motif with the highest prevalence (# of hits found), from which the sliding window
with the highest IC is selected
"""

import os
import math
import pandas as pd
from Bio import motifs
from collections import defaultdict
import sys

def load_meme_motifs(xml_path):
    with open(xml_path) as handle:
        return motifs.parse(handle, "meme")

def parse_maxsites_from_txt(txt_path):
    """input: .txt file output from MEME motif search pipeline
    fetches the total number of query sites for which the motif search attempted to look for matches
    """
    if not os.path.exists(txt_path):
        return None
    with open(txt_path) as f:
        for line in f:
            if line.strip().startswith("nsites:"):
                try:
                    parts = line.strip().split()
                    maxsites = int(parts[4])
                    return maxsites
                except Exception:
                    return None

def compute_ic(position_counts, pseudocount=0.01):
    """
    position counts: dictionary with base name : count mapping, for a single position within the motif
    returns: information content for the position specified by `position_counts` 
    """
    total = sum(position_counts.values()) + 4 * pseudocount
    ic = 2.0 # max entropy value 
    for base in 'ACGT':
        freq = (position_counts.get(base, 0) + pseudocount) / total
        if freq > 0:
            ic += freq * math.log2(freq)
    return ic

def extract_kmers_from_consensus(motif, k=6):
    """given a Bio.motif motif, perform sliding window to extract kmers of length `k`
    returns: dictionary mapping {kmer : (start_index, average IC)} for each kmer along the window
    
    * skips monomeric kmers as possible motifs
    """
    motif_counts = motif.counts
    consensus = str(motif.consensus).upper()

    def avg_ic_for_window(start, k):
        """computes IC for each position within the kmer in the sliding window starting at `start` index
        then averages over the length of window `k`"""
        return sum(compute_ic(motif_counts[:, pos]) for pos in range(start, start + k)) / k

    kmer_to_pos_ic = {}
    # loop through the consensus sequence as sliding window of length k
    for i in range(len(consensus) - k + 1):
        avg_ic = avg_ic_for_window(i, k)
        kmer = consensus[i:i+k]
        # skip if kmer is monomeric
        if len(set(kmer)) == 1:
            continue
        kmer_to_pos_ic[kmer] = (i, avg_ic)
    return kmer_to_pos_ic

def score_kmers_consensus_kmers(kmer_to_pos_ic, instances, total_input_seqs):
    scores = []
    for kmer, (pos, avg_ic) in kmer_to_pos_ic.items():
        prevalence = len(instances) / total_input_seqs
        scores.append({
            'kmer': kmer,
            'motif_pos': pos,
            'avg_ic': avg_ic,
            'prevalence': prevalence})
    return pd.DataFrame(scores)

def process_species_directory(species_dir):
    """
    input: species_dir contains the outputs of the MEME motif search pipeline
    * specifically, there should be files ending in *meme.xml and *meme.txt which contain the results of the MEME algorithm

    species_dir will also be where the output scored kmers are stored
    """
    xml_path = os.path.join(species_dir, "meme.xml")
    txt_path = os.path.join(species_dir, "meme.txt")
    species = os.path.basename(species_dir.rstrip("/"))

    if not os.path.exists(xml_path):
        print(f"[SKIP] No meme.xml in {species_dir}")
        return

    motifs_list = load_meme_motifs(xml_path)
    maxsites = parse_maxsites_from_txt(txt_path) # total number of query sites

    all_kmer_results = []
    for i, motif in enumerate(motifs_list): # loop through each "found" motif
        max_input_seqs = maxsites if maxsites else len(motif.instances)
        print(f"[{species}] Motif {i+1}: width={motif.length}, max_input_seqs={max_input_seqs}")

        # sliding window over motif
        for k in [6, 7]:
            if k > len(motif.consensus):
                continue
            kmer_to_pos_ic = extract_kmers_from_consensus(motif, k=k)
            df = score_kmers_consensus_kmers(kmer_to_pos_ic, motif.instances, max_input_seqs)
            df['motif_index'] = i + 1
            df['k'] = k
            all_kmer_results.append(df)

    if all_kmer_results:
        final_df = pd.concat(all_kmer_results).reset_index(drop=True)

        # separate out kmers of length 6 and 7
        k6 = final_df[final_df['k'] == 6]
        k7 = final_df[final_df['k'] == 7]

        def get_best_kmer(df):
            """get "best" kmer out of possible kmer motifs
            here "best" is defined by (1) first selecting the found motif with the highest prevalence 
            among the total tRNAs, then (2) selecting the sliding window kmer with the highest information content (IC)
            """
            # only keep rows with highest prevalence
            max_prevalence = df['prevalence'].max()
            df_max_prevalence = df[df['prevalence'] == max_prevalence]

            # Of those, keep the row where avg_ic is the greatest
            result = pd.DataFrame(df_max_prevalence.loc[[df_max_prevalence['avg_ic'].idxmax()]])
            return result
        
        k6 = get_best_kmer(k6)
        k7 = get_best_kmer(k7)

        k6.to_csv(os.path.join(species_dir, f"{species}_consensus_meme_k6_scores.csv"), index=False)
        k7.to_csv(os.path.join(species_dir, f"{species}_consensus_meme_k7_scores.csv"), index=False)

def main():

    if len(sys.argv) != 2:
        print("Usage: python3 extract_kmer_from_meme.py <base_dir>")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    if not os.path.exists(base_dir):
        print(f"Base directory {base_dir} does not exist.")
        sys.exit(1)

    for species_name in os.listdir(base_dir):
        species_dir = os.path.join(base_dir, species_name)
        if os.path.isdir(species_dir):
            process_species_directory(species_dir)

if __name__ == "__main__":
    main()

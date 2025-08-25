"""
Author: Julie Chen
cds_randomizer.py

Main workhorse script that performs the Monte Carlo Markov chain randomization swapping codons 
"""

import math
import argparse
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt

import re
import copy

import numpy as np
from collections import defaultdict

from numba import njit, types
from numba.typed import Dict

class CdsRandomizer():
    codon_to_int = {}
    int_to_codon = {}

    def __init__(self):
        self.build_codon_maps()
        self.cds_entries = []  # Each entry: dict with keys: id, codons (str), aa (str)
        self.aa_to_int = {}
        self.int_to_aa = {}

    @staticmethod
    def build_codon_maps():
        """
        Builds codon to integer and integer to codon mappings.
        This is a static method that initializes the codon mappings used for encoding and decoding.
        """
        bases = ['T', 'C', 'A', 'G']
        codons = [a + b + c for a in bases for b in bases for c in bases]
        CdsRandomizer.codon_to_int = {codon: i for i, codon in enumerate(codons)}
        CdsRandomizer.int_to_codon = {i: codon for i, codon in enumerate(codons)}

    def build_aa_map(self):
        """
        Builds amino acid to integer and integer to amino acid mappings.
        """
        aa_set = set(''.join(entry['aa'] for entry in self.cds_entries))
        self.aa_to_int = {aa: i for i, aa in enumerate(sorted(aa_set))}
        self.int_to_aa = {i: aa for aa, i in self.aa_to_int.items()}

    def add_entry(self, refseq_id, seq, aa):
        """ Adds a CDS entry with the given reference sequence ID, nucleotide sequence, and amino acid sequence."""
        codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
        self.cds_entries.append({'id': refseq_id, 'codons': codons, 'aa': aa})

    def read_cds_entries(self, seq_file, aa_file):
        """ Reads CDS entries from the given files and populates the cds_entries list.
        The seq file should contain nucleotide sequences, and the translated file should contain amino acid sequences.
        Each entry is identified by a unique RefSeq ID, which is extracted from the header lines of the files, allowing for
        the association of nucleotide and amino acid sequences.

        Currently, if coding segments are annotated as pseudogenes `pseudo=true` flag present, then we skip those entries
        Adds a stop codon to the translated amino acid sequence
        """

        def validate_and_fix_aa_seq(refseq_id, nt_seq, aa_seq):
            """
            Validates and adjusts sequences:
            - If nt_seq is divisible by 3, only accept if aa_seq is exactly one short (missing stop codon).
            - If nt_seq is not divisible by 3, trim to full codons, then ensure codon count matches aa length.
            """
            if len(nt_seq) % 3 != 0:
                # Case: Needs trimming
                remainder = len(nt_seq) % 3
                nt_seq = nt_seq[:-remainder]
                codon_len = len(nt_seq) // 3

                if len(aa_seq) > codon_len:
                    aa_seq = aa_seq[:codon_len]  # Trim amino acids to match
                    return aa_seq, nt_seq
                elif len(aa_seq) == codon_len:
                    return aa_seq, nt_seq
                else:
                    raise ValueError(
                        f"[TRIMMED] Amino acid length too short for {refseq_id}: "
                        f"{len(aa_seq)} amino acids vs {codon_len} codons"
                    )
            else:
                # Case: Perfectly divisible
                codon_len = len(nt_seq) // 3
                if len(aa_seq) == codon_len - 1:
                    return aa_seq + '*', nt_seq 
                elif len(aa_seq) == codon_len:
                    return aa_seq, nt_seq 
                else:
                    raise ValueError(
                        f"[DIV3] Expected {codon_len - 1} or {codon_len} amino acids for {refseq_id}, got {len(aa_seq)}"
                    )
                
        def read_refseq_cds_translated(cds_file, translated_file):
            sequences = {}  # dictionary to hold both nucleotide and amino acid sequences with RefSeq ID as key

            # Read CDS sequences
            current_id = None
            current_sequence = ""
            skip_entry = False
            with open(cds_file, 'r') as cds_fh:
                for line in cds_fh:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_id is not None and not skip_entry:
                            sequences[current_id] = {'seq': current_sequence, 'aa': ''}
                        # Check for pseudogene flag
                        skip_entry = '[pseudo=true]' in line.lower()
                        current_id = line[1:].split('|')[1].split(' ')[0].split('_cds')
                        current_id = "".join(current_id)
                        current_sequence = ""
                    else:
                        if not skip_entry:
                            current_sequence += line
                if current_id is not None and not skip_entry:
                    sequences[current_id] = {'seq': current_sequence, 'aa': ''}

            # Read translated amino acid sequences
            with open(translated_file, 'r') as translated_fh:
                current_id = None
                current_aa_sequence = ""
                for line in translated_fh:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_id is not None and current_id in sequences:
                            nt_seq = sequences[current_id]['seq']
                            aa_seq, nt_seq = validate_and_fix_aa_seq(current_id, nt_seq, current_aa_sequence)
                            sequences[current_id]['aa'] = aa_seq
                        current_id = line[1:].split('|')[1].split(' ')[0].split('_prot')
                        current_id = "".join(current_id)
                        current_aa_sequence = ""
                    else:
                        current_aa_sequence += line

                # Handle final entry
                if current_id is not None and current_id in sequences:
                    nt_seq = sequences[current_id]['seq']
                    aa_seq, nt_seq = validate_and_fix_aa_seq(current_id, nt_seq, current_aa_sequence)
                    sequences[current_id]['aa'] = aa_seq

            return sequences

        sequences = read_refseq_cds_translated(seq_file, aa_file)
        #print(len(sequences), "sequences read from files")
        for refseq_id, data in sequences.items():
            self.add_entry(refseq_id, data['seq'], data['aa'])

        self.build_aa_map()

    def write_cds_entries(self, out_seq_filename, aa_filename, entries=None):
        """
        Write out the nucleotide and amino acid sequences from the randomized entries to fasta files
        each entry in the output file should consist of a header line of the form `>{entry.id}`
        followed by either the nucleotide (entry.codons decoded) or amino acid sequence (entry.aa) in the next line
        """
        if entries == None:
            entries = self.cds_entries
        with open(out_seq_filename, 'w') as nt_file, open(aa_filename, 'w') as aa_file:
            for entry in entries:
                # Build nucleotide and amino acid sequences
                nt_seq = ''.join(entry.codons)
                aa_seq = entry.aa

                # Write nucleotide FASTA
                nt_file.write(f">{entry.id}\n")
                nt_file.write(f"{nt_seq}\n")

                # Write amino acid FASTA
                aa_file.write(f">{entry.id}\n")
                aa_file.write(f"{aa_seq}\n")

    def estimate_max_positions_per_context(self, codon_arr, aa_arr):
        counts = defaultdict(int)
        n_codons = len(self.codon_to_int)
        n_aas = len(self.aa_to_int)

        for i in range(1, len(codon_arr) - 1):
            prev = codon_arr[i - 1]
            aa = aa_arr[i]
            next_ = codon_arr[i + 1]
            ctx_id = prev * n_aas * n_codons + aa * n_codons + next_
            counts[ctx_id] += 1
        max_count = max(counts.values()) if counts else 0
        return max_count

    def flatten_with_tracking(self):
        codon_ints = []
        aa_ints = []
        index_map = []

        entry_lengths = {}

        # First pass: collect codons, a.a.s, and index mapping
        for entry_index, entry in enumerate(self.cds_entries):
            for codon_index, (codon, aa) in enumerate(zip(entry['codons'], entry['aa'])):
                if codon not in self.codon_to_int or aa not in self.aa_to_int:
                    continue
                codon_int = self.codon_to_int[codon]
                aa_int = self.aa_to_int[aa]
                codon_ints.append(codon_int)
                aa_ints.append(aa_int)
                index_map.append((entry_index, codon_index))
                entry_lengths.setdefault(entry_index, 0)
                entry_lengths[entry_index] = max(entry_lengths[entry_index], codon_index + 1)
        
        # Second pass: compute valid swap positions which are not at the edges of the codon sequence
        valid_swap_pos = []
        for (entry_idx, codon_idx) in index_map:
            is_edge = codon_idx == 0 or codon_idx == (entry_lengths[entry_idx] - 1)
            valid_swap_pos.append(not is_edge)

        #print(len(codon_ints), "codons found in total")
        #print(len(aa_ints), "amino acids found in total")
        return (
            np.array(codon_ints, dtype=np.int32),
            np.array(aa_ints, dtype=np.int32),
            np.array(index_map, dtype=np.int32),
            np.array(valid_swap_pos, dtype=bool)
        )

    def run_mcmc(self, seeded=False, max_swaps=None, checkpoints=None):
        codon_arr, aa_arr, index_map, valid_swap_pos = self.flatten_with_tracking()

        n_codons = 64
        n_aas = len(self.aa_to_int)
        max_count = self.estimate_max_positions_per_context(codon_arr, aa_arr)
        max_pos_per_context = max(500, int(max_count * 3))

        n = 3 * len(codon_arr)
        if max_swaps is None:
            max_swaps = int(3 * n * math.log(n))

        if checkpoints is not None:
            checkpoints = np.array(checkpoints, dtype=np.int32)
            results = _mcmc_codons_with_checkpoints(
                codon_arr, aa_arr, n_codons, n_aas, checkpoints, max_pos_per_context, valid_swap_pos, seeded
            )
            randomized_checkpoints = {}
            for swap_count, codon_arr_result in results.items():
                entries = copy.deepcopy(self.cds_entries)
                for flat_index, codon_int in enumerate(codon_arr_result):
                    entry_i, codon_i = index_map[flat_index]
                    codon_str = self.int_to_codon[codon_int]
                    entries[entry_i]['codons'][codon_i] = codon_str
                randomized_checkpoints[swap_count] = entries
            return randomized_checkpoints

        else:
            n_swaps = max_swaps
            new_codon_arr = _mcmc_codons(
                codon_arr, aa_arr, n_codons, n_aas, n_swaps, max_pos_per_context, valid_swap_pos, seeded
            )
            randomized_entries = copy.deepcopy(self.cds_entries)
            for flat_index, codon_int in enumerate(new_codon_arr):
                entry_i, codon_i = index_map[flat_index]
                codon_str = self.int_to_codon[codon_int]
                randomized_entries[entry_i]['codons'][codon_i] = codon_str
            return randomized_entries

    def get_dicodon_frequency(self, randomized_entries):
        """
        Compares the frequency of dicodons in the original and randomized CDS entries.
        Returns a DataFrame with the frequency of each dicodon in both sets.
        """
        dicodon_freq = defaultdict(int)
        # Count dicodon frequency in the original entries
        for entry in self.cds_entries:
            codons = entry['codons']
            for i in range(len(codons) - 1):
                dicodon = codons[i] + codons[i + 1]
                dicodon_freq[dicodon] += 1
        
        # Create a DataFrame to store the frequencies
        dicodon_df = pd.DataFrame.from_dict(dicodon_freq, orient='index', columns=['original_freq'])
        # Count dicodon frequency in the randomized entries
        randomized_dicodon_freq = defaultdict(int)
        for entry in randomized_entries:
            codons = entry['codons']
            for i in range(len(codons) - 1):
                dicodon = codons[i] + codons[i + 1]
                randomized_dicodon_freq[dicodon] += 1
        
        rand_df = pd.DataFrame.from_dict(randomized_dicodon_freq, orient='index', columns=['randomized_freq'])
        # Join the original and randomized frequencies
        dicodon_df = dicodon_df.join(rand_df, how='outer')
        dicodon_df.fillna(0, inplace=True)  # Fill NaN values with
        dicodon_df = dicodon_df.astype(int)
        dicodon_df.index.name = 'dicodon'  # Set the index name
        dicodon_df.reset_index(inplace=True)

        # print whether the dicodon frequencies match
        if dicodon_df['original_freq'].equals(dicodon_df['randomized_freq']):
            print("Dicodon frequencies match between original and randomized entries.")
        else:
            print("Dicodon frequencies do not match between original and randomized entries.")
            # print out rows where frequencies differ
            diff_df = dicodon_df[dicodon_df['original_freq'] != dicodon_df['randomized_freq']]
            print("Differences:")
            print(diff_df)
            # save differences to a file
            diff_df.to_csv('dicodon_differences.csv', index=False)
        return dicodon_df        

    def count_kmer_frequency(self, k, entries=None):
        """
        Counts the frequency of k-mers in the CDS entries.
        
        Parameters:
            k (int): The length of the k-mers to count.
        
        Returns:
            dict: A dictionary mapping each k-mer to its frequency.
        """
        kmer_freq = defaultdict(int)
        if entries is None:
            entries = self.cds_entries
        for entry in entries:
            seq = ''.join(entry['codons'])
            for i in range(len(seq) - k + 1):
                if k % 3 == 0 and i % 3 == 0: # only count out of frame kmers if k = 6
                    continue
                kmer = seq[i:i + k]# Only keep kmers with valid bases (A, T, C, G)
                if re.fullmatch(r"[ATCG]+", kmer):
                    kmer_freq[kmer] += 1
        return dict(kmer_freq)
    
    def get_kmer_freq_randomized(self, k, num_randomizations=20, max_swaps=None):
        """Generate a table of k-mer frequences where rows are kmers
        and columns are the k-mer frequencies in the original entries and
        in the randomized entries."""

        # count kmer frequency in the original entries
        original_kmer_freq = self.count_kmer_frequency(k)
        kmer_df = pd.DataFrame.from_dict(original_kmer_freq, orient='index', columns=['original'])

        # randomize the entries and count kmer frequency
        for i in range(num_randomizations):
            randomized_entries = self.run_mcmc(max_swaps=max_swaps)
            randomized_kmer_freq = self.count_kmer_frequency(k, randomized_entries)
            rand_df = pd.DataFrame.from_dict(randomized_kmer_freq, orient='index', columns=[f'random_{i+1}'])
            kmer_df = kmer_df.join(rand_df, how='outer')  # outer join to keep all k-mers

        kmer_df.fillna(0, inplace=True)
        kmer_df = kmer_df.astype(int)
        kmer_df.index.name = 'kmer'
        kmer_df.reset_index(inplace=True)

        return kmer_df
    
    def plot_kmer_freq_vs_nswaps(self, k, output_file, convergence_threshold=0.1):
        n_codons = sum(len(entry['codons']) for entry in self.cds_entries)
        n = 3 * n_codons
        max_swaps = int(3 * n * math.log(n))

        # Your defined checkpoints
        checkpoints = sorted(set([
            1, int(n), int(n * math.log(math.log(n))), int(n * math.log(n)/5),
            int(n * math.log(n)/4), int(n * math.log(n)/3), int(n * math.log(n)/2),
            int(n * math.log(n) * 0.75), 
            int(n * math.log(n)), int(1.5 * n * math.log(n)), int(n * math.log(n)) * 2, max_swaps
        ]))
        checkpoints = [c for c in checkpoints if c > 0 and c <= max_swaps]

        # Run MCMC
        randomized_entries_by_checkpoint = self.run_mcmc(seeded=True, checkpoints=checkpoints)

        # Compute k-mer frequencies
        freq_snapshots = []
        for swap_count in checkpoints:
            #print(f"Counting k-mers after {swap_count} swaps")
            entries = randomized_entries_by_checkpoint[swap_count]
            freq = self.count_kmer_frequency(k, entries)
            freq_snapshots.append(freq)

        kmers = sorted(set().union(*[f.keys() for f in freq_snapshots]))
        
        # Calculate absolute delta in k-mer frequencies
        y_vals = []

        for i in range(1, len(freq_snapshots)):
            prev = freq_snapshots[i - 1]
            curr = freq_snapshots[i]
            prev_arr = np.array([prev.get(kmer, 0) for kmer in kmers])
            curr_arr = np.array([curr.get(kmer, 0) for kmer in kmers])
            delta = np.sum(np.abs(curr_arr - prev_arr))
            y_vals.append(delta)

        # Compute relative change in delta values
        rel_delta_changes = []
        for i in range(1, len(y_vals)):
            delta_prev = y_vals[i - 1]
            delta_curr = y_vals[i]
            rel_change = abs(delta_curr - delta_prev) / max(1, delta_prev)
            rel_delta_changes.append(rel_change)

        # Determine convergence as point after which rel_change stays below threshold
        converged_at = None
        for i in range(len(rel_delta_changes)):
            if all(rc < convergence_threshold for rc in rel_delta_changes[i:]):
                converged_at = checkpoints[i + 2]  # +2 as before
                break
        
        # Return symbolic label for convergence
        def label_checkpoint(val):
            candidates = {
                '1': 1,
                'n': int(n),
                'nlog(log(n))': int(n * math.log(math.log(n))),
                'nlog(n)/5': int(n * math.log(n) / 5),
                'nlog(n)/4': int(n * math.log(n) / 4),
                'nlog(n)/3': int(n * math.log(n) / 3),
                'nlog(n)/2': int(n * math.log(n) / 2),
                '3nlog(n)/4': int(0.75 * n * math.log(n)),
                'nlog(n)': int(n * math.log(n)),
                '3nlog(n)/2': int(1.5 * n * math.log(n)),
                '2nlog(n)': int(2 * n * math.log(n)),
                'max_swaps': max_swaps
            }
            for label, ref in candidates.items():
                if abs(val - ref) / max(1, ref) < 0.05:  # within 5%
                    return label
            return str(val)

        if converged_at is not None:
            #print(f"Delta stabilized at {label_checkpoint(converged_at)} swaps")
            pass
        else:
            print("Delta never fully stabilized")

        # Plot
        plt.figure(figsize=(8, 5))
        plt.plot(checkpoints[1:], y_vals, marker='o')
        plt.xlabel("Number of swaps")
        plt.ylabel("Δ total k-mer counts vs previous")
        plt.title(f"Convergence of k-mer frequencies (k={k})")
        plt.grid(True)

        if converged_at is not None:
            plt.axvline(converged_at, color='red', linestyle='--', label=f'Converged at {converged_at}')
            plt.legend()

        plt.savefig(output_file, dpi=1200)
        plt.close()

        return label_checkpoint(converged_at) if converged_at is not None else None

@njit(inline="always")
def encode_context(prev, aa, next_, n_codons, n_aas):
    """
    Hashes codon context (previous codon, amino acid translation, next codon)
    """
    return prev * n_aas * n_codons + aa * n_codons + next_

@njit
def _build_context_map(codon_seq, aa_seq, n_codons, n_aas, max_pos_per_context, valid_swap_pos):
    """
    Constructs a codon context map (max possible contexts, max positions per context)
    * max_pos_per_context is a heuristic estimate which is several fold greater than the initial max positions per context
    * context idx (row index) is generated through the hashed encode_context function
    * loops through all codons across all segments and updates context map at the context row with the position index of the matching codon

    Maintains a context_counts (max possible contexts, 1) array which tracks how many positions are recorded per context in the map  
    """
    max_contexts = n_codons * n_aas * n_codons
    context_map_positions = -1 * np.ones((max_contexts, max_pos_per_context), dtype=np.int32)
    context_counts = np.zeros(max_contexts, dtype=np.int32)

    # go through all positions in the 1d codon array, covers all coding segments
    for i in range(1, len(codon_seq) - 1):
        if not valid_swap_pos[i]:
            continue
        prev = codon_seq[i - 1]
        aa = aa_seq[i]
        next_ = codon_seq[i + 1]
        ctx_id = encode_context(prev, aa, next_, n_codons, n_aas)

        cnt = context_counts[ctx_id]
        if cnt < max_pos_per_context:
            context_map_positions[ctx_id, cnt] = i
            context_counts[ctx_id] = cnt + 1
        
    return context_map_positions, context_counts

@njit(inline="always")
def _remove_position_from_context_by_ctx(pos, ctx_id, context_map_positions, context_counts, valid_swap_pos=None):
    """
    Given a context and a position, searches for the entry in context map mapping the context to position and removes it
    Updates context counts accordingly
    """
    cnt = context_counts[ctx_id]
    for idx in range(cnt):
        if context_map_positions[ctx_id, idx] == pos:
            last_idx = cnt - 1
            # swap the position with the last one
            context_map_positions[ctx_id, idx] = context_map_positions[ctx_id, last_idx]
            # remove the last position after swapping
            context_map_positions[ctx_id, last_idx] = -1
            context_counts[ctx_id] = cnt - 1
            return
        
    print("Warning: position not found in context during removal.")

    # search for the position in the entire context map
    for ctx in range(len(context_map_positions)):
        for idx in range(max(context_counts[ctx], 0)):
            if context_map_positions[ctx, idx] == pos:
                print("Found position in context.")
                return
              
@njit(inline="always")
def _add_position_to_context_by_ctx(pos, ctx_id, context_map_positions, context_counts, max_pos_per_context):
    """
    Adds a new context-position pair to the map and updates context_counts
    """
    cnt = context_counts[ctx_id]
    if cnt < max_pos_per_context:
        context_map_positions[ctx_id, cnt] = pos
        context_counts[ctx_id] = cnt + 1
    
@njit
def _mcmc_codons(
    codon_seq, aa_seq, n_codons, n_aas, n_swaps, max_pos_per_context, valid_swap_pos, seeded=False
):
    """
    Runs the MCMC simulation for `n_swaps` iterations, does not swap adjacent codons or start/end codons
    `codon_seq` is a (1, total # of codons) array which concatenates the codons across all the segments continuously
    `aa_seq` is constructed similarly

    valid_swap_pos is a Boolean array of shape (1, total # of codons) which indicates which codons are eligible for swapping. 
    * codons that are at the start or end of the coding segment are marked as invalid for swapping
    """
    if seeded:
        np.random.seed(42)
    codon_seq = codon_seq.copy()
    aa_seq = aa_seq.copy()

    # Only use valid positions to build the context map
    context_map_positions, context_counts = _build_context_map(
        codon_seq, aa_seq, n_codons, n_aas, max_pos_per_context, valid_swap_pos
    )
    n_contexts = n_codons * n_aas * n_codons

    for swap_i in range(n_swaps):
        ctx = np.random.randint(n_contexts)
        if context_counts[ctx] < 2:
            continue

        cnt = context_counts[ctx]
        pos_list = context_map_positions[ctx, :cnt]

        # Choose two positions from valid list in this context
        i, j = np.random.choice(pos_list, 2, replace=False)

        # Ensure positions are valid for swapping, i.e. not directly adjacent
        if abs(i - j) <= 1 or not valid_swap_pos[i]:
            continue

        positions = [i - 1, i, i + 1, j - 1, j, j + 1]
        visited_positions = []
        # Remove old context for surrounding positions, if within bounds
        for pos in positions:
            if pos in visited_positions:
                continue
            if pos <= 0 or pos >= len(codon_seq) - 1:
                continue
            if not valid_swap_pos[pos]:
                continue
            prev = codon_seq[pos - 1]
            aa = aa_seq[pos]
            next_ = codon_seq[pos + 1]
            old_ctx = encode_context(prev, aa, next_, n_codons, n_aas)
            _remove_position_from_context_by_ctx(pos, old_ctx, context_map_positions, context_counts, valid_swap_pos=valid_swap_pos)

            visited_positions.append(pos)

        # Perform the codon swap
        codon_seq[i], codon_seq[j] = codon_seq[j], codon_seq[i]

        visited_positions = []
        # Add new context for surrounding positions
        for pos in positions:
            if pos in visited_positions:
                continue
            if pos <= 0 or pos >= len(codon_seq) - 1:
                continue
            if not valid_swap_pos[pos]:
                continue
            prev = codon_seq[pos - 1]
            aa = aa_seq[pos]
            next_ = codon_seq[pos + 1]
            new_ctx = encode_context(prev, aa, next_, n_codons, n_aas)
            _add_position_to_context_by_ctx(
                pos,
                new_ctx,
                context_map_positions,
                context_counts,
                max_pos_per_context,
            )
            visited_positions.append(pos)

    return codon_seq

@njit
def _mcmc_codons_with_checkpoints(
    codon_seq, aa_seq, n_codons, n_aas, checkpoints, max_pos_per_context, valid_swap_pos, seeded=False
):
    checkpoints = np.sort(checkpoints)
    
    codon_seq = codon_seq.copy()
    aa_seq = aa_seq.copy()

    context_map_positions, context_counts = _build_context_map(
        codon_seq, aa_seq, n_codons, n_aas, max_pos_per_context, valid_swap_pos
    )
    n_contexts = n_codons * n_aas * n_codons

    # Declare results as a typed dict: keys are int, values are int arrays
    results = Dict.empty(key_type=types.int64, value_type=types.int32[:])
    if seeded:
        np.random.seed(42)

    checkpoint_index = 0

    for swap_i in range(checkpoints[-1] + 1):

        # check if we reached a checkpoint before attempting a swap - this avoids missing checkpoints if there is no viable swap
        if checkpoint_index < len(checkpoints) and swap_i == checkpoints[checkpoint_index]:
            # Save a copy of the current codon sequence
            #print(f"Checkpoint {checkpoint_index}: {swap_i} swaps")
            results[swap_i] = codon_seq.copy()
            checkpoint_index += 1

        ctx = np.random.randint(n_contexts)
        if context_counts[ctx] < 2:
            continue

        cnt = context_counts[ctx]
        pos_list = context_map_positions[ctx, :cnt]
        i, j = np.random.choice(pos_list, 2, replace=False)
        if abs(i - j) <= 1 or not valid_swap_pos[i]:
            continue

        positions = [i - 1, i, i + 1, j - 1, j, j + 1]
        visited_positions = set()
        for pos in positions:
            if pos in visited_positions or pos <= 0 or pos >= len(codon_seq) - 1 or not valid_swap_pos[pos]:
                continue
            prev = codon_seq[pos - 1]
            aa = aa_seq[pos]
            next_ = codon_seq[pos + 1]
            old_ctx = encode_context(prev, aa, next_, n_codons, n_aas)
            _remove_position_from_context_by_ctx(pos, old_ctx, context_map_positions, context_counts)
            visited_positions.add(pos)

        codon_seq[i], codon_seq[j] = codon_seq[j], codon_seq[i]

        visited_positions = set()
        for pos in positions:
            if pos in visited_positions or pos <= 0 or pos >= len(codon_seq) - 1 or not valid_swap_pos[pos]:
                continue
            prev = codon_seq[pos - 1]
            aa = aa_seq[pos]
            next_ = codon_seq[pos + 1]
            new_ctx = encode_context(prev, aa, next_, n_codons, n_aas)
            _add_position_to_context_by_ctx(pos, new_ctx, context_map_positions, context_counts, max_pos_per_context)
            visited_positions.add(pos)

    return results

def main():
    parser = argparse.ArgumentParser(description="Calculate kmer frequencies in ground-truth and randomized genomes.")

    parser.add_argument('cds_file', type=str, help='Path to FASTA file containing coding segment sequences.')
    parser.add_argument('translated_file', type=str, help='Path to FASTA file containing AA translations of coding segment sequences.')
    parser.add_argument('output_dir', type=str, help='Path to output directory.')

    # add mode indicating how to run the randomizer, either to generate plots of freq vs nswaps or to generate kmer frequencies
    parser.add_argument('--mode', type=str, choices=['plot', 'kmer', 'kmer_limited'], default='kmer',
                        help='Mode of operation: "plot" to generate plots of kmer frequency vs number of swaps, "kmer" to generate kmer frequencies.')
    args = parser.parse_args()

    cds_file = args.cds_file
    translated_file = args.translated_file

    cds_randomizer = CdsRandomizer()
    cds_randomizer.read_cds_entries(cds_file, translated_file)
    #print(f"Read {len(cds_randomizer.cds_entries)} CDS entries from {cds_file} and {translated_file}")

    # print out some example entries
    #for entry in cds_randomizer.cds_entries[:3]:
    #    print(f"id: {entry['id']} | seq: {entry['codons']} | AA: {entry['aa']}")

    species_name = cds_file.split('/')[-2]  # assuming the species name is the parent directory of the cds_file
    if args.mode == 'plot':
        converged_at = cds_randomizer.plot_kmer_freq_vs_nswaps(k=6, output_file=args.output_dir + f'{species_name}_kmer_freq_vs_nswaps_k6.png')
        print(f"{species_name} k=6 converged at {converged_at} swaps")
    else:
        randomized_entries = cds_randomizer.run_mcmc()

        # qc randomized genome to check that dicodon frequencies are not changed
        dicodon_df = cds_randomizer.get_dicodon_frequency(randomized_entries)
        print(dicodon_df.head())
        
        #kmer_freq = cds_randomizer.get_kmer_freq_randomized(k=5, num_randomizations=40)
        #kmer_freq.to_csv(args.output_dir + 'kmer_frequencies_k5.csv', index=False)

        if args.mode == 'kmer_limited':
            n_codons = sum(len(entry['codons']) for entry in cds_randomizer.cds_entries)
            n = 3 * n_codons
            max_swaps = int(0.5 * n * math.log(n))
        else: 
            max_swaps = None

        kmer_freq = cds_randomizer.get_kmer_freq_randomized(k=6, num_randomizations=20, max_swaps=max_swaps)
        kmer_freq.to_csv(args.output_dir + f'{species_name}_kmer_frequencies_k6.csv', index=False)

        #kmer_freq = cds_randomizer.get_kmer_freq_randomized(k=7, num_randomizations=5)
        #kmer_freq.to_csv(args.output_dir + f'{species_name}_kmer_frequencies_k7.csv', index=False)

if __name__ == "__main__":
    main()
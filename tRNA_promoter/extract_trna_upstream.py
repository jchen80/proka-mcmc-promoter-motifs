"""
Author: Julie Chen
extract_trna_upstream.py

this script takes in the output of the tRNAscan-SE results and extracts the putative promoter regions
for the detected tRNAs, which is defined as the 100bp upstream of the gene
"""

from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

UPSTREAM_LENGTH = 100  # bp upstream to extract

def overlaps(interval1, interval2):
    """Check if two [start, end] intervals overlap"""
    return max(interval1[0], interval2[0]) <= min(interval1[1], interval2[1])


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 extract_trna_upstream.py <base_dir> <output_dir>")
        sys.exit(1)

    base_dir = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])
    output_dir.mkdir(exist_ok=True)

    for species_dir in base_dir.iterdir():
        if not species_dir.is_dir():
            continue

        trnascan_files = list(species_dir.glob("*.trnascan.tsv"))
        fna_files = [
            f for f in species_dir.glob("*.fna")
            # filter out files that contain "cds_from" in their name
            if "cds_from" not in f.name
        ]
        if not trnascan_files or not fna_files:
            continue

        trna_file = trnascan_files[0]
        genome_file = fna_files[0]
        species_name = species_dir.name

        genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
        upstream_records = []

        # collect all tRNA intervals
        trna_intervals = {}  # {contig_id: list of (start, end)}
        trna_entries = []
        with open(trna_file) as f:
            for line in f:
                if line.startswith("Name") or line.startswith("Sequence") or line.startswith("--------") or line.strip() == "":
                    continue
                if "pseudo" in line: # skip pseudo tRNAs
                    continue

                fields = line.strip().split()
                contig_id = fields[0]
                tRNA_type = fields[1]
                start = int(fields[2])
                end = int(fields[3])

                # keeps the start and end in same order in order to differentiate strand later
                trna_entries.append((contig_id, tRNA_type, start, end))

                if contig_id not in trna_intervals:
                    trna_intervals[contig_id] = []
                trna_intervals[contig_id].append((min(start, end), max(start, end)))  # normalized interval for finding overlaps

        num_minus = 0
        num_skipped_overlap = 0
        
        # loop through tRNA gene intervals, check if there are overlaps with other tRNAs, if not, extract upstream sequence
        for contig_id, tRNA_type, start, end in trna_entries: 
            strand = "+" if start < end else "-"
            trna_start, trna_end = min(start, end), max(start, end)

            if strand == "+":
                upstream_start = max(0, trna_start - UPSTREAM_LENGTH - 1)
                upstream_end = trna_start - 1
            else:
                upstream_start = trna_end
                upstream_end = min(len(genome[contig_id]) - 1, trna_end + UPSTREAM_LENGTH)
                num_minus += 1

            upstream_interval = (upstream_start + 1, upstream_end + 1)  # convert to 1-based for comparison

            # Check for overlap with other tRNAs
            overlaps_any = any(
                overlaps(upstream_interval, other)
                for other in trna_intervals[contig_id]
                if other != (trna_start, trna_end)
            )
            if overlaps_any:
                num_skipped_overlap += 1
                continue

            if strand == "+":
                upstream_seq = genome[contig_id].seq[upstream_start:upstream_end]
            else:
                upstream_seq = genome[contig_id].seq[upstream_start:upstream_end].reverse_complement()

            if len(upstream_seq) < UPSTREAM_LENGTH:
                continue  # Skip incomplete upstream sequences

            record_id = f"{species_name}|{contig_id}|{tRNA_type}|{strand}|{start}-{end}"
            upstream_records.append(SeqRecord(upstream_seq, id=record_id, description=""))

        if upstream_records:
            out_path = output_dir / f"{species_name}_tRNA_promoter.fasta"
            SeqIO.write(upstream_records, out_path, "fasta")
            print(f"Wrote {len(upstream_records)} upstreams → {out_path}")
            print(f"Number of minus strand tRNAs: {num_minus}")
            print(f"Number skipped due to overlap: {num_skipped_overlap}")
        else:
            print(f"No upstreams found for {species_name}")


if __name__ == "__main__":
    main()


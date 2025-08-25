"""
Author: Julie Chen
run_trnascan_all_species.py

this script runs tRNAscan-SE on all archaeal species in the 'archaea' directory
uses the *_genomic.fna file as the genomic sequence input
"""

import subprocess
from pathlib import Path
import sys

def main():

    if len(sys.argv) != 2:
        print("Usage: python3 run_trnascan_all_species.py <parent_dir>")
        sys.exit(1)

    base_dir = Path(sys.argv[1])
    if not base_dir.is_dir():
        print(f"Error: {base_dir} is not a valid directory.")
        sys.exit(1)

    for species_dir in base_dir.iterdir():
        if not species_dir.is_dir():
            continue

        # Look for the correct *_genomic.fna.gz file only
        fna_files = [
            f for f in species_dir.glob("*_genomic.fna")
            if "cds_from" not in f.name
        ]
        if not fna_files:
            print(f"No .fna file found in {species_dir}")
            continue
        fna_file = fna_files[0]  # Assuming one per species

        # Run tRNAscan
        trna_output = species_dir / f"{fna_file.stem}.trnascan.tsv"
        # Check if output already exists
        if trna_output.exists():
            print(f"tRNA scan output already exists for {species_dir}, skipping.")
            continue

        print(f"Running tRNAscan-SE on {fna_file.name}")
        
        try:
            if base_dir.name == "archaea":
                subprocess.run(
                    ["tRNAscan-SE", "-A", "-o", str(trna_output.with_suffix(".tsv")), str(fna_file)],
                check=True
                )
            elif base_dir.name == "bacteria":
                subprocess.run(
                    ["tRNAscan-SE", "-B", "-o", str(trna_output.with_suffix(".tsv")), str(fna_file)],
                check=True
                )
        except subprocess.CalledProcessError:
            print(f"tRNAscan failed on {fna_file}")


if __name__ == "__main__":
    main()
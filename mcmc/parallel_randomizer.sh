#!/bin/bash

BASE_DIR="paxdb_bacteria"

run_for_species() {
  SPECIES_DIR="$1"

  # Skip if kmer_frequencies_k6.csv already exists
  if compgen -G "$SPECIES_DIR/*kmer_frequencies_k6.csv" > /dev/null; then
    echo "Skipping $SPECIES_DIR — kmer_frequencies_k6.csv already exists"
    return
  fi

  CDS_FILE=$(find "$SPECIES_DIR" -name '*cds_from_genomic.fna' | head -n 1)
  TRANS_FILE=$(find "$SPECIES_DIR" -name '*translated_cds.faa' | head -n 1)

  if [[ -f "$CDS_FILE" && -f "$TRANS_FILE" ]]; then
    echo "Running: $SPECIES_DIR"
    python3 mcmc/cds_randomizer.py "$CDS_FILE" "$TRANS_FILE" "$SPECIES_DIR/" --mode kmer
  else
    echo "Skipping $SPECIES_DIR due to missing files"
  fi
}
export -f run_for_species

export BASE_DIR
find "$BASE_DIR" -mindepth 1 -maxdepth 1 -type d | parallel --ungroup -j 16 run_for_species

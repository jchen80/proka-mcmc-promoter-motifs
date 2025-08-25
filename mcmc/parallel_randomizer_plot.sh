#!/bin/bash

BASE_DIR="archaea"

run_for_species() {
  SPECIES_DIR="$1"
  CDS_FILE=$(find "$SPECIES_DIR" -name '*_cds_from_genomic.fna' | head -n 1)
  TRANS_FILE=$(find "$SPECIES_DIR" -name '*_translated_cds.faa' | head -n 1)

  if [[ -f "$CDS_FILE" && -f "$TRANS_FILE" ]]; then
    #echo "Running: $SPECIES_DIR"
    python3 mcmc/cds_randomizer.py "$CDS_FILE" "$TRANS_FILE" "$SPECIES_DIR/" --mode plot
  else
    echo "Skipping $SPECIES_DIR due to missing files"
  fi
}
export -f run_for_species

export BASE_DIR
find "$BASE_DIR" -mindepth 1 -maxdepth 1 -type d | parallel --ungroup -j 16 run_for_species

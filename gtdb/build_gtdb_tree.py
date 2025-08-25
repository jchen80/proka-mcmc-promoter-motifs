"""
Author: Julie Chen
build_gtdb_tree.py

This file constructs a phylogenetic tree of the given species by fetching the relevant entries 
from either the ar53 or bac120 GTDB reference trees and trimming accordingly
"""

import argparse
import pandas as pd
import sys
from ete3 import Tree
from tqdm import tqdm

def load_list(species_csv, assembly_col, species_col):
    """
    loads in df which contains Refseq assembly accessions and species' names
    returns a dictionary mapping assembly accession ID to species name, and aforementioned dataframe
    """
    df = pd.read_csv(species_csv)
    df["assembly_id"] = "RS_" + df[assembly_col]
    assembly_to_species = dict(zip(df["assembly_id"], df[species_col]))
    return assembly_to_species, df

def load_gtdb_metadata(gtdb_metadata_path):
    meta = pd.read_csv(gtdb_metadata_path, sep='\t', low_memory=False)
    return meta

def match_metadata(assembly_to_species, gtdb_meta):
    assembly_list = list(assembly_to_species.keys())
    matched_gtdb = gtdb_meta[gtdb_meta["accession"].isin(assembly_list)].copy()
    matched_gtdb["input_species"] = matched_gtdb["accession"].map(assembly_to_species)

    unmatched = set(assembly_list) - set(matched_gtdb["accession"])
    if unmatched:
        print("Warning: Could not match the following assemblies to GTDB metadata:")
        for s in unmatched:
            print(f"  - {s}")
    
    return matched_gtdb

def build_tree(matched_metadata):
    """
    build phylogenetic tree given the relevant trimmed metadata information about taxonomy
    """
    tree = Tree()
    nodes = {"root": tree}

    for _, row in tqdm(matched_metadata.iterrows()):
        lineage = row["gtdb_taxonomy"].split(";")
        lineage = [x.split("__")[1] for x in lineage if x]
        
        parent = tree
        path = []
        for taxon in lineage:
            path.append(taxon)
            key = "/".join(path)
            if key not in nodes:
                nodes[key] = parent.add_child(name=taxon)
            parent = nodes[key]

        # Add species tip with user-provided species name
        parent.add_child(name=row["input_species"])

    return tree

def main():
    parser = argparse.ArgumentParser(description="Build phylogenetic tree from GTDB metadata and species list, matching GTDB and NCBI taxonomy species.")
    parser.add_argument("--csv_file", required=True, help="CSV file with species data.")
    parser.add_argument("--assembly_col", default="assembly", help="Column name with assembly accession")
    parser.add_argument("--species_col", default="species_reformatted", help="Column name with species (e.g. F. placidus)")
    parser.add_argument("--gtdb_metadata", required=True, help="GTDB metadata TSV (ar122 or bac120).")
    
    args = parser.parse_args()

    assembly_list, df = load_list(args.csv_file, args.assembly_col, args.species_col)
    gtdb_meta = load_gtdb_metadata(args.gtdb_metadata)
    matched = match_metadata(assembly_list, gtdb_meta)

    if matched.empty:
        print("No species matched between your CSV and GTDB metadata.")
        sys.exit(1)

    output_metadata = args.csv_file.replace('.csv', '_matched_metadata.tsv')
    matched.to_csv(output_metadata, sep='\t', index=False)

    tree = build_tree(matched)
    output_tree = args.csv_file.replace('.csv', '_tree.nwk')
    tree.write(outfile=output_tree)

    print(f"Tree saved to {output_tree}")
    print(f"Matched metadata saved to {output_metadata}")

if __name__ == "__main__":
    main()


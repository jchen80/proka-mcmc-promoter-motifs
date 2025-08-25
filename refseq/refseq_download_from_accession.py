"""
Author: Julie Chen
refseq_download_from_url.py

Given an input csv file that contains refseq accession codes, 
download the cds_from_genomic, genomic, and translated_cds files
"""

import os
import sys
import pandas as pd
import requests
import gzip
import shutil

def download_and_extract(url, out_path):
    """Download and unzip a .gz file from url into out_path."""
    filename = os.path.basename(url)
    gz_path = os.path.join(out_path, filename)
    out_file = gz_path[:-3]  # remove .gz extension

    # Skip if already downloaded and extracted
    if os.path.exists(out_file):
        print(f"Already exists: {out_file}")
        return

    # Download
    print(f"Downloading: {url}")
    r = requests.get(url, stream=True)
    if r.status_code != 200:
        print(f"Failed to download: {url}")
        return
    with open(gz_path, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)

    # Unzip
    with gzip.open(gz_path, 'rb') as f_in, open(out_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(gz_path)  # delete the .gz file

    print(f"Extracted to: {out_file}")

def main(csv_file, output_root="output"):
    """Main function to process the CSV and download/extract files into output_root/species_name/"""
    # Make the output root directory
    os.makedirs(output_root, exist_ok=True)

    df = pd.read_csv(csv_file)

    for _, row in df.iterrows():
        acc_code = row['Assembly']
        species = row['Specie'].replace(" ", "_")  # Clean species name

        # Path to species-specific directory inside the output root
        species_dir = os.path.join(output_root, species)
        os.makedirs(species_dir, exist_ok=True)

        base_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/{acc_code[4:7]}/{acc_code[7:10]}/{acc_code[10:13]}"
        print(f"Base URL: {base_url}")
        
        code_no_version = acc_code[4:].split(".")[0].lstrip("0")[:-1]
        version = acc_code.split(".")[1]
        base_name = f"GCF_{acc_code[4:]}_ASM{code_no_version}v{version}"
        cds_url = f"{base_url}/{base_name}/{base_name}_cds_from_genomic.fna.gz"
        genomic_url = f"{base_url}/{base_name}/{base_name}_genomic.fna.gz"
        translated_url = f"{base_url}/{base_name}/{base_name}_translated_cds.faa.gz"

        download_and_extract(cds_url, species_dir)
        download_and_extract(genomic_url, species_dir)
        download_and_extract(translated_url, species_dir)

if __name__ == "__main__":
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python3 refseq_download_from_accession.py <input_csv_file> [output_directory]")
        sys.exit(1)

    csv_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) == 3 else "output"
    main(csv_file, output_dir)

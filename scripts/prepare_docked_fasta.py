#!/usr/bin/env python3
"""
Extract FASTA sequences for docked candidates from filtered_hmmsearch.faa
"""

import pandas as pd
from Bio import SeqIO
import os

BASE_DIR = os.getcwd()

# Read VINA summary (docked candidates)
vina_summary = os.path.join(BASE_DIR, "output_files/vina_summary.csv")
df_vina = pd.read_csv(vina_summary)

if df_vina.empty:
    print("No docked candidates found. Creating empty FASTA file.")
    output_fasta = os.path.join(BASE_DIR, "data/fasta_files/docked_hits.faa")
    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
    open(output_fasta, 'w').close()
    exit(0)

# Get list of docked proteins
docked_proteins = df_vina['protein'].tolist()
print(f"Found {len(docked_proteins)} docked proteins")

# Extract sequences from filtered_hmmsearch.faa
input_fasta = os.path.join(BASE_DIR, "data/fasta_files/filtered_hmmsearch.faa")
output_fasta = os.path.join(BASE_DIR, "data/fasta_files/docked_hits.faa")

os.makedirs(os.path.dirname(output_fasta), exist_ok=True)

if not os.path.exists(input_fasta):
    print(f"WARNING: Input FASTA not found: {input_fasta}")
    print("Creating empty output FASTA")
    open(output_fasta, 'w').close()
    exit(0)

# Extract sequences
extracted = 0
with open(output_fasta, 'w') as out_handle:
    for record in SeqIO.parse(input_fasta, "fasta"):
        if record.id in docked_proteins:
            SeqIO.write(record, out_handle, "fasta")
            extracted += 1

print(f"Extracted {extracted} sequences to {output_fasta}")

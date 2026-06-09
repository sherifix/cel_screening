#!/usr/bin/env python3
"""
Calculate isoelectric point (pI) for protein sequences and merge with EpHod predictions
"""

import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
import os

BASE_DIR = os.getcwd()

# Input files
input_fasta = os.path.join(BASE_DIR, "data/fasta_files/docked_hits.faa")
input_top_hits = os.path.join(BASE_DIR, "output_files/top_hits.csv")
input_ephod = os.path.join(BASE_DIR, "output_files/opt_ph_prediction.csv")
output_csv = os.path.join(BASE_DIR, "output_files/ph_pi_summary.csv")

# Check if input FASTA exists
if not os.path.exists(input_fasta):
    print(f"ERROR: FASTA file not found: {input_fasta}")
    print("Creating empty output file")
    pd.DataFrame(columns=["accession", "pI", "length", "molecular_weight", "optimal_ph"]).to_csv(output_csv, index=False)
    exit(0)

# Calculate pI for each sequence
print("Calculating pI for each candidate...")
pi_results = []
for record in SeqIO.parse(input_fasta, "fasta"):
    accession = record.id.split()[0]
    sequence = str(record.seq)
    
    try:
        analysis = ProteinAnalysis(sequence)
        pi = analysis.isoelectric_point()
        length = len(sequence)
        mw = analysis.molecular_weight()
        pi_results.append({
            "accession": accession,
            "pI": round(pi, 2),
            "length": length,
            "molecular_weight": round(mw, 2)
        })
    except Exception as e:
        print(f"Error calculating pI for {accession}: {e}")
        pi_results.append({
            "accession": accession,
            "pI": None,
            "length": len(sequence),
            "molecular_weight": None
        })

df_pi = pd.DataFrame(pi_results)
print(f"Calculated pI for {len(df_pi)} sequences")

# Load EpHod predictions if available
df_ephod = None
if os.path.exists(input_ephod):
    df_ephod = pd.read_csv(input_ephod)
    print(f"Loaded EpHod predictions for {len(df_ephod)} sequences")
else:
    print("Warning: EpHod output not found. Will create without optimal pH.")

# Load top hits if available
df_top_hits = None
if os.path.exists(input_top_hits):
    df_top_hits = pd.read_csv(input_top_hits)
    print(f"Loaded top hits for {len(df_top_hits)} sequences")

# Start with pI results
df_merged = df_pi

# Merge with top hits (add family, TM-score, etc.)
if df_top_hits is not None:
    # Select relevant columns from top_hits
    top_hits_cols = ['accession', 'family']
    available_cols = [col for col in top_hits_cols if col in df_top_hits.columns]
    df_merged = df_merged.merge(df_top_hits[available_cols], on="accession", how="left")
    print(f"Merged with top_hits.csv")

# Merge with EpHod predictions (optimal pH)
if df_ephod is not None:

    if 'Header' in df_ephod.columns:
        df_ephod['accession'] = df_ephod['Header'].str.split().str[0]
    elif 'accession' not in df_ephod.columns:
        # Try to extract from first column
        first_col = df_ephod.columns[0]
        df_ephod['accession'] = df_ephod[first_col].str.split().str[0]
    
    df_ephod = df_ephod.round(2)

    # Select pH-related columns (adjust column names as needed)
    ph_cols = ['accession']
    for col in df_ephod.columns:
        if 'ensemble' in col.lower() or 'svr' in col.lower() or 'rlattr' in col.lower():
            ph_cols.append(col)
    
    df_merged = df_merged.merge(df_ephod[ph_cols], on="accession", how="left")
    print(f"Merged with EpHod predictions")

# Reorder columns for better readability
desired_order = ['accession', 'family', 'pI', 'Ensemble', 'SVR', 'RLATtr','molecular_weight', 'length']
existing_cols = [col for col in desired_order if col in df_merged.columns]
other_cols = [col for col in df_merged.columns if col not in desired_order]
df_merged = df_merged[existing_cols + other_cols]

# Save final summary
df_merged.to_csv(output_csv, index=False)
print(f"\n=== Summary ===")
print(f"Total candidates: {len(df_merged)}")
print(f"Saved to: {output_csv}")


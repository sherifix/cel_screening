#!/usr/bin/env python3
"""
Parse AutoDock VINA output logs and generate summary CSV
"""

import re
import pandas as pd
import os
from pathlib import Path

BASE_DIR = os.getcwd()

# Paths
log_dir = Path(BASE_DIR) / "results/autodock_vina/docking_results"
output_file = Path(BASE_DIR) / "output_files/vina_summary.csv"

# Check if log directory exists
if not log_dir.exists():
    print(f"WARNING: Log directory not found: {log_dir}")
    print("Creating empty output file")
    pd.DataFrame(columns=["protein", "vina_score", "log_file"]).to_csv(output_file, index=False)
    exit(0)

# Find all log files
log_files = list(log_dir.glob("*.log"))
print(f"Found {len(log_files)} log files")

if len(log_files) == 0:
    print("No log files found. Creating empty output file")
    pd.DataFrame(columns=["protein", "vina_score", "log_file"]).to_csv(output_file, index=False)
    exit(0)

rows = []

for log_file in log_files:
    # Extract protein name from filename (e.g., AAA34210.1.log -> AAA34210.1)
    protein = log_file.stem
    
    # Read log file
    try:
        with open(log_file, 'r') as f:
            text = f.read().splitlines()
    except Exception as e:
        print(f"Error reading {log_file}: {e}")
        continue
    
    best_score = None

    for line in text:
        # Match lines starting with mode number (1, 2, 3, etc.)
        match = re.match(r"^\s*1\s+", line)
        if match:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    best_score = float(parts[1])
                except ValueError:
                    continue
            break
    
    rows.append({
        "protein": protein,
        "vina_score": best_score,
        "log_file": str(log_file)
    })

# Create DataFrame
df = pd.DataFrame(rows)

# Filter out None scores
df = df[df['vina_score'].notna()]

# Sort by score (more negative = better binding)
if len(df) > 0:
    df = df.sort_values(by='vina_score', ascending=True)
    
    # Keep best score per protein (in case of multiple pockets)
    df = df.drop_duplicates(subset='protein', keep='first').reset_index(drop=True)

# Save output
df.to_csv(output_file, index=False)
print(f" Saved {len(df)} docking results to {output_file}")

# Print summary
print(f"\n======== VINA Docking Summary ========")
print(f"Total proteins docked: {len(df)}")
if len(df) > 0:
    print(f"\nTop 5 binding affinities:")
    print(df[['protein', 'vina_score']].head(10).to_string(index=False))
else:
    print("No docking results found")

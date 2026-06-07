import pandas as pd
from Bio import PDB
import numpy as np
import os
from pathlib import Path

BASE_DIR = os.getcwd()

# Paths
pdb_list_file = os.path.join(BASE_DIR, "data/trimmed_pdb_list.ds")
trimming_summary = os.path.join(BASE_DIR, "data/trimmed_structures/trimming_summary.csv")
predictions_dir = Path(BASE_DIR) / "results/p2rank_predictions"
output_file = os.path.join(BASE_DIR, "output_files/pockets_summary.tsv")

# Check if files exist
if not os.path.exists(trimming_summary):
    print(f"ERROR: Trimming summary not found: {trimming_summary}")
    exit(1)

if not predictions_dir.exists():
    print(f"ERROR: P2Rank predictions not found: {predictions_dir}")
    exit(1)

# Read trimming summary
df_trimming = pd.read_csv(trimming_summary)

# Parser for PDB files
parser = PDB.PDBParser(QUIET=True)

results = []

# Read PDB paths from list file
with open(pdb_list_file, 'r') as f:
    pdb_paths = [line.strip() for line in f if line.strip()]

print(f"Found {len(pdb_paths)} PDB files in list")

for i, pdb_path in enumerate(pdb_paths):
    print(f"\n[{i+1}/{len(pdb_paths)}] Processing: {os.path.basename(pdb_path)}")
    
    if not os.path.exists(pdb_path):
        print(f"  PDB file not found: {pdb_path}")
        continue
    
    basename = os.path.basename(pdb_path)
    stem = os.path.splitext(basename)[0]
    protein_name = stem.replace("_trimmed", "")
    
    pred_csv = predictions_dir / f"{stem}.pdb_predictions.csv"
    
    if not pred_csv.exists():
        print(f"  No prediction CSV found: {pred_csv.name}")
        continue

    try:
        df = pd.read_csv(pred_csv)
        df.columns = df.columns.str.strip()
    except Exception as e:
        print(f"  Error reading predictions: {e}")
        continue

    if df.empty:
        print(f"  No pockets predicted")
        continue

    # Take only the first pocket (rank 0)
    pocket = df.iloc[0]

    try:
        center = np.array([pocket['center_x'], pocket['center_y'], pocket['center_z']])
    except KeyError:
        print(f"  Center not found")
        continue

    probability = pocket.get('probability', 0)
    if probability < 0.5:
        print(f"  Pocket probability too low: {probability}")
        continue

    surf_data = pocket.get('surf_atom_ids', '')
    atom_ids = [int(x) for x in surf_data.split() if x.strip()] if isinstance(surf_data, str) else []

    if not atom_ids:
        print(f"  No atom IDs")
        continue

    # Get pocket residues
    residues = pocket.get('residue_ids', '').strip().replace("A_", "").split(' ')
    residues_in_pocket = [int(r) for r in residues if r.isdigit()]
    residues_str = ",".join(str(r) for r in residues_in_pocket) if residues_in_pocket else "NA"
    pocket_start_end = f"{min(residues_in_pocket)}_{max(residues_in_pocket)}" if residues_in_pocket else "NA"

    # Get N-terminal trimmed residues
    n_term_matches = df_trimming[df_trimming['structure'] == protein_name]
    N_residues = int(n_term_matches['n_trimmed'].iloc[0]) if len(n_term_matches) > 0 else 0

    try:
        structure = parser.get_structure(protein_name, pdb_path)
        serial_to_coord = {atom.get_serial_number(): atom.get_coord() 
                         for atom in structure.get_atoms()}
        
        coords = [serial_to_coord[aid] for aid in atom_ids if aid in serial_to_coord]
        if not coords:
            print(f"  No coordinates found")
            continue
        
        coords = np.array(coords)
        furthest_idx = np.argmax(np.linalg.norm(coords - center, axis=1))
        furthest_coord = coords[furthest_idx]
        box_size = 2 * (furthest_coord - center) + 5.0
        box_size = np.maximum(box_size, 22.5)
        
        results.append({
            "protein": protein_name,
            "pocket_score": pocket.get('score', 'N/A'),
            "probability": round(probability, 3),
            "pocket_start_end": pocket_start_end,
            "residues": residues_str,
            "center_x": round(center[0], 3),
            "center_y": round(center[1], 3),
            "center_z": round(center[2], 3),
            "size_x": round(box_size[0], 3),
            "size_y": round(box_size[1], 3),
            "size_z": round(box_size[2], 3),
            "trimmed_n": N_residues,
            "num_atoms": len(atom_ids),
            "pdb_path": pdb_path
        })
        print(f"  Added: {protein_name} | residues {pocket_start_end} | prob={probability:.3f}")
        
    except Exception as e:
        print(f"  Error: {e}")
        continue

# Save results
print(f"\n{'='*60}")
print(f"RESULTS")
print(f"{'='*60}")
print(f"Total pockets (top 1 per protein): {len(results)}")

if results:
    df_results = pd.DataFrame(results)
    df_results.to_csv(output_file, sep='\t', index=False)
    print(f" Saved to: {output_file}")
    print(f"\n Pockets:")
    print(df_results[['protein', 'pocket_start_end', 'probability']].to_string(index=False))
else:
    print(" No results generated!")
    pd.DataFrame(columns=["protein", "pocket_score", "probability", "pocket_start_end", "residues", "center_x", "center_y", "center_z", "size_x", "size_y", "size_z", "trimmed_n", "num_atoms", "pdb_path"]).to_csv(output_file, sep='\t', index=False)
    print(f"Created empty file: {output_file}")

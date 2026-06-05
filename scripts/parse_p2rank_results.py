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
output_file = os.path.join(BASE_DIR, "output_files/docking_boxes.csv")

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

# Process each PDB file
for i, pdb_path in enumerate(pdb_paths):
    print(f"\n[{i+1}/{len(pdb_paths)}] Processing: {os.path.basename(pdb_path)}")
    
    if not os.path.exists(pdb_path):
        print(f"  PDB file not found: {pdb_path}")
        continue
    
    # Extract protein name
    basename = os.path.basename(pdb_path)
    stem = os.path.splitext(basename)[0]
    protein_name = stem.replace("_trimmed", "")
    
    # Look for prediction CSV
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

    # Get N-terminal trimmed residues
    n_term_matches = df_trimming[df_trimming['structure'] == protein_name]
    if len(n_term_matches) > 0:
        N_residues = int(n_term_matches['n_trimmed'].iloc[0])
    else:
        N_residues = 0

    # Take only the top pocket (rank 1)
    pocket = df.iloc[0]

    # Get center
    try:
        center = np.array([pocket['center_x'], pocket['center_y'], pocket['center_z']])
    except KeyError:
        print(f"  Center not found")
        continue

    # Get surface atom IDs
    if 'surf_atom_ids' not in pocket:
        print(f"  No surf_atom_ids")
        continue
    
    surf_data = pocket['surf_atom_ids']
    
    # Parse atom IDs
    atom_ids = []
    if isinstance(surf_data, str):
        atom_ids = [int(x) for x in surf_data.split() if x.strip()]
    elif isinstance(surf_data, (int, np.integer)):
        atom_ids = [int(surf_data)]
    
    if not atom_ids:
        print(f"  No atom IDs")
        continue
    
    # Load PDB and get coordinates
    try:
        structure = parser.get_structure(protein_name, pdb_path)
        serial_to_coord = {atom.get_serial_number(): atom.get_coord() 
                         for atom in structure.get_atoms()}
        
        coords = []
        for atom_id in atom_ids:
            if atom_id in serial_to_coord:
                coords.append(serial_to_coord[atom_id])
        
        if not coords:
            print(f"  No coordinates found")
            continue
            
        coords = np.array(coords)
        
        # Calculate box size
        furthest_idx = np.argmax(np.linalg.norm(coords - center, axis=1))
        furthest_coord = coords[furthest_idx]
        box_size = 2 * (furthest_coord - center) + 5.0
        box_size = np.maximum(box_size, 22.5)
        
        results.append({
            "protein": protein_name,
            "center_x": center[0],
            "center_y": center[1],
            "center_z": center[2],
            "size_x": box_size[0],
            "size_y": box_size[1],
            "size_z": box_size[2],
            "pdb_path": pdb_path
        })
        print(f"  Added docking box for {protein_name}")
        
    except Exception as e:
        print(f"  Error processing: {e}")
        continue

# Save results
print(f"\n{'='*60}")
print(f"RESULTS")
print(f"{'='*60}")
print(f"Total docking boxes: {len(results)}")

if results:
    df_results = pd.DataFrame(results)
    df_results.to_csv(output_file, index=False)
    print(f"\n Saved to: {output_file}")
else:
    print(" No results generated!")
    pd.DataFrame(columns=["protein", "center_x", "center_y", "center_z", "size_x", "size_y", "size_z"]).to_csv(output_file, index=False)
    print(f"Created empty file: {output_file}")

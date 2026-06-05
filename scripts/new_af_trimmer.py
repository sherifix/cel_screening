import os
import pandas as pd
from Bio import PDB

# =========================
# USER PATHS
# =========================
BASE_INPUT_DIR = "/home/ubuntu/storage/thesis/data/hits_3D_structures"
BASE_OUTPUT_DIR = "/home/ubuntu/storage/thesis/data/hits_3D_trimmed_structures"
PDB_LIST_FILE = "/home/ubuntu/storage/thesis/results/pdb_list.ds"  # File with absolute paths

PLDDT_THRESHOLD = 50.0
MIN_TERMINAL_LENGTH = 3

os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

parser = PDB.PDBParser(QUIET=True)
io = PDB.PDBIO()

summary = []

def residue_plddt(residue):
    atoms = list(residue.get_atoms())
    return sum(a.get_bfactor() for a in atoms) / len(atoms) if atoms else 0.0

# Read the list of PDB files to process
print(f"Reading PDB list from {PDB_LIST_FILE}...")
with open(PDB_LIST_FILE, 'r') as f:
    pdb_paths = [line.strip() for line in f if line.strip() and line.strip().endswith('.pdb')]

print(f"Found {len(pdb_paths)} PDB files to process")

# Group PDB paths by GH family for organized output
for pdb_path in pdb_paths:
    # Extract GH family and PDB ID from the absolute path
    # Expected path format: .../thermophiles_structures/GH_family/xxxx.pdb
    path_parts = pdb_path.split(os.sep)
    
    # Find the GH family (should be the directory after thermophiles_structures)
    try:
        pdb_file = path_parts[-1]
        pdb_id = pdb_file.replace('.pdb', '')
    except (ValueError, IndexError):
        print(f"Warning: Could not parse GH family from path: {pdb_path}")
        continue

    if not os.path.exists(pdb_path):
        print(f"Warning: PDB file not found: {pdb_path}")
        continue

    print(f"Processing {pdb_file}...")
    

    try:
        structure = parser.get_structure(pdb_id, pdb_path)
        model = structure[0]
        chain = list(model.get_chains())[0]

        residues = [r for r in chain if PDB.is_aa(r, standard=True)]
        total_residues = len(residues)

        if total_residues < 20:
            print(f"  Skipping very small structure: {pdb_id} ({total_residues} residues)")
            summary.append({
                "structure": pdb_id,
                "total_aa_residues": total_residues,
                "n_terminal_trimmed": 0,
                "c_terminal_trimmed": 0,
                "remaining_aa_residue_count": total_residues,
                "n_terminal_range": "None",
                "c_terminal_range": "None",
                "internal_low_confidence": None,
                "output_path": "SKIPPED - too small",
                "status": "skipped_too_small"
            })
            continue

        res_nums = [r.get_id()[1] for r in residues]
        plddts = [residue_plddt(r) for r in residues]

        # ---- N-terminal trimming ----
        n_trim = 0
        for p in plddts:
            if p < PLDDT_THRESHOLD:
                n_trim += 1
            else:
                break
        
        # Only trim if we meet the minimum length requirement
        # AND we haven't reached the end of the protein
        if n_trim >= total_residues or n_trim < MIN_TERMINAL_LENGTH:
            n_trim = 0

        # ---- C-terminal trimming ----
        c_trim = 0
        for p in reversed(plddts):
            if p < PLDDT_THRESHOLD:
                c_trim += 1
            else:
                break
        
        # Only trim if we meet the minimum length requirement
        # AND we haven't reached the beginning of the protein
        if c_trim >= total_residues or c_trim < MIN_TERMINAL_LENGTH:
            c_trim = 0

        # Ensure we don't overlap trimming
        if n_trim + c_trim > total_residues:
            c_trim = total_residues - n_trim
            if c_trim < 0:
                c_trim = 0

        remaining_residues = total_residues - n_trim - c_trim

        # ---- Check for internal low-confidence regions (for analysis only) ----
        internal_low_conf = False
        if remaining_residues > 0:
            internal_region = plddts[n_trim:total_residues - c_trim]
            internal_low_conf = any(p < PLDDT_THRESHOLD for p in internal_region)

        # ---- Residues to keep ----
        keep_indices = list(range(n_trim, total_residues - c_trim))
        keep_res = set(res_nums[i] for i in keep_indices)

        class TrimSelect(PDB.Select):
            def accept_residue(self, residue):
                return residue.get_id()[1] in keep_res

        out_pdb = os.path.join(BASE_OUTPUT_DIR, f"{pdb_id}_trimmed.pdb")
        io.set_structure(structure)
        io.save(out_pdb, TrimSelect())

        # Safely create range strings
        n_term_range = "None"
        if n_trim > 0 and n_trim <= len(res_nums):
            n_term_range = f"{res_nums[0]}-{res_nums[min(n_trim-1, len(res_nums)-1)]}"
        
        c_term_range = "None"
        if c_trim > 0 and c_trim <= len(res_nums):
            c_term_range = f"{res_nums[max(len(res_nums)-c_trim, 0)]}-{res_nums[-1]}"

        summary.append({
            "structure": pdb_id,
            "total_aa_residues": total_residues,
            "n_terminal_trimmed": n_trim,
            "c_terminal_trimmed": c_trim,
            "remaining_aa_residue_count": remaining_residues,
            "n_terminal_range": n_term_range,
            "c_terminal_range": c_term_range,
            "internal_low_confidence": internal_low_conf,
            "output_path": out_pdb,
            "status": "processed"
        })

        print(f"  Trimmed: N={n_trim}, C={c_trim}, Remaining={remaining_residues}")

    except Exception as e:
        print(f"  Error processing {pdb_path}: {str(e)}")
        summary.append({
            "structure": pdb_id,
            "total_aa_residues": None,
            "n_terminal_trimmed": 0,
            "c_terminal_trimmed": 0,
            "remaining_aa_residue_count": None,
            "n_terminal_range": "None",
            "c_terminal_range": "None",
            "internal_low_confidence": None,
            "output_path": "ERROR",
            "status": f"error: {str(e)}"
        })

# ---- Save global summary ----
df = pd.DataFrame(summary)
summary_path = os.path.join(BASE_OUTPUT_DIR, "trimming_summary.csv")
df.to_csv(summary_path, index=False)

# Print summary statistics
processed = len(df[df['status'] == 'processed'])
skipped = len(df[df['status'] == 'skipped_too_small'])
errors = len(df[df['status'].str.startswith('error', na=False)])

print("\n" + "="*50)
print("PROCESSING COMPLETE")
print("="*50)
print(f"Total files in list: {len(pdb_paths)}")
print(f"Successfully processed: {processed}")
print(f"Skipped (too small): {skipped}")
print(f"Errors: {errors}")
print(f"Summary written to: {summary_path}")

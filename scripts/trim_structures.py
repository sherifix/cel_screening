import os
import pandas as pd
from Bio import PDB
from pathlib import Path

BASE_DIR = os.getcwd()

# Paths
input_dir = Path(BASE_DIR) / "data/query_structures"
output_dir = Path(BASE_DIR) / "data/trimmed_structures"
output_dir.mkdir(parents=True, exist_ok=True)

PLDDT_THRESHOLD = 50.0
MIN_TERMINAL_LENGTH = 3

parser = PDB.PDBParser(QUIET=True)
io = PDB.PDBIO()

summary = []

def residue_plddt(residue):
    atoms = list(residue.get_atoms())
    return sum(a.get_bfactor() for a in atoms) / len(atoms) if atoms else 0.0

# Get all PDB files from query_structures
pdb_files = list(input_dir.glob("*.pdb"))
print(f"Found {len(pdb_files)} structures to trim")

for pdb_path in pdb_files:
    pdb_id = pdb_path.stem
    print(f"Processing {pdb_id}...")
    
    try:
        structure = parser.get_structure(pdb_id, str(pdb_path))
        model = structure[0]
        chain = list(model.get_chains())[0]
        
        residues = [r for r in chain if PDB.is_aa(r, standard=True)]
        total = len(residues)
        
        if total < 20:
            print(f"  Skipping: too small ({total} residues)")
            summary.append({
                "structure": pdb_id,
                "total_residues": total,
                "status": "skipped_too_small"
            })
            continue
        
        # Get residue numbers and pLDDT scores
        res_nums = [r.get_id()[1] for r in residues]
        plddts = [residue_plddt(r) for r in residues]
        
        # N-terminal trim
        n_trim = 0
        for p in plddts:
            if p < PLDDT_THRESHOLD:
                n_trim += 1
            else:
                break
        if n_trim < MIN_TERMINAL_LENGTH or n_trim >= total:
            n_trim = 0
        
        # C-terminal trim
        c_trim = 0
        for p in reversed(plddts):
            if p < PLDDT_THRESHOLD:
                c_trim += 1
            else:
                break
        if c_trim < MIN_TERMINAL_LENGTH or c_trim >= total:
            c_trim = 0
        
        # Ensure we don't overlap
        if n_trim + c_trim > total:
            c_trim = total - n_trim
            if c_trim < 0:
                c_trim = 0
        
        remaining = total - n_trim - c_trim
        
        # Check for internal low-confidence regions (pLDDT < 50 in the middle)
        internal_low_conf = False
        if remaining > 0:
            internal_region = plddts[n_trim:total - c_trim]
            internal_low_conf = any(p < PLDDT_THRESHOLD for p in internal_region)
            if internal_low_conf:
                print(f"  ==== Internal low-confidence region detected (pLDDT < {PLDDT_THRESHOLD})")
        
        # Select residues to keep
        keep_indices = list(range(n_trim, total - c_trim))
        keep_res = set(res_nums[i] for i in keep_indices)
        
        class TrimSelect(PDB.Select):
            def accept_residue(self, residue):
                return residue.get_id()[1] in keep_res
        
        out_file = output_dir / f"{pdb_id}_trimmed.pdb"
        io.set_structure(structure)
        io.save(str(out_file), TrimSelect())
        
        # Create range strings for summary
        n_term_range = "None"
        if n_trim > 0 and n_trim <= len(res_nums):
            n_term_range = f"{res_nums[0]}-{res_nums[min(n_trim-1, len(res_nums)-1)]}"
        
        c_term_range = "None"
        if c_trim > 0 and c_trim <= len(res_nums):
            c_term_range = f"{res_nums[max(len(res_nums)-c_trim, 0)]}-{res_nums[-1]}"
        
        summary.append({
            "structure": pdb_id,
            "total_residues": total,
            "n_trimmed": n_trim,
            "c_trimmed": c_trim,
            "remaining_residues": remaining,
            "n_terminal_range": n_term_range,
            "c_terminal_range": c_term_range,
            "internal_low_confidence": internal_low_conf,
            "output": str(out_file),
            "status": "processed"
        })
        print(f"  Trimmed: N={n_trim}, C={c_trim}, Remaining={remaining}, Internal_low_conf={internal_low_conf}")
        
    except Exception as e:
        print(f"  Error: {e}")
        summary.append({"structure": pdb_id, "status": f"error: {e}"})

# Save summary
df = pd.DataFrame(summary)
output_summary = output_dir / "trimming_summary.csv"
df.to_csv(output_summary, index=False)

print(f"\n Trimmed structures saved to: {output_dir}")
print(f"Summary saved to: {output_summary}")

# Print warning for structures with internal low confidence
internal_conf = df[df.get('internal_low_confidence', False) == True]
if len(internal_conf) > 0:
    print(f"\n WARNING: {len(internal_conf)} structures have internal low-confidence regions:")
    for _, row in internal_conf.iterrows():
        print(f"  - {row['structure']}")
    print("\n  These may affect docking results. Consider manual inspection.")

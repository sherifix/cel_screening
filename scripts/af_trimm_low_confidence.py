import os
import pandas as pd
from Bio import PDB

# =========================
# USER PATHS
# =========================
BASE_INPUT_DIR = "/home/ubuntu/storage/second_run/thesis/data/structure_hits/thermophiles_structures"
BASE_OUTPUT_DIR = "/home/ubuntu/storage/second_run/thesis/data/structure_hits/thermophiles_structures_trimmed"

PLDDT_THRESHOLD = 50.0
MIN_TERMINAL_LENGTH = 3

os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

parser = PDB.PDBParser(QUIET=True)
io = PDB.PDBIO()

summary = []

def residue_plddt(residue):
    atoms = list(residue.get_atoms())
    return sum(a.get_bfactor() for a in atoms) / len(atoms) if atoms else 0.0

for gh_family in sorted(os.listdir(BASE_INPUT_DIR)):
    gh_path = os.path.join(BASE_INPUT_DIR, gh_family)
    if not os.path.isdir(gh_path):
        continue

    print(f"Processing {gh_family}...")
    out_family_dir = os.path.join(BASE_OUTPUT_DIR, gh_family)
    os.makedirs(out_family_dir, exist_ok=True)

    for pdb_file in os.listdir(gh_path):
        if not pdb_file.endswith(".pdb"):
            continue

        pdb_id = pdb_file.replace(".pdb", "")
        pdb_path = os.path.join(gh_path, pdb_file)

        structure = parser.get_structure(pdb_id, pdb_path)
        model = structure[0]
        chain = list(model.get_chains())[0]

        residues = [r for r in chain if PDB.is_aa(r, standard=True)]
        total_residues = len(residues)

        if total_residues < 20:
            print(f"Skipping very small structure: {pdb_id}")
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
        if n_trim < MIN_TERMINAL_LENGTH:
            n_trim = 0

        # ---- C-terminal trimming ----
        c_trim = 0
        for p in reversed(plddts):
            if p < PLDDT_THRESHOLD:
                c_trim += 1
            else:
                break
        if c_trim < MIN_TERMINAL_LENGTH:
            c_trim = 0

        # ---- Internal low-confidence detection ----
        internal_low_conf = any(
            p < PLDDT_THRESHOLD
            for p in plddts[n_trim:total_residues - c_trim]
        )

        remaining_residues = total_residues - n_trim - c_trim

        # ---- Residues to keep ----
        keep_res = set(res_nums[n_trim:total_residues - c_trim])

        class TrimSelect(PDB.Select):
            def accept_residue(self, residue):
                return residue.get_id()[1] in keep_res

        out_pdb = os.path.join(out_family_dir, f"{pdb_id}_trimmed.pdb")
        io.set_structure(structure)
        io.save(out_pdb, TrimSelect())

        summary.append({
            "gh_family": gh_family,
            "structure": pdb_id,
            "total_aa_residues": total_residues,
            "n_terminal_trimmed": n_trim,
            "c_terminal_trimmed": c_trim,
            "remaining_aa_residue_count": remaining_residues,
            "n_terminal_range": f"{res_nums[0]}-{res_nums[n_trim-1]}" if n_trim > 0 else "None",
            "c_terminal_range": f"{res_nums[-c_trim]}-{res_nums[-1]}" if c_trim > 0 else "None",
            "internal_low_confidence": internal_low_conf,
            "output_path": out_pdb
        })

# ---- Save global summary ----
df = pd.DataFrame(summary)
summary_path = os.path.join(BASE_OUTPUT_DIR, "trimming_summary.csv")
df.to_csv(summary_path, index=False)

print("\n✔ All GH families processed")
print(f"✔ Summary written to {summary_path}")

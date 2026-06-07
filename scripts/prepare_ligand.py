import os
import sys
import subprocess
from pathlib import Path

BASE_DIR = os.getcwd()

# Output directory for docking
output_dir = Path(BASE_DIR) / "results/autodock_vina"
output_dir.mkdir(parents=True, exist_ok=True)

# Find the SDF file in results/autodock_vina/
sdf_files = list(output_dir.glob("*.sdf"))

if len(sdf_files) == 0:
    print(f"ERROR: No SDF file found in {output_dir}")
    print("\nPlease place your ligand file (SDF format) in:")
    print(f"  {output_dir}/ (e.g., {output_dir}/cellotetraose.sdf)")
    sys.exit(1)
elif len(sdf_files) > 1:
    print(f"WARNING: Multiple SDF files found in {output_dir}")
    print("Using the first one:", sdf_files[0].name)

input_sdf = sdf_files[0]

# Prepare PDBQT file
output_pdbqt = output_dir / f"{input_sdf.stem}.pdbqt"
print(f"Processing ligand from: {input_sdf}")
print(f"Output PDBQT will be saved to: {output_pdbqt}")

# Run mk_prepare_ligand.py
cmd = ["mk_prepare_ligand.py", "-i", str(input_sdf), "-o", str(output_pdbqt)]

try:
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    print(f" Ligand prepared and saved to: {output_pdbqt}")
except subprocess.CalledProcessError as e:
    print(f"ERROR: mk_prepare_ligand.py failed")
    print(f"STDERR: {e.stderr}")
    sys.exit(1)

# Create completion flag
flag_file = Path(BASE_DIR) / "data/raw/.ligand_prepared_complete"
flag_file.touch()
print(f" Completion flag created: {flag_file}")

print(f"\n Docking directory: {output_dir}")
print(f"   - Original SDF: {input_sdf.name}")
print(f"   - Prepared PDBQT: {output_pdbqt.name}")

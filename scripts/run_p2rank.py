import subprocess
import os
import sys
from pathlib import Path

BASE_DIR = os.getcwd()

# Paths
p2rank_bin = os.path.expanduser("~/tools/p2rank_2.5.1/prank")
pdb_list = os.path.join(BASE_DIR, "data/trimmed_pdb_list.ds")
output_dir = os.path.join(BASE_DIR, "results/p2rank_predictions")

# Check if prank exists
if not os.path.exists(p2rank_bin):
    print(f"ERROR: P2Rank not found at {p2rank_bin}")
    print("Please install P2Rank or update the path")
    sys.exit(1)

# Check if pdb_list exists
if not os.path.exists(pdb_list):
    print(f"ERROR: PDB list not found at {pdb_list}")
    print("Run prepare_pdb_lists.py first")
    sys.exit(1)

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Run P2Rank
print(f"Running P2Rank on {pdb_list}...")
print(f"Command: {p2rank_bin} predict -c alphafold -o {output_dir} {pdb_list}")

cmd = [p2rank_bin, "predict", "-c", "alphafold", "-o", output_dir, pdb_list]

try:
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0:
        print(" P2Rank completed successfully")
        print(f"Output saved to: {output_dir}")
        # List output files
        print("\n Generated pocket prediction files:")
        for f in Path(output_dir).glob("*.pdb"):
            print(f"  - {f.name}")
        for f in Path(output_dir).glob("*.csv"):
            print(f"  - {f.name}")
    else:
        print(f" P2Rank failed with error:")
        print(result.stderr)
        sys.exit(1)
except Exception as e:
    print(f" Error running P2Rank: {e}")
    sys.exit(1)

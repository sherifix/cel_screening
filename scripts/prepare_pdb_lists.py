import pandas as pd
from pathlib import Path
import os

BASE_DIR = os.getcwd()

# Read top hits
top_hits_file = Path(BASE_DIR) / "output_files/top_hits.csv"
df_hits = pd.read_csv(top_hits_file)

# Get candidate accessions
candidates = set(df_hits['accession'].str.upper())
print(f"Loaded {len(candidates)} top hits")

# Create data directory for .ds files
data_dir = Path(BASE_DIR) / "data"
data_dir.mkdir(exist_ok=True)

#Create list for untrimmed structures
untrimmed_dir = Path(BASE_DIR) / "data/query_structures"
untrimmed_list = data_dir / "pdb_list.ds"

pdb_paths = []
for pdb_file in untrimmed_dir.glob("*.pdb"):
    stem = pdb_file.stem.upper()
    if stem in candidates:
        pdb_paths.append(str(pdb_file.resolve()))

with open(untrimmed_list, 'w') as f:
    for path in pdb_paths:
        f.write(f"{path}\n")

print(f" Untrimmed list: {len(pdb_paths)} PDBs → {untrimmed_list}")


# Create list for trimmed structures
trimmed_dir = Path(BASE_DIR) / "data/trimmed_structures"
trimmed_list = data_dir / "trimmed_pdb_list.ds"

if trimmed_dir.exists():
    trimmed_paths = []
    for pdb_file in trimmed_dir.glob("*_trimmed.pdb"):
        stem = pdb_file.stem.replace("_trimmed", "").upper()
        if stem in candidates:
            trimmed_paths.append(str(pdb_file.resolve()))
    
    with open(trimmed_list, 'w') as f:
        for path in trimmed_paths:
            f.write(f"{path}\n")
    
    print(f" Trimmed list: {len(trimmed_paths)} PDBs → {trimmed_list}")
else:
    print(f" trimmed structures directory not found: {trimmed_dir}")

# Summary
print(f"\n Summary:")
print(f"  Total top hits: {len(candidates)}")
print(f"  Untrimmed structures: {len(pdb_paths)}")
print(f"  Trimmed structures: {len(trimmed_paths) if trimmed_dir.exists() else 0}")

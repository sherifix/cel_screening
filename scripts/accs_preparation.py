import pandas as pd
import os

GH_FAMILIES = [
    'GH5_1', 'GH5_5', 'GH5_2', 'GH5_4', 'GH5_9',
    'GH5_25', 'GH5_26', 'GH5_37', 'GH5_38', 'GH45',
    'GH6', 'GH7', 'GH9', 'GH44', 'GH12', 'GH48'
]

os.chdir("./data/raw")

for gh_fam in GH_FAMILIES:
    csv_file = f'{gh_fam}.csv'
    
    # Skip if file is empty
    if os.path.getsize(csv_file) == 0:
        print(f"Skipping {csv_file}: empty file")
        continue
    
    df = pd.read_csv(csv_file)
    
    if df["GenBank_accession"].empty:
        print(f"Skipping {csv_file}: no accessions")
        continue
    
    df["GenBank_accession"].to_csv(
        f"{gh_fam}.txt",
        index=False,
        header=False
    )
import requests
import pandas as pd
import os
from pathlib import Path
import time
import sys

BASE_DIR = os.getcwd()

# Use hits_thermo_sp.csv
input_csv = os.path.join(BASE_DIR, "output_files/hits_thermo_sp.csv")
df_hits = pd.read_csv(input_csv)
print(f"Loaded {len(df_hits)} hits from {input_csv}")

# Get unique accessions
accessions = df_hits['accession'].dropna().unique()
print(f"Processing {len(accessions)} unique accessions")

# Create output directories
struct_dir = os.path.join(BASE_DIR, "data/query_structures")
os.makedirs(struct_dir, exist_ok=True)

def get_uniprot(genbank_id):
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {"query": genbank_id, "format": "json", "size": 1}
    try:
        r = requests.get(url, params=params, timeout=30)
        if r.status_code == 200 and r.json()["results"]:
            return r.json()["results"][0]["primaryAccession"]
    except:
        pass
    return None

results = []
missing_list = []

for i, acc in enumerate(accessions):
    print(f"Processing {i+1}/{len(accessions)}: {acc}")
    uniprot = get_uniprot(acc)
    results.append({"genbank": acc, "uniprot": uniprot})
    
    if uniprot:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v6.pdb"
        out_file = os.path.join(struct_dir, f"{acc}.pdb")
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 200:
                with open(out_file, "wb") as f:
                    f.write(r.content)
                print(f"  === Downloaded: {uniprot}")
            else:
                print(f"  --> No structure for {uniprot}")
                missing_list.append(acc)
        except Exception as e:
            print(f"  >>> Error: {e}")
            missing_list.append(acc)
    else:
        print(f"  --> No Uniprot ID for {acc}")
        missing_list.append(acc)
    
    time.sleep(0.3)

# Save mapping
df_map = pd.DataFrame(results)
df_map.to_csv(os.path.join(BASE_DIR, "output_files/genbank_uniprot_map.csv"), index=False)

# Report missing structures
if missing_list:
    print("\n" + "="*60)
    print(f"===== MISSING STRUCTURES: {len(missing_list)} out of {len(accessions)}")
    print("="*60)
    for acc in missing_list:
        print(f"  - {acc}")
    print("="*60)
    
    # Ask user what to do
    while True:
        answer = input("\nDo you want to:\n  [1] Continue with available structures only\n  [2] Stop and add missing structures manually, then rerun\n\nEnter Continue or Stop: ")
        
        if answer.upper() == "CONTINUE" or answer == "1":
            print("\n===== Continuing with available structures only...")
            print(f"  Available: {len(accessions) - len(missing_list)} structures")
            print(f"  Missing: {len(missing_list)} structures (will be skipped in alignment)")
            break
        elif answer.upper() == "STOP" or answer == "2":
            print("\n===== Stopping pipeline. To continue:")
            print("   1. Predict structures for missing accessions using a 3D structure prediction tool")
            print("   2. Save as PDB files in: data/query_structures/")
            print("   3. Use filename: {accession}.pdb (e.g., AAA23226.1.pdb)")
            print("   4. Re-run: snakemake --cores 8")
            sys.exit(1)
        else:
            print("Invalid input. Please enter Continue or Stop.")
else:
    print(f"\n===== All {len(accessions)} structures downloaded successfully!")

print(f"\nStructures saved to: {struct_dir}")
print(f"Mapping saved to: output_files/genbank_uniprot_map.csv")

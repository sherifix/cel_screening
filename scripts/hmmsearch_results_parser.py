import os
import pandas as pd
from Bio import SearchIO

# --- Paths ---
hmmer_results = "/home/ubuntu/storage/thesis/test/results/hmmsearch_results"

blacklist_file = "/home/ubuntu/storage/thesis/data/all_characterized_raw_data/all_genbank_accessions.txt"

# --- Load blacklist accessions ---
with open(blacklist_file) as f:
    remove_accessions = set(line.strip() for line in f)

# --- Initialize storage ---
rows = []

# --- Walk through all .tbl files ---
for root, dirs, files in os.walk(hmmer_results):
    for file in files:
        if file.endswith(".domtbl"):
            tbl_path = os.path.join(root, file)
            
            # Extract taxID from filename (before first "_")
            taxid = file.split("_")[0]

            # Parse hmmscan --tblout
            for qresult in SearchIO.parse(tbl_path, "hmmsearch3-domtab"):
                family = qresult.id.split("_trimmed")[0]

                # Skip blacklisted accessions
                

                for hit in qresult:
                    for hsp in hit:
                        accession  = hsp.hit_id
                        if accession in remove_accessions:
                            print(f"{accession} deleted")
                            continue
                        else:    
                            rows.append({"taxID": taxid,
                                "accession": accession,
                                "family": family,
                                "evalue": hsp.evalue,                                         
                                "evalue_cond": hsp.evalue_cond,
                                "score": hsp.bitscore,
                                 "domain_coord": f"{hsp.env_start}:{hsp.env_end}",
                                 "domain_len": hsp.env_end-hsp.env_start
                                })

# --- Build final DataFrame ---
df = pd.DataFrame(rows)
df = df[(df['evalue'] <= 1e-5) &
    (df['evalue_cond'] <= 1e-5)
    ]

# --- Quick summary ---
print(f"Total hits : {len(df)}")
print(df.head())

df.style.format({"evalue": "{:.2e}"})

output_csv = "/home/ubuntu/storage/thesis/test/results/raw_results_hmmsearch.csv"
df.to_csv(output_csv, index=False)
print(f"Saved {len(df)} hits to {output_csv}")

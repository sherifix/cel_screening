import os
import pandas as pd
from Bio import SearchIO


hmmer_results = "results/hmmsearch_results"
output_csv = "results/raw_results_hmmsearch.csv"
blacklist_file = "data/blacklist_accessions.txt"


# Load blacklist accessions if file exists
remove_accessions = set()
if os.path.exists(blacklist_file):
    with open(blacklist_file) as f:
        remove_accessions = set(line.strip() for line in f)
    print(f"Loaded {len(remove_accessions)} blacklisted accessions")
else:
    print(f"Warning: {blacklist_file} not found. No accessions will be blacklisted.")


rows = []

#Walk through all .domtbl files
for root, dirs, files in os.walk(hmmer_results):
    for file in files:
        if file.endswith(".domtbl"):
            tbl_path = os.path.join(root, file)

            # Extract proteome name from filename
            proteome = file.replace(".domtbl", "")
            # Remove family suffix to get proteome name
            parts = proteome.split("_")
            if len(parts) > 1:
                proteome_name = "_".join(parts[:-1])
            else:
                proteome_name = proteome

            try:
                # Parse hmmscan --domtblout
                for qresult in SearchIO.parse(tbl_path, "hmmsearch3-domtab"):
                    family = qresult.id.replace("_trimmed", "")

                    for hit in qresult:
                        for hsp in hit:
                            accession = hit.id

                            # Skip blacklisted accessions
                            if accession in remove_accessions:
                                print(f"Skipping blacklisted: {accession}")
                                continue

                            rows.append({
                                "accession": accession,
                                "family": family,
                                "evalue": hsp.evalue,
                                "evalue_cond": hsp.evalue_cond,
                                "score": hsp.bitscore,
                                "domain_coord": f"{hsp.env_start}:{hsp.env_end}",
                                "domain_len": hsp.env_end - hsp.env_start
                            })
            except Exception as e:
                print(f"Error parsing {tbl_path}: {e}")
                continue


if rows:
    df = pd.DataFrame(rows)
    df = df[(df['evalue'] <= 1e-5) & (df['evalue_cond'] <= 1e-5)]

    # quick summary
    print(f"Total hits after filtering: {len(df)}")

    # Create output directory if needed
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    df.to_csv(output_csv, index=False)
    print(f"Saved {len(df)} hits to {output_csv}")
else:
    print("No hits found!")

import pandas as pd
import re
import subprocess
from Bio import SeqIO
import os
from pathlib import Path


os.makedirs("home/ubuntu/storage/thesis/test/output_files", exist_ok=True)
os.makedirs("/home/ubuntu/storage/thesis/test/data/fasta_files", exist_ok=True)

df = pd.read_csv("../results/raw_results_hmmsearch.csv")
df["evalue"] = pd.to_numeric(df["evalue"])


#filter candidates with overlapping domains
coords = df["domain_coord"].astype(str).str.extract(r"(\d+):(\d+)")
df["start"] = coords[0].astype(int)
df["end"] = coords[1].astype(int)

df = df.sort_values(["accession", "evalue"], ascending=[True, True])

kept_rows = []

for acc, sub in df.groupby("accession", sort=False):
    sub = sub.sort_values("evalue", ascending=True)
    selected = []
    for _, row in sub.iterrows():
        overlap = False
        for s in selected:
            if not (row["end"] < s["start"] or row["start"] > s["end"]):
                overlap = True
                break
        if not overlap:
            selected.append(row)
    kept_rows.extend(selected)

df_filtered = pd.DataFrame(kept_rows).drop(columns=["start", "end"])



# get domain of entries: bacteria, archaea or eukaryote
def get_domain(taxid):
    cmd = subprocess.run(f'efetch -db taxonomy -id {taxid} -mode xml | xtract -pattern LineageEx -block Taxon -if Rank -equals domain -element ScientificName' , shell = True, capture_output=True, 
        text=True, 
        check=True
    )
    return cmd.stdout.strip()

#get domain for only unique taxIDs
taxids = df_filtered["taxID"].dropna().tolist()
unique_taxids = set(taxids)

domain_mapping = {}
for taxid in unique_taxids:
    domain_mapping[f'{int(taxid)}'] = get_domain(int(taxid))

df_filtered["domain"] = df_filtered["taxID"].astype(str).map(domain_mapping)

#save first .csv filtered file
df_filtered.to_csv("../output_files/filtered_hmmsearch_results.csv", index=False)


#prepare fasta file with sequences of the filtered accessions for next steps of analysis
accessions = set(df_filtered["accession"].tolist())

output_fasta = "../data/fasta_files/filtered_hmmsearch.faa"
with open(output_fasta, "w") as f:
    for proteome_dir in Path("../data/thermobase_data/proteomes").glob("*"):
        if proteome_dir.is_dir():
            fasta_file = proteome_dir / f"{proteome_dir.name}.faa"
            if fasta_file.exists():
                print(f"scanning {fasta_file.name}..")
                for record in SeqIO.parse(fasta_file, "fasta"):
                    if record.id in accessions:
                        SeqIO.write(record, f, "fasta")

print(f" Extracted sequences to {output_fasta}")


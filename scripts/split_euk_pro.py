import pandas as pd
from Bio import SeqIO


df = pd.read_csv("./output_files/filtered_hmmsearch_results.csv")


df_prok = df[df['domain'].isin(['Bacteria', 'Archaea'])]
df_euk =df[df['domain'].isin(['Eukaryota'])]


target_pro = set(df_prok['accession'].tolist())
target_euk = set(df_euk['accession'].tolist())


with open("./data/fasta_files/prokaryotic.faa", "w") as out_fh:
    for record in SeqIO.parse("./data/fasta_files/filtered_hmmsearch.faa", "fasta"):
        if record.description.split(' ')[0] in target_pro:
            SeqIO.write(record, out_fh, "fasta")

with open("./data/fasta_files/eukaryotic.faa", "w") as out_fh:
    for record in SeqIO.parse("./data/fasta_files/filtered_hmmsearch.faa", "fasta"):
        if record.description.split(' ')[0] in target_euk:
            SeqIO.write(record, out_fh, "fasta")
            
import pandas as pd
from Bio import SeqIO


df_char = pd.read_csv("/home/ubuntu/storage/thesis/test/data/CAZy_all_GenBank_accessions.txt", header=None)
df_filtered = pd.read_csv("/home/ubuntu/storage/thesis/test/output_files/hits_thermo_sp.csv")


df_filtered['is_characterized'] = df_filtered['accession'].isin(df_char.iloc[:, 0])

df_cryptic = df_filtered[~df_filtered['is_characterized']]

print(f"Cryptic: {len(df_cryptic)} / {len(df_filtered)} total")

df_cryptic['domain'].value_counts()


protein_lengths = {}

for record in SeqIO.parse("/home/ubuntu/storage/thesis/test/data/fasta_files/filtered_hmmsearch.faa", "fasta"):
    header = record.description
    first_token = header.split()[0]
    
    # store by accession-like IDs
    protein_lengths[first_token] = len(record.seq)
    
    # also store by record.id
    protein_lengths[record.id] = len(record.seq)

# map length using accession first, then header-derived ID
df_cryptic['total_aa'] = df_cryptic['accession'].map(protein_lengths)

if 'Header' in df_cryptic.columns:
    df_cryptic['total_aa'] = df_cryptic['total_aa'].fillna(
        df_cryptic['Header'].astype(str).str.split().str[0].map(protein_lengths)
    )

df_cryptic = df_cryptic.drop(columns="is_characterized")

df_cryptic.to_csv("/home/ubuntu/storage/thesis/test/output_files/filtered_cryptic_hits.csv", index=False)

print("Characterized data excluded.\nProtein length Added")
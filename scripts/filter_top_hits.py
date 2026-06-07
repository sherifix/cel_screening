import pandas as pd
import os
from Bio import SeqIO

BASE_DIR = os.getcwd()

# Read the US-align results
usalign_file = os.path.join(BASE_DIR, "output_files/usalign_summary.csv")
df_usalign = pd.read_csv(usalign_file)
print(f"Loaded {len(df_usalign)} US-align results")

# Read the hits file
hits_file = os.path.join(BASE_DIR, "output_files/hits_thermo_sp.csv")
df_hits = pd.read_csv(hits_file)
print(f"Loaded {len(df_hits)} hits from hits_thermo_sp.csv")

# Check if total_aa exists, if not calculate from fasta
if 'total_aa' not in df_hits.columns:
    print("total_aa column not found. Calculating from sequences...")
    fasta_file = os.path.join(BASE_DIR, "data/fasta_files/filtered_hmmsearch.faa")
    if os.path.exists(fasta_file):
        seq_lengths = {}
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_lengths[record.id] = len(record.seq)
        df_hits['total_aa'] = df_hits['accession'].map(seq_lengths)
        print(f"Added total_aa for {df_hits['total_aa'].notna().sum()} sequences")
    else:
        print("Warning: No fasta file found. total_aa will be NaN")

# Merge on family + accession
merged = pd.merge(df_hits, df_usalign, 
                  left_on=['family', 'accession'], 
                  right_on=['family', 'candidate'],
                  how='inner')
print(f"Merged {len(merged)} rows")

# Remove unnecessary columns
columns_to_drop = ['txt_file', 'reference', 'candidate', 'TM_score_candidate', 'Coverage_cand_percent']
merged = merged.drop(columns=[c for c in columns_to_drop if c in merged.columns])

df_filtered = merged

#delete under me 

print(df_filtered)
'''
# Apply filters: TM-score >= 0.75, RMSD <= 2.06
df_filtered = merged[
    (merged['TM_score_ref'] >= 0.75) & 
    (merged['RMSD'] <= 2.06)
]

# Apply total_aa filter if column exists
if 'total_aa' in df_filtered.columns:
    df_filtered = df_filtered[df_filtered['total_aa'] <= 800]
    print(f"After length filter (<=800 aa): {len(df_filtered)}")
'''
# Sort and deduplicate
df_top_hits = (
    df_filtered.sort_values(by="TM_score_ref", ascending=False)
               .drop_duplicates(subset="accession", keep="first")
)

print(f"\n Filtering results:")
print(f"  After TM-score/RMSD filter: {len(df_filtered)}")
print(f"  Top unique hits: {len(df_top_hits)}")

# Save top hits
output_file = os.path.join(BASE_DIR, "output_files/top_hits.csv")
df_top_hits.to_csv(output_file, index=False)
print(f"\n✓ Saved {len(df_top_hits)} top hits to {output_file}")

# Display top hits
if len(df_top_hits) > 0:
    print("\n TOP HITS:")
    cols = ['accession', 'family', 'TM_score_ref', 'RMSD']
    if 'total_aa' in df_top_hits.columns:
        cols.append('total_aa')
    print(df_top_hits[cols].head(10).to_string(index=False))
else:
    print("\n No hits passed the filters. Try adjusting thresholds.")

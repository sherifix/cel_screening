import pandas as pd
import os

BASE_DIR = os.getcwd()

# Read thermoprot predictions and filtered hmmsearch results
df_thermo = pd.read_csv(os.path.join(BASE_DIR, "output_files/thermoprot_prediction.csv"))
df_filtered = pd.read_csv(os.path.join(BASE_DIR, "output_files/filtered_hmmsearch_results.csv"))

# Add accession column and merge both dataframes
df_thermo['accession'] = df_thermo['Header'].str.split(' ').str[0]

df_merged = pd.merge(
    left=df_filtered,
    right=df_thermo,
    how='outer',
    left_on='accession',
    right_on='accession'
)

# Read signal predictions for prokaryotes and eukaryotes
signalp_results_dir = os.path.join(BASE_DIR, "results/signalp")
output_dir = os.path.join(BASE_DIR, "output_files")
os.makedirs(output_dir, exist_ok=True)

df_sp_pro = pd.read_csv(
    os.path.join(signalp_results_dir, "pro/prediction_results.txt"),
    sep='\t', skiprows=1
)
df_sp_euk = pd.read_csv(
    os.path.join(signalp_results_dir, "euk/prediction_results.txt"),
    sep='\t', skiprows=1
)

# Clean column and save a copy
df_sp_euk['CS_Prob'] = df_sp_euk['CS Position'].str.extract(r'Pr:\s*([\d.]+)').astype(float)
df_sp_pro['CS_Prob'] = df_sp_pro['CS Position'].str.extract(r'Pr:\s*([\d.]+)').astype(float)

df_sp_euk.to_csv(os.path.join(output_dir, "signalp_eukaryotic.csv"), index=False)
df_sp_pro.to_csv(os.path.join(output_dir, "signalp_prokaryotic.csv"), index=False)

# Keep only hits with SP prediction
df_sp_euk = df_sp_euk[
    (df_sp_euk['Prediction'] == 'SP') &
    (df_sp_euk['SP(Sec/SPI)'] > 0.55)
]
df_sp_pro = df_sp_pro[
    (df_sp_pro['Prediction'] == 'SP') &
    (df_sp_pro['SP(Sec/SPI)'] > 0.55)
]

# Concatenate both dataframes
df_sp = pd.concat([df_sp_euk, df_sp_pro])
df_sp = df_sp.drop(columns=['LIPO(Sec/SPII)', 'TAT(Tat/SPI)', 'TATLIPO(Tat/SPII)', 'PILIN(Sec/SPIII)'], errors='ignore')

# Merge thermoprot predictions with signalp predictions
df_final = pd.merge(
    left=df_merged,
    right=df_sp,
    how='outer',
    left_on='Header',
    right_on='# ID'
)

# Keep entries which are secreted and thermophilic
df_final = df_final[df_final['Prediction_y'] == 'SP']
df_final = df_final[df_final['Prediction_x'] == 'Thermophilic/Hyperthermophilic']

# Drop unnecessary columns
columns_to_drop = ['Class', '# ID', 'OTHER', 'Prediction_y', 'CS_Prob']
df_final = df_final.drop(columns=[c for c in columns_to_drop if c in df_final.columns], errors='ignore')
df_final = df_final.reset_index(drop=True)

# Save thermophilic, secreted hits
output_csv = os.path.join(output_dir, "hits_thermo_sp.csv")
df_final.to_csv(output_csv, index=False)

print(f"Saved {len(df_final)} thermophilic and secreted hits to {output_csv}")
print(f"Columns saved: {list(df_final.columns)}")

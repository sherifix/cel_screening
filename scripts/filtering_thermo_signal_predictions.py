import pandas as pd


# read prediction of thermoprot and filtered hmmsearch results
df_thermo = pd.read_csv('/home/ubuntu/storage/thesis/test/output_files/thermoprot_prediction.csv')
df_filtered = pd.read_csv('/home/ubuntu/storage/thesis/test/output_files/filtered_hmmsearch_results.csv')


# add accession column and merge both dataframes
df_thermo['accession'] = df_thermo['Header'].str.split(' ').str[0]

df_new = pd.merge(left = df_filtered,
                          right= df_thermo,
                          how= 'outer',
                          left_on = 'accession',
                           right_on='accession')


# read signal prediction for prokaryotes and eukaryotes
df_sp_pro = pd.read_csv('/home/ubuntu/storage/thesis/test/results/signalp/pro/prediction_results.txt', sep = '\t', skiprows=1)
df_sp_euk = pd.read_csv('/home/ubuntu/storage/thesis/test/results/signalp/euk/prediction_results.txt', sep = '\t', skiprows=1)

# clean column and save a copy
df_sp_euk['CS_Prob'] = df_sp_euk['CS Position'].str.extract(r'Pr:\s*([\d.]+)').astype(float)
df_sp_pro['CS_Prob'] = df_sp_pro['CS Position'].str.extract(r'Pr:\s*([\d.]+)').astype(float)

df_sp_euk.to_csv("/home/ubuntu/storage/thesis/test/output_files/signalp_eukaryotic.csv", index=False)
df_sp_pro.to_csv("/home/ubuntu/storage/thesis/test/output_files/signalp_prokaryotic.csv", index = False)


# keep only hits with SP prediction
df_sp_euk = df_sp_euk[
    (df_sp_euk['Prediction'] == 'SP') &
    (df_sp_euk['SP(Sec/SPI)'] > 0.55) 
    ]
df_sp_pro = df_sp_pro[
    (df_sp_pro['Prediction'] == 'SP') &
    (df_sp_pro['SP(Sec/SPI)'] > 0.55)
    ]


# concatenate both dataframes and delete unneeded columns
frames = [df_sp_euk, df_sp_pro]
df_sp = pd.concat(frames)

df_sp = df_sp.drop(columns=['LIPO(Sec/SPII)', 'TAT(Tat/SPI)', 'TATLIPO(Tat/SPII)', 	'PILIN(Sec/SPIII)'])


#merge thermoprot predictions with signalp predictions
df_filtered = pd.merge(left = df_new,
                          right= df_sp,
                          how= 'outer',
                          left_on = 'Header',
                           right_on='# ID')

# keep entries which are secreted and thermophilic
df_filtered = df_filtered[df_filtered['Prediction_y'] == 'SP']
df_filtered = df_filtered[df_filtered['Prediction_x'] == 'Thermophilic/Hyperthermophilic']

# drop some columns
df_filtered = df_filtered.drop(columns=['Class', '# ID', 'OTHER', 'Prediction_y', 'CS_Prob'])

df_filtered = df_filtered.reset_index(drop=True)



# save thermophilic, secreted hits
df_filtered.to_csv('/home/ubuntu/storage/thesis/test/output_files/hits_thermo_sp.csv', index = False)














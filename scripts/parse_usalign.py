#!/usr/bin/env python3

import glob
import pandas as pd
import re
import os

BASE_DIR = os.getcwd()

rows = []

# Look for USalign output files in results/usalign_results
for txt in glob.glob(os.path.join(BASE_DIR, "results/usalign_results/*.txt")):
    filename = os.path.basename(txt)

    # Parse filename: accession_vs_reference.txt
    name = filename.replace(".txt", "")

    if "_vs_" not in name:
        continue

    candidate, reference = name.split("_vs_")

    tm1 = tm2 = rmsd = aln_len = seqid = coverage_cand = coverage_ref = None

    with open(txt) as f:
        for line in f:
            # Aligned length, RMSD, Seq_ID (ALL IN ONE LINE)
            m = re.search(
                r"Aligned length=\s*(\d+),\s*RMSD=\s*([\d\.]+),\s*Seq_ID=.*=\s*([\d\.]+)",
                line
            )
            if m:
                aln_len = int(m.group(1))
                rmsd = float(m.group(2))
                seqid = float(m.group(3))

            # TM-score normalized by structure 1 (candidate)
            m1 = re.search(r"TM-score=\s*([\d\.]+).*Structure_1:\s*L=(\d+)", line)
            if m1:
                tm1 = float(m1.group(1))
                if aln_len:
                    coverage_cand = aln_len / int(m1.group(2)) * 100

            # TM-score normalized by structure 2 (reference)
            m2 = re.search(r"TM-score=\s*([\d\.]+).*Structure_2:\s*L=(\d+)", line)
            if m2:
                tm2 = float(m2.group(1))
                if aln_len:
                    coverage_ref = aln_len / int(m2.group(2)) * 100

    rows.append({
        "candidate": candidate,
        "reference": reference,
        "TM_score_candidate": tm1,
        "TM_score_ref": tm2,
        "Coverage_cand_percent": coverage_cand,
        "Coverage_ref_percent": coverage_ref,
        "RMSD": rmsd,
        "aligned_length": aln_len,
        "seq_identity": seqid,
        "txt_file": txt
    })

if rows:
    df = pd.DataFrame(rows)

    # Add family information from hits_thermo_sp.csv
    hits_file = os.path.join(BASE_DIR, "output_files/hits_thermo_sp.csv")
    if os.path.exists(hits_file):
        df_hits = pd.read_csv(hits_file)
        acc_to_family = dict(zip(df_hits['accession'], df_hits['family']))
        df['family'] = df['candidate'].map(acc_to_family)

    output_file = os.path.join(BASE_DIR, "output_files/usalign_summary.csv")
    df.to_csv(output_file, index=False)
    print(f"[DONE] Parsed {len(df)} alignments → {output_file}")

    # Print summary
    print(f"\n ======== Summary:")
    print(f"  Total alignments: {len(df)}")
    print(f"  Unique candidates: {df['candidate'].nunique()}")
    if 'family' in df.columns:
        print(f"\n  By family:")
        print(df.groupby('family').size())
else:
    print("No alignment files found in results/usalign_results/")
    # Create empty file to avoid snakemake errors
    output_file = os.path.join(BASE_DIR, "results/usalign_summary.csv")
    pd.DataFrame(columns=["candidate", "reference", "TM_score_candidate", "TM_score_ref", "RMSD"]).to_csv(output_file, index=False)
    print(f"Created empty file: {output_file}")

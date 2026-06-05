#!/usr/bin/env python3

import glob
import pandas as pd
import re

rows = []

for txt in glob.glob("/home/ubuntu/storage/thesis/data/hits_3D_structures/**/*.txt", recursive=True):
    parts = txt.split("/")
    family = parts[-2]
    name = parts[-1].replace(".txt", "")

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

            # TM-score normalized by structure 1
            m1 = re.search(r"TM-score=\s*([\d\.]+).*Structure_1:\s*L=(\d+)", line)
            if m1:
                tm1 = float(m1.group(1))
                coverage_cand = aln_len / int(m1.group(2)) * 100
                
            # TM-score normalized by structure 2
            m2 = re.search(r"TM-score=\s*([\d\.]+).*Structure_2:\s*L=(\d+)", line)
            if m2:
                tm2 = float(m2.group(1))
                coverage_ref = aln_len / int(m2.group(2)) * 100


    rows.append({
        "family": family,
        "candidate": candidate,
        "reference": reference,
        "TM_score_candidate": tm1,
        "TM_score_ref": tm2,
        "Coverage_cand_len (%)": coverage_cand,
        "Coverage_ref_len (%)": coverage_ref,
        "RMSD": rmsd,
        "aligned_length": aln_len,
        "seq_identity": seqid,
        "txt_file": txt
    })

df = pd.DataFrame(rows)
df.to_csv("../results/USalign_summary.csv", index=False)

print(f"[DONE] Parsed {len(df)} alignments → USalign_summary.csv")


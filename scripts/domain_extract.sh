#!/bin/bash
set -euo pipefail

THREADS=8
HMM_FILE="../../dbcan/dbCAN.hmm"  

# parent directory containing family/subfamily folders
PARENT_DIR="data/GH_families"

# output directory
OUTPUT_DIR="$PWD/data/extracted_domains"

mkdir -p $OUTPUT_DIR

# Main Loop
cd "$PARENT_DIR"

for fasta in GH*.faa ; do
    family="${fasta%.faa}"
    echo 'processing ${family}'

    hmmscan --domtblout $OUTPUT_DIR/${family}_dom.tbl --cpu $THREADS "$HMM_FILE" "$fasta"
done

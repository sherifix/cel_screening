#!/bin/bash

THREADS=8
HMM_FILE="$1"
PARENT_DIR="data/GH_families"
OUTPUT_DIR="data/extracted_domains"

mkdir -p "$OUTPUT_DIR"

if [[ ! -f "$HMM_FILE" ]]; then
    echo "ERROR: HMM file not found: $HMM_FILE"
    exit 1
fi

for fasta in "$PARENT_DIR"/GH*.faa; do
    if [[ ! -f "$fasta" ]]; then
        echo "No GH*.faa files found in $PARENT_DIR"
        exit 0
    fi

    family=$(basename "$fasta" .faa)
    echo "Processing $family"

    hmmscan --domtblout "$OUTPUT_DIR/${family}_dom.tbl" --cpu "$THREADS" "$HMM_FILE" "$fasta"
done

echo "Domain extraction complete"

#!/bin/bash

# Activate ephod environment
source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate ephod

BASE_DIR="$(pwd)"

FASTA="${BASE_DIR}/data/fasta_files/docked_hits.faa"
OUT_DIR="${BASE_DIR}/output_files"
CSV='opt_ph_prediction.csv'

mkdir -p "$OUT_DIR"

echo "================="
echo "Running EpHod"
echo "================="

# Check if FASTA file exists and is not empty
if [ ! -s "$FASTA" ]; then
    echo "WARNING: FASTA file is empty or missing: $FASTA"
    echo "Creating empty output files"
    touch "${OUT_DIR}/${CSV}"
    exit 0
fi

cd ~/tools/EpHod

python ./ephod/run.py \
    --fasta_path "$FASTA" \
    --save_dir "$OUT_DIR" \
    --csv_name "$CSV" \
    --verbose 1 \
    --save_attention_weights 0 \
    --save_embeddings 0

conda deactivate

echo "================="
echo "EpHod completed"
echo "Output saved to: $OUT_DIR"
echo "================="

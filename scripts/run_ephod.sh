#!/bin/bash

FASTA="/home/ubuntu/storage/thesis/results/docked_hits.faa"
OUT_DIR="/home/ubuntu/storage/thesis/results/ephod"
CSV='opt_ph_prediction.csv'

mkdir -p "$OUT_DIR"


echo "================="
echo "Running EpHod"
echo "================="

cd ~/tools/EpHod

python ./ephod/run.py \
    --fasta_path "$FASTA" \
    --save_dir "$OUT_DIR" \
    --csv_name "$CSV" \
    --verbose 1 \
    --save_attention_weights 0 \
    --save_embeddings 0
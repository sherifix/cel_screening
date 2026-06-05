#!/bin/bash
set -e


python helper_script/split_euk_pro.py

# Initialize conda for bash
source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

conda activate signalp6

BASE="$HOME/storage/thesis/test"
FASTA_DIR="$BASE/data/fasta_files"
OUT_BASE="$BASE/results/signalp"
MODELS="$HOME/tools/signalp/signalp6_slow_sequential/signalp-6-package/models"


for group in pro euk; do
    if [[ "$group" == "pro" ]]; then
        fasta="$FASTA_DIR/prokaryotic.faa"
        organism="other"
        outdir="$OUT_BASE/pro"
    else
        fasta="$FASTA_DIR/eukaryotic.faa"
        organism="euk"
        outdir="$OUT_BASE/euk"
    fi

    mkdir -p "$outdir"

    echo "Running SignalP6 for $group ..."
    signalp6 \
        --fastafile "$fasta" \
        --organism "$organism" \
        --output_dir "$outdir" \
        --format txt \
        --mode slow-sequential \
        --model_dir "$MODELS"
done

conda deactivate
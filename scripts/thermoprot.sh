#!/bin/bash
set -e

# Initialize conda for bash
source /home/ubuntu/miniconda3/etc/profile.d/conda.sh

conda activate thermoprot

python ~/tools/ThermoProt/thermoprot/thermoprot.py \
  --infile "$HOME/storage/thesis/test/data/fasta_files/filtered_hmmsearch.faa" \
  --outfile "$HOME/storage/thesis/test/output_files/thermoprot_prediction.csv" \
  --modelname MTH

conda deactivate

echo "Done!"
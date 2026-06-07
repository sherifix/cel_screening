#!/bin/bash


# Activate conda environment for vina tools
source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate vina


BASE_DIR="$(pwd)"
INPUT_CSV="${BASE_DIR}/output_files/pockets_summary.tsv"

# Skip header and loop through each line
tail -n +2 "$INPUT_CSV" | while IFS=$'\t' read -r protein pocket_score probability pocket_start_end residues center_x center_y center_z size_x size_y size_z trimmed_n num_atoms pdb_path; do
    # Create main directory for this protein (full path)
    PROT_DIR="${BASE_DIR}/results/autodock_vina/receptors/${protein}"
    mkdir -p "$PROT_DIR"

    # Compose box size and center strings
    BOX_CENTER="$center_x $center_y $center_z"
    BOX_SIZE="$size_x $size_y $size_z"

    echo "Processing ${protein} - Pocket ${pocket_start_end}..."

    # Run receptor preparation
    mk_prepare_receptor.py -i "$pdb_path" -o "${PROT_DIR}/${protein}_receptor" -p -v \
        --box_size $BOX_SIZE --box_center $BOX_CENTER
done

conda deactivate

echo "==============="
echo " Receptor preparation complete!"
echo "Outputs saved to: results/autodock_vina/receptors/"

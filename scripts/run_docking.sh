#!/bin/bash

# Activate conda environment for vina
source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
conda activate vina

BASE_DIR="$(pwd)"

# Ligand file (already prepared by prepare_ligand rule)
LIGAND="${BASE_DIR}/results/autodock_vina/cellotetraose.pdbqt"

# Exhaustiveness (higher = more accurate but slower)
EXHAUST=32

# Output folder
OUTPUT_DIR="${BASE_DIR}/results/autodock_vina/docking_results"
mkdir -p "$OUTPUT_DIR"

# Check if ligand exists
if [ ! -f "$LIGAND" ]; then
    echo "ERROR: Ligand file not found: $LIGAND"
    echo "Run prepare_ligand first"
    exit 1
fi

# Loop over all receptor directories
for RECEPTOR_DIR in "${BASE_DIR}/results/autodock_vina/receptors"/*; do
    [ -d "$RECEPTOR_DIR" ] || continue
    
    PROTEIN_NAME=$(basename "$RECEPTOR_DIR")
    echo "========================================="
    echo "Processing $PROTEIN_NAME..."
    echo "========================================="
    
    # Find receptor PDBQT file
    RECEPTOR_FILE="${RECEPTOR_DIR}/${PROTEIN_NAME}_receptor.pdbqt"
    
    if [ ! -f "$RECEPTOR_FILE" ]; then
        echo "  Warning: No receptor file found: $RECEPTOR_FILE"
        continue
    fi
    
    # Find corresponding box.txt file
    BOX_FILE="${RECEPTOR_DIR}/${PROTEIN_NAME}_receptor.box.txt"
    
    if [ ! -f "$BOX_FILE" ]; then
        echo "  Warning: No box file found: $BOX_FILE"
        continue
    fi
    
    # Output files
    OUT_FILE="${OUTPUT_DIR}/${PROTEIN_NAME}_out.pdbqt"
    LOG_FILE="${OUTPUT_DIR}/${PROTEIN_NAME}.log"
    
    # Run Vina
    echo "  Docking $LIGAND -> $PROTEIN_NAME..."
    echo "  Box: $(cat $BOX_FILE | head -1)"
    
    vina --receptor "$RECEPTOR_FILE" \
         --ligand "$LIGAND" \
         --config "$BOX_FILE" \
         --exhaustiveness $EXHAUST \
         --out "$OUT_FILE" > "$LOG_FILE" 2>&1
    
    if [ $? -eq 0 ]; then
        echo "  Docking Done → $OUT_FILE"
        # Extract binding affinity from log
        AFFINITY=$(grep "Affinity" "$LOG_FILE" | head -1 | awk '{print $4}')
        echo "  Binding affinity: $AFFINITY kcal/mol"
    else
        echo "   Docking Failed for $PROTEIN_NAME"
    fi
done

conda deactivate

echo ""
echo "All docking runs completed!"
echo "Results saved to: $OUTPUT_DIR"

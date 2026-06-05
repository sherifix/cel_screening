#!/bin/bash
set -e

CAND_DIR="data/query_structures"
REF_DIR="data/reference_structures"
OUTPUT_DIR="results/usalign_results"
HITS_CSV="output_files/hits_thermo_sp.csv"

mkdir -p "$OUTPUT_DIR"

# Read CSV and create accession -> family mapping
# Skip header, format: accession,family,...
tail -n +2 "$HITS_CSV" | while IFS=',' read -r accession family rest; do
    # Remove quotes if present
    accession=$(echo "$accession" | tr -d '"')
    family=$(echo "$family" | tr -d '"')
    
    cand_file="$CAND_DIR/${accession}.pdb"
    
    # Check if candidate structure exists
    if [ ! -f "$cand_file" ]; then
        echo "WARNING: No structure found for $accession"
        continue
    fi
    
    # Check if reference folder exists for this family
    ref_folder="$REF_DIR/$family"
    if [ ! -d "$ref_folder" ]; then
        echo "WARNING: No reference folder for family $family"
        continue
    fi
    
    echo "========================================="
    echo "Processing: $accession (Family: $family)"
    echo "========================================="
    
    # Align against all references in the family folder
    for ref in "$ref_folder"/*.pdb; do
        [ -f "$ref" ] || continue
        ref_base=$(basename "$ref" .pdb)
        
        out_prefix="$OUTPUT_DIR/${accession}_vs_${ref_base}"
        out_txt="${out_prefix}.txt"
        
        if [ -f "$out_txt" ]; then
            echo "  [SKIP] $accession vs $ref_base"
            continue
        fi
        
        echo "  [ALIGN] $accession vs $ref_base"
        USalign "$cand_file" "$ref" -o "$out_prefix" > "$out_txt" 2>/dev/null
    done
done

echo "[DONE] USalign alignments finished. Results in: $OUTPUT_DIR"

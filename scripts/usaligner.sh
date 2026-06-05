#!/bin/bash

CAND_DIR="/home/ubuntu/storage/thesis/data/hits_3d_structures"
REF_DIR="/home/ubuntu/storage/thesis/data/structure_hits/ref_structure"

for fam in "$CAND_DIR"/*; do
    [ -d "$fam" ] || continue
    family=$(basename "$fam")
    ref_folder="$REF_DIR/$family"

    if [ ! -d "$ref_folder" ]; then
        echo "!!!! No reference folder for $family"
        continue
    fi

    echo ">>>> Processing family: $family"

    for cand in "$fam"/*.pdb; do
        [ -f "$cand" ] || continue
        cand_base=$(basename "$cand" .pdb)

        for ref in "$ref_folder"/*.pdb; do
            [ -f "$ref" ] || continue
            ref_base=$(basename "$ref" .pdb)

            out_prefix="${fam}/${cand_base}_vs_${ref_base}"
            out_txt="${out_prefix}.txt"

            # Skip if already done
            if [ -f "$out_txt" ]; then
                echo "[SKIP] $cand_base vs $ref_base"
                continue
            fi

            echo "[ALIGN] $cand_base vs $ref_base"

            USalign "$cand" "$ref" -o "$out_prefix" \
                > "$out_txt"
        done
    done
done

echo "[DONE] USalign alignments finished."

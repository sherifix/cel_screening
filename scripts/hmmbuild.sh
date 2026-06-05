#!/bin/bash

INPUT_DIR="data/trimmed"
OUTPUT_DIR="data/hmmer_profiles"

mkdir -p "$OUTPUT_DIR"

for faa in "$INPUT_DIR"/*_trimmed.faa; do

    # Get basename without path: e.g. GH5_trimmed.faa
    basename_file=$(basename "$faa")

    # Strip _trimmed.faa → GH5
    family="${basename_file/_trimmed.faa/}"

    # Output profile path
    hmm_out="$OUTPUT_DIR/${family}.hmm"

    echo "Building profile for: $family"
    hmmbuild --amino "$hmm_out" "$faa"

done

echo "All profiles built"



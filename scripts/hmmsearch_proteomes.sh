#!/usr/bin/env bash
set -euo pipefail

PROFILES="data/hmmer_profiles"
PROTEOMES="data/proteomes"
OUTDIR="results/hmmsearch_results"

mkdir -p "$OUTDIR"

# Find all .faa and .fasta files recursively
find "$PROTEOMES" -type f \( -name "*.faa" -o -name "*.fasta" \) | while read -r faa; do
    # Get relative path without extension
    rel_path="${faa#$PROTEOMES/}"
    base="${rel_path%.faa}"
    base="${base%.fasta}"
    # Replace slashes with underscores for output directory name
    base_dir=$(echo "$base" | tr '/' '_')
    
    echo "Processing $rel_path"

    # Make output subdirectory for this proteome
    org_out="$OUTDIR/$base_dir"
    mkdir -p "$org_out"

    # Run hmmsearch for all profiles
    for hmm in "$PROFILES"/*.hmm; do
        if [[ ! -f "$hmm" ]]; then
            echo "No .hmm files found in $PROFILES"
            exit 0
        fi

        hmmbase=$(basename "$hmm" .hmm)
        hmmsearch --cpu 8 --domtblout "$org_out/${base_dir}_${hmmbase}.domtbl" \
            "$hmm" "$faa" > /dev/null 2>&1
    done
done

echo "hmmsearch complete"

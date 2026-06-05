#!/bin/bash
set -e

ACC_DIR="data/raw"
OUTPUT_DIR="data/GH_families"
MAX_RETRIES=3
RETRY_DELAY=2

mkdir -p "$OUTPUT_DIR"

# Function to download with retry
download_protein() {
    local acc=$1
    local outfile=$2
    local retry=0
    
    while [ $retry -lt $MAX_RETRIES ]; do
        if efetch -db protein -id "$acc" -format fasta >> "$outfile" 2>/dev/null; then
            return 0
        fi
        retry=$((retry + 1))
        echo "Retry $retry for $acc..." >&2
        sleep $RETRY_DELAY
    done
    echo "Failed to download $acc after $MAX_RETRIES attempts" >&2
    return 1
}

for file in "$ACC_DIR"/GH*.txt; do
    if [[ ! -e "$file" ]]; then
        echo "No accession files found in $ACC_DIR"
        break
    fi

    family=$(basename "$file" .txt)

    echo "Processing $family ..."

    sed 's/\r//g; s/"//g' "$file" | grep -v '^$' > acc_clean.txt

    num_acc=$(wc -l < acc_clean.txt)
    echo "Found $num_acc accessions for $family"

    if [[ $num_acc -eq 0 ]]; then
        echo "Skipping $family: no valid accessions"
        rm -f acc_clean.txt
        continue
    fi

    OUTFILE="$OUTPUT_DIR/${family}.faa"
    > "$OUTFILE"

    while read -r acc; do
        [[ -z "$acc" ]] && continue
        download_protein "$acc" "$OUTFILE"
    done < acc_clean.txt

    rm -f acc_clean.txt
    echo "Finished $family, saved to $OUTFILE"
done

echo "All done!"

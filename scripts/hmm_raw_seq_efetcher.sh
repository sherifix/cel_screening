#!/bin/bash
set -e

# Directory containing accession files
ACC_DIR="data/raw" 
OUTPUT_DIR="data/GH_families"

# Create parent output directory
mkdir -p "$OUTPUT_DIR"

# Loop over accession files
for file in "$ACC_DIR"/GH*.txt; do
    # Skip if no files found
    [[ ! -e "$file" ]] && echo "No accession files found in $ACC_DIR" && break

    # Extract family name: accs_GH44.txt -> GH44
    family=$(basename "$file" .txt)
    family=${family#accs_}

    echo "Processing $family ..."

    # Clean file: remove \r, quotes, empty lines
    sed 's/\r//g; s/"//g' "$file" | grep -v '^$' > acc_clean.txt

    # Count valid accessions
    num_acc=$(wc -l < acc_clean.txt)
    if [[ $num_acc -lt 10 ]]; then
        echo "Skipping $family: only $num_acc accessions"
        rm -f acc_clean.txt
        continue
    fi


    # Download proteins using efetch one accession at a time
    echo "Downloading $num_acc proteins for $family..."
    OUTFILE="$OUTPUT_DIR/${family}.faa"
    > "$OUTFILE"   # empty file before appending

    while read -r acc; do
        # skip empty lines
        [[ -z "$acc" ]] && continue
        efetch -db protein -id "$acc" -format fasta >> "$OUTFILE"
    done < acc_clean.txt

    # Clean up
    rm -f acc_clean.txt

    echo "Finished $family, saved to $OUTFILE"
done

echo "All GH families processed!"

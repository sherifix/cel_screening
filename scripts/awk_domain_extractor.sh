#!/bin/bash
set -e

FAA_DIR="data/GH_families"
DOMTBL_DIR="data/extracted_domains"
OUTPUT_DIR="data/domains"

mkdir -p "$OUTPUT_DIR"

MIN_LEN=110
EVALUE_MAX=1e-5

for domtbl in "$DOMTBL_DIR"/*_dom.tbl; do
    basename_file=$(basename "$domtbl")
    family="${basename_file/_dom.tbl/}"
    
    echo "Processing $family ..."

    FAA_FILE="$FAA_DIR/${family}.faa"
    if [[ ! -f "$FAA_FILE" ]]; then
        echo "ERROR: $FAA_FILE not found. Skipping."
        continue
    fi

    # Create deduplicated version of FASTA
    DEDUP_FASTA="${FAA_FILE}.dedup"
    awk '/^>/ {if (seen[$1]++) next} 1' "$FAA_FILE" > "$DEDUP_FASTA"
    
    # Ensure SSI index exists for dedup file
    if [[ ! -f "${DEDUP_FASTA}.ssi" ]]; then
        esl-sfetch --index "$DEDUP_FASTA"
    fi

    # Extract domains using deduped FASTA (no -s flag)
    awk -v fam="${family}.hmm" -v minlen="$MIN_LEN" -v evalue_max="$EVALUE_MAX" '
    $1 !~ /^#/ && $1 == fam && $13 < evalue_max && $20 > 0 && $21 > 0 && ($21 - $20 + 1) >= minlen {
        print $4 "/" $20 "-" $21, $20, $21, $4, $13
    }' "$domtbl" | sort -u | esl-sfetch -Cf "$DEDUP_FASTA" - > "$OUTPUT_DIR/${family}_domains.faa" 2>/dev/null
    
    # Clean up dedup file
    rm -f "$DEDUP_FASTA" "${DEDUP_FASTA}.ssi"
    
    if [[ ! -s "$OUTPUT_DIR/${family}_domains.faa" ]]; then
        echo "WARNING: No domains extracted for $family"
    else
        echo "OK: Extracted $(grep -c '^>' $OUTPUT_DIR/${family}_domains.faa) domains for $family"
    fi
done

echo "Domain extraction complete!"

#!/bin/bash
set -e

cd ~/storage/thesis/GH_families/preparing_hmmer_profiles

MIN_LEN=110
CDHIT_ID=0.8
THREADS=8

for domfile in *_dom.tbl; do
    family="${domfile%_dom.tbl}"
    echo "Processing $family ..."

    if [[ ! -f "${family}.faa" ]]; then
        echo "ERROR: ${family}.faa not found. Skipping."
        continue
    fi

    # Ensure SSI index exists
    if [[ ! -f "${family}.faa.ssi" ]]; then
        esl-sfetch --index "${family}.faa"
    fi

    # Extract exact catalytic domains safely
    awk -v fam="${family}.hmm" -v minlen="$MIN_LEN" -v evalue_max=1e-5 '
    $1 !~ /^#/ && $1 == fam && $13 < evalue_max && $20 > 0 && $21 > 0 && ($21 - $20 + 1) >= minlen {
        print $4 "/" $20 "-" $21, $20, $21, $4, $13
    }' "${domfile}" | esl-sfetch -Cf "${family}.faa" - > "${family}_domains.faa"
    
    # Skip if no domains extracted
    if [[ ! -s "${family}_domains.faa" ]]; then
        echo "No valid domains found for $family. Skipping."
        continue
    fi

    # Cluster with CD-HIT
    cd-hit -i "${family}_domains.faa" \
           -o "clustered_${family}.faa" \
           -c $CDHIT_ID -n 5 -d 0 -T $THREADS

    # Multiple sequence alignemnt
    mafft --maxiterate 1000 --localpair "clustered_${family}.faa" > "clustered_${family}_aligned.faa"

    #trim alignment
    trimal -in "clustered_${family}_aligned.faa" -gappyout -out "${family}_trimmed.faa"
    


    echo "Completed $family."
done

echo "All families processed!"

#!/usr/bin/env bash
set -euo pipefail

# Absolute paths
#ACCESSIONS="/home/ubuntu/storage/thesis/data/all_characterized_raw_data/all_genbank_accessions.txt"
PROFILES="/home/ubuntu/storage/thesis/test/data/hmmer_profiles"
PROTEOMES="/home/ubuntu/storage/thesis/test/data/proteomes"
OUTDIR="/home/ubuntu/storage/thesis/test/results/hmmsearch_results_tester"



# Loop over each TAXID_ACCESSION folder
for dir in "$PROTEOMES"/*/; do
    faa_file=("$dir"/*.faa)
    [[ ! -f "${faa_file[0]}" ]] && continue
    faa="${faa_file[0]}"
    base=$(basename "$faa" .faa)
    echo "Processing $base"

    # Make output subdirectory for this organism
    org_out="$OUTDIR/$base"
    mkdir -p "$org_out"

    # Run hmmsearch for all profiles     
    for hmm in "$PROFILES"/*.hmm; do
        [[ ! -f "${hmm}.h3f" ]] && hmmpress "$hmm"
        hmmbase=$(basename "$hmm" .hmm)
        hmmsearch --cpu 8 --domtblout "$org_out/${base}_${hmmbase}.domtbl" \
            "$hmm" "${faa}" > /dev/null
    done
    
done

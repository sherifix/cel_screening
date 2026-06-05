import glob


#use families with sequence >= 10
def get_valid_families():
    """Return families with at least 10 sequences in their original .faa file"""
    import glob
    valid = []
    for faa_file in glob.glob("data/GH_families/GH*.faa"):
        family = faa_file.replace("data/GH_families/", "").replace(".faa", "")

        with open(faa_file, 'r') as f:
            count = sum(1 for line in f if line.startswith('>'))
        if count >= 10:
            valid.append(family)
            print(f"Keeping {family}: {count} sequences")
        else:
            print(f"Filtering out {family}: only {count} sequences (need >=10)")
    return valid

VALID_FAMILIES = get_valid_families()


rule all:
    input:
        expand("data/trimmed/{family}_trimmed.faa", family=VALID_FAMILIES),
        "data/raw/.hmm_profiles_built_complete",
        "data/raw/.hmmsearch_complete",
        "data/raw/.hmmsearch_parsing_complete",
        "data/raw/.filtering_hmmsearch_complete",
        "data/raw/.thermostability_predicted",
        "data/raw/.signalp_prediction",
        "data/raw/.secreted_thermostable_results",
        "output_files/genbank_uniprot_map.csv",
        "data/raw/.usalign_complete",
        "output_files/usalign_summary.csv",
        "output_files/top_hits.csv",
        "data/trimmed_structures/trimming_summary.csv",
        "data/pdb_list.ds",
        "data/trimmed_pdb_list.ds",
        "data/raw/.p2rank_complete",
        "output_files/docking_boxes.csv"


rule download_cazy_data:
    input:
        "data/GH_families.txt"
    output:
        touch("data/raw/.cazy_download_complete")
    log:
        "logs/download_cazy_data.log"
    shell:
        """
        python scripts/download_cazy_characterized.py --families data/GH_families.txt 2> {log}
        """


rule prepare_accessions:
    input:
        "data/raw/.cazy_download_complete"
    output:
        touch("data/raw/.accs_prepared_complete")
    log:
        "logs/prepare_accessions.log"
    shell:
        """
        python scripts/accs_preparation.py 2> {log}
        """

rule fetch_sequences:
    input:
        "data/raw/.accs_prepared_complete"
    output:
        touch("data/raw/.sequences_fetched_complete")
    log:
        "logs/fetch_sequences.log"
    shell:
        """
        bash scripts/hmm_raw_seq_efetcher.sh 2> {log}
        """

rule extract_domains:
    input:
        "data/raw/.sequences_fetched_complete",
        "dbcan/dbCAN.hmm"
    output:
        touch("data/raw/.domains_extracted_complete")
    log:
        "logs/domain_extract.log"
    shell:
        """
        bash scripts/domain_extract.sh dbcan/dbCAN.hmm 2> {log}
        """

rule extract_domains_only:
    input:
        "data/raw/.domains_extracted_complete"
    output:
        touch("data/raw/.domains_extracted_complete_v2")
    log:
        "logs/extract_domains_only.log"
    shell:
        """
        bash scripts/awk_domain_extractor.sh 2> {log}
        """

rule cluster_domains:
    input:
        "data/domains/{family}_domains.faa"
    output:
        "data/clustered/{family}_clustered.faa"
    log:
        "logs/cluster_{family}.log"
    params:
        cdhit_id = 0.8,
        threads = 8
    run:
        import os
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        
        shell("""
        if [ ! -s {input} ]; then
            echo "WARNING: {input} is empty, skipping clustering" >> {log}
            touch {output}
        else
            cd-hit -i {input} -o {output} -c {params.cdhit_id} -n 5 -d 0 -T {params.threads} 2>> {log}
        fi
        """)


rule msa_domains:
    input:
        "data/clustered/{family}_clustered.faa"
    output:
        "data/alignments/{family}_aligned.faa"
    log:
        "logs/msa_{family}.log"
    params:
        maxiterate = 1000
    run:
        import os
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        
        shell("""
        if [ ! -s {input} ]; then
            echo "WARNING: {input} is empty, skipping alignment" >> {log}
            touch {output}
        else
            mafft --maxiterate {params.maxiterate} --localpair {input} > {output} 2>> {log}
        fi
        """)


rule trim_alignment:
    input:
        "data/alignments/{family}_aligned.faa"
    output:
        "data/trimmed/{family}_trimmed.faa"
    log:
        "logs/trim_{family}.log"
    run:
        import os
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        
        shell("""
        if [ ! -s {input} ]; then
            echo "WARNING: {input} is empty, skipping trimming" >> {log}
            touch {output}
        else
            trimal -in {input} -gappyout -out {output} 2>> {log}
        fi
        """)


rule build_hmm_profiles:
    input:
        expand("data/trimmed/{family}_trimmed.faa", family=VALID_FAMILIES)
    output:
        touch("data/raw/.hmm_profiles_built_complete")
    log:
        "logs/hmmbuild.log"
    shell:
        """
        bash scripts/hmmbuild.sh >> {log} 2>&1
        """

rule hmmsearch_proteomes:
    input:
        "data/raw/.hmm_profiles_built_complete"  
    output:
        touch("data/raw/.hmmsearch_complete")
    log:
        "logs/hmmsearch.log"
    shell:
        """
        bash scripts/hmmsearch_proteomes.sh 2> {log}
        """

rule download_blacklist:
    output:
        "data/blacklist_accessions.txt"
    log:
        "logs/download_blacklist.log"
    shell:
        """
        python scripts/download_blacklist.py 2> {log}
        """

rule hmmsearch_parsing:
    input:
        "data/raw/.hmmsearch_complete",
        "data/blacklist_accessions.txt"
    output:
        touch("data/raw/.hmmsearch_parsing_complete")
    log:
        "logs/hmmsearch_parse.log"
    shell:
        """
        python scripts/hmmsearch_results_parser.py 2> {log}
        """

rule filter_overlapping_domains:
    input:
        "data/raw/.hmmsearch_parsing_complete"
    output:
        touch("data/raw/.filtering_hmmsearch_complete")
    log:
        "logs/filter_hmmsearch_results.log"
    shell:
        """
        python scripts/filtering_overlapping_domains.py 2> {log}
        """

rule thermoprot_prediction:
    input:
        "data/raw/.filtering_hmmsearch_complete"
    output:
        touch("data/raw/.thermostability_predicted")
    log:
        "logs/thermoprot_prediction.log"
    shell:
        """
        bash scripts/thermoprot.sh 2> {log}
        """

rule signalp_prediction:
    input:
        "data/raw/.thermostability_predicted"
    output:
        touch("data/raw/.signalp_prediction")
    log:
        "logs/signalp_prediction.log"
    shell:
        """
        bash scripts/signalp.sh 2> {log}
        """

rule filtering_thermo_signal_results:
    input:
        "data/raw/.signalp_prediction"
    output:
        touch("data/raw/.secreted_thermostable_results")
    log:
        "logs/secreted_thermostable.log"
    shell:
        """
        python scripts/filtering_thermo_signal_predictions.py 2> {log}
        """

rule prepare_structures:
    input:
        "output_files/hits_thermo_sp.csv"
    output:
        "output_files/genbank_uniprot_map.csv"
    log:
        "logs/prepare_structures.log"
    shell:
        """
        python scripts/prepare_structures.py 2> {log}
        """

rule run_usalign:
    input:
        "output_files/genbank_uniprot_map.csv",
        "data/reference_structures"
    output:
        touch("data/raw/.usalign_complete")
    log:
        "logs/usalign.log"
    shell:
        """
        bash scripts/usaligner.sh 2> {log}
        """

rule parse_usalign:
    input:
        "data/raw/.usalign_complete"
    output:
        "output_files/usalign_summary.csv"
    log:
        "logs/parse_usalign.log"
    shell:
        """
        python scripts/parse_usalign.py 2> {log}
        """

rule filter_top_hits:
    input:
        "output_files/usalign_summary.csv",
        "output_files/hits_thermo_sp.csv"
    output:
        "output_files/top_hits.csv"
    log:
        "logs/filter_top_hits.log"
    shell:
        """
        python scripts/filter_top_hits.py 2> {log}
        """

rule trim_structures:
    input:
        "data/query_structures"
    output:
        "data/trimmed_structures/trimming_summary.csv"
    log:
        "logs/trim_structures.log"
    shell:
        """
        python scripts/trim_structures.py 2> {log}
        """

rule prepare_pdb_lists:
    input:
        "output_files/top_hits.csv",
        "data/query_structures",
        "data/trimmed_structures"
    output:
        "data/pdb_list.ds",
        "data/trimmed_pdb_list.ds"
    log:
        "logs/prepare_pdb_lists.log"
    shell:
        """
        python scripts/prepare_pdb_lists.py 2> {log}
        """

rule run_p2rank:
    input:
        "data/trimmed_pdb_list.ds"
    output:
        touch("data/raw/.p2rank_complete")
    log:
        "logs/p2rank.log"
    shell:
        """
        python scripts/run_p2rank.py 2> {log}
        """

rule parse_p2rank_results:
    input:
        "data/trimmed_pdb_list.ds",
        "data/trimmed_structures/trimming_summary.csv",
        "data/raw/.p2rank_complete"
    output:
        "output_files/docking_boxes.csv"
    log:
        "logs/parse_p2rank.log"
    shell:
        """
        python scripts/parse_p2rank_results.py 2> {log}
        """

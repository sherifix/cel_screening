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

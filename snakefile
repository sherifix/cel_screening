import glob
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

def get_families():
    """Get list of GH families from .faa files in GH_families directory"""
    faa_files = glob.glob("data/GH_families/GH*.faa")
    families = [f.replace("data/GH_families/", "").replace(".faa", "") for f in faa_files]
    return families

def get_all_trimmed_files():
    """Return all trimmed files for all families"""
    families = get_families()
    return [f"data/trimmed/{fam}_trimmed.faa" for fam in families]

def get_all_trimmed_files():
    """Return all trimmed files for all families"""
    return glob.glob("data/trimmed/GH*_trimmed.faa")

rule all:
    input:
        get_all_trimmed_files(),
	"data/raw/.hmm_profiles_built_complete",
	"data/raw/.hmmsearch_complete",
	"data/raw/.hmmsearch_parsing_complete",
	"data/raw/.filtering_hmmsearch_complete",
	"data/raw/.thermostability_predicted",
	"data/raw/.signalp_prediction",
	"data/raw/.secreted_thermostable_results",
	"data/raw/.cazy_filtered_complete"
	

rule download_cazy_data:
    output:
        touch("data/raw/.cazy_download_complete")
    log:
        "logs/download_cazy_data.log"
    shell:
        """
        python helper_script/characterized_cazy_data_downloader.py 2> {log}
        """

rule prepare_accessions:
    input:
        "data/raw/.cazy_download_complete"
    output:
        touch("data/raw/.accs_prepared_complete")
    log:
        "logs/accs_preparation.log"
    shell:
        """
        python helper_script/accs_preparation.py 2> {log}
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
        bash helper_script/hmm_raw_seq_efetcher.sh 2> {log}
        """

rule extract_domains:
    input:
        "data/raw/.sequences_fetched_complete"
    output:
        touch("data/raw/.domains_extracted_complete")
    log:
        "logs/domain_extract.log"
    shell:
        """
        bash helper_script/domain_extract.sh 2> {log}
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
        bash helper_script/awk_domain_extractor.sh 2> {log}
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
        expand("data/trimmed/{family}_trimmed.faa", family=get_families())
    output:
        touch("data/raw/.hmm_profiles_built_complete")
    log:
        "logs/hmmbuild.log"
    shell:
        """
        bash helper_script/hmmbuild.sh >> {log} 2>&1
        """


rule hmmsearch_proteomes:
    input:
        "data/raw/.hmm_profiles_built_complete"  # Wait for HMM profiles to be built
    output:
        touch("data/raw/.hmmsearch_complete")
    log:
        "logs/hmmsearch.log"
    shell:
        """
        bash helper_script/hmmsearch_proteomes.sh 2> {log}
        """

rule hmmsearch_parsing:
    input:
        "data/raw/.hmmsearch_complete"
    output:
        touch("data/raw/.hmmsearch_parsing_complete")
    log:
        "logs/hmmsearch_parse.log"
    shell:
        """
        python helper_script/hmmsearch_results_parser.py 2> {log}
        """


rule filter_overlapping_domains:
    input:
        "data/raw/.hmmsearch_parsing_complete"
    output:
        touch("data/raw/.filtering_hmmsearch_complete")
    log:
        "logs/filter_hmmsearch_results.log"
    shell:
        "python helper_script/filtering_overlapping_domains.py 2> {log}"


rule thermoprot_prediction:
    input:
        "data/raw/.filtering_hmmsearch_complete"
    output:
        touch("data/raw/.thermostability_predicted")
    log:
        "logs/thermoprot_prediction.log"
    shell:
        """
        bash helper_script/thermoprot.sh 2> {log}
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
        bash helper_script/signalp.sh 2> {log}
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
        python helper_script/filtering_thermo_signal_predictions.py
        """

rule download_cazy_filtered:
    input:
        "data/raw/.secreted_thermostable_results"
    output:
        touch("data/raw/.cazy_filtered_complete")
    log:
        "logs/cazy_filtered.log"
    shell:
        """
        bash helper_script/excluding_char_adding_lengths.sh 2> {log}
        """

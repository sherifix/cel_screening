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
        "data/raw/.structures_prepared",
        "output_files/opt_ph_prediction.csv",
        "output_files/ph_pi_summary.csv"


rule download_cazy_data:
    input:
        "data/GH_families.txt"
    output:
        expand("data/raw/{family}.csv", family=VALID_FAMILIES)
    log:
        "logs/download_cazy_data.log"
    shell:
        """
        python scripts/download_cazy_characterized.py --families data/GH_families.txt 2> {log}
        """


rule prepare_accessions:
    input:
        expand("data/raw/{family}.csv", family=VALID_FAMILIES)
    output:
        expand("data/raw/{family}.txt", family=VALID_FAMILIES)
    log:
        "logs/prepare_accessions.log"
    shell:
        """
        python scripts/accs_preparation.py 2> {log}
        """

rule fetch_sequences:
    input:
        expand("data/raw/{family}.txt", family=VALID_FAMILIES)
    output:
        expand("data/GH_families/{family}.faa", family=VALID_FAMILIES)
    log:
        "logs/fetch_sequences.log"
    shell:
        """
        bash scripts/hmm_raw_seq_efetcher.sh 2> {log}
        """

rule extract_domains:
    input:
        expand("data/GH_families/{family}.faa", family=VALID_FAMILIES),
        "dbcan/dbCAN.hmm"
    output:
        expand("data/extracted_domains/{family}_dom.tbl", family=VALID_FAMILIES)
    log:
        "logs/domain_extract.log"
    shell:
        """
        bash scripts/domain_extract.sh dbcan/dbCAN.hmm 2> {log}
        """

rule extract_domains_only:
    input:
        expand("data/extracted_domains/{family}_dom.tbl", family=VALID_FAMILIES)
    output:
        expand("data/domains/{family}_domains.faa", family=VALID_FAMILIES)
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
        expand("data/hmmer_profiles/{family}.hmm", family=VALID_FAMILIES)
    log:
        "logs/hmmbuild.log"
    shell:
        """
        bash scripts/hmmbuild.sh >> {log} 2>&1
        """

rule hmmsearch_proteomes:
    input:
        expand("data/hmmer_profiles/{family}.hmm", family=VALID_FAMILIES)
    output:
        directory("results/hmmsearch_results")
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
        "results/hmmsearch_results",
        "data/blacklist_accessions.txt"
    output:
        "results/raw_results_hmmsearch.csv"
    log:
        "logs/hmmsearch_parse.log"
    shell:
        """
        python scripts/hmmsearch_results_parser.py 2> {log}
        """

rule filter_overlapping_domains:
    input:
        "results/raw_results_hmmsearch.csv"
    output:
        "output_files/filtered_hmmsearch_results.csv"
    log:
        "logs/filter_hmmsearch_results.log"
    shell:
        """
        python scripts/filtering_overlapping_domains.py 2> {log}
        """

rule thermoprot_prediction:
    input:
        "output_files/filtered_hmmsearch_results.csv"
    output:
        "output_files/thermoprot_prediction.csv"
    log:
        "logs/thermoprot_prediction.log"
    shell:
        """
        bash scripts/thermoprot.sh 2> {log}
        """

rule signalp_prediction:
    input:
        "output_files/thermoprot_prediction.csv"
    output:
        directory("results/signalp")
    log:
        "logs/signalp_prediction.log"
    shell:
        """
        bash scripts/signalp.sh 2> {log}
        """

rule filtering_thermo_signal_results:
    input:
        "output_files/thermoprot_prediction.csv",
        "results/signalp"
    output:
        "output_files/hits_thermo_sp.csv"
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
        "output_files/genbank_uniprot_map.csv",
        touch("data/raw/.structures_prepared")
    log:
        "logs/prepare_structures.log"
    shell:
        """
        python scripts/prepare_structures.py 2> {log}
        """

rule run_usalign:
    input:
        "output_files/genbank_uniprot_map.csv",
        "data/raw/.structures_prepared"
    output:
        directory("results/usalign_results")
    log:
        "logs/usalign.log"
    shell:
        """
        bash scripts/usaligner.sh 2> {log}
        """

rule parse_usalign:
    input:
        "results/usalign_results"
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
        "data/raw/.structures_prepared"
    output:
        "data/trimmed_structures/trimming_summary.csv",
        directory("data/trimmed_structures")
    log:
        "logs/trim_structures.log"
    shell:
        """
        python scripts/trim_structures.py 2> {log}
        """

rule prepare_pdb_lists:
    input:
        "output_files/top_hits.csv",
        "data/raw/.structures_prepared",
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
        directory("results/p2rank_predictions")
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
        "results/p2rank_predictions"
    output:
        "output_files/pockets_summary.tsv"
    log:
        "logs/parse_p2rank.log"
    shell:
        """
        python scripts/parse_p2rank_results.py 2> {log}
        """

rule prepare_ligand:
    input:
        "results/autodock_vina"  
    output:
        "results/autodock_vina/cellotetraose.pdbqt"
    log:
        "logs/prepare_ligand.log"
    shell:
        """
        python scripts/prepare_ligand.py 2> {log}
        """

rule prepare_receptors:
    input:
        "output_files/pockets_summary.tsv"
    output:
        directory("results/autodock_vina/receptors")
    log:
        "logs/prepare_receptors.log"
    shell:
        """
        bash scripts/prepare_receptors.sh 2> {log}
        """

rule run_docking:
    input:
        "results/autodock_vina/receptors",
        "results/autodock_vina/cellotetraose.pdbqt"
    output:
        directory("results/autodock_vina/docking_results")
    log:
        "logs/run_docking.log"
    shell:
        """
        bash scripts/run_docking.sh 2> {log}
        """

rule parse_vina_results:
    input:
        "results/autodock_vina/docking_results"
    output:
        "output_files/vina_summary.csv"
    log:
        "logs/parse_vina.log"
    shell:
        """
        python scripts/parse_vina_results.py 2> {log}
        """


rule prepare_docked_fasta:
    input:
        "output_files/vina_summary.csv"
    output:
        "data/fasta_files/docked_hits.faa"
    log:
        "logs/prepare_docked_fasta.log"
    shell:
        """
        python scripts/prepare_docked_fasta.py 2> {log}
        """

rule run_ephod:
    input:
        "data/fasta_files/docked_hits.faa"
    output:
        "output_files/opt_ph_prediction.csv"
    log:
        "logs/run_ephod.log"
    shell:
        """
        bash scripts/run_ephod.sh 2> {log}
        """

rule calculate_pi:
    input:
        "data/fasta_files/docked_hits.faa",
        "output_files/top_hits.csv",
        "output_files/opt_ph_prediction.csv"
    output:
        "output_files/ph_pi_summary.csv"
    log:
        "logs/calculate_pi.log"
    shell:
        """
        python scripts/calculate_pi.py 2> {log}
        """

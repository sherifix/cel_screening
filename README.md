
# Cel_screening

**Pipeline for screening thermostable cellulases from GH family sequences**.

This pipeline was developed as a part of master thesis.

## Overview

This Snakemake pipeline processes GH family protein sequences to identify 
potential thermostable cellulases through:

- HMM profiling and domain annotation
- Thermostability prediction (ThermoProt)
- Signal peptide prediction (SignalP6)
- Structure prediction and analysis (AlphaFold DB, US-align)
- Pocket prediction (P2Rank)
- Molecular docking (AutoDock VINA)
- Predict optimum pH (EpHod)


The Carbohydrate Active enZYmes (CAZy) database serves as a source of the sequences
used for construction of HMMER profiles.


## System requirements 

- Linux OS (tested on Ubuntu)
- Minimum 16 GB RAM
- 30 GB free disk space
- Conda (Miniconda or Anaconda)
- Git
- Internet connection
- Available GUI for visualizations


## Installation

### 1. Clone the repository

```bash
cd 
git clone https://github.com/sherifix/cel_screening
cd cel_screening
```

### 2. Create conda environments 

```bash
cd environments
conda env create -f thesis.yml
conda env create -f vina.yml
cd ..

# install AutoDock Vina

conda activate vina
pip install vina
conda deactivate 
```

### 3. Install external tools separately

These tools should be installed manually. The pipeline expects them at the locations below 
with the exact environment names.

| Tool | Installation |  Expected location  | Environment name |
|------|--------------|---------------------|------------------|
| SignalP6 | https://github.com/fteufel/signalp-6.0 | ~/tools/signalp/ | signalp6
| ThermoProt | https://github.com/jafetgado/ThermoProt | ~/tools/ThermoProt | thermoprot
| EpHod | https://github.com/beckham-lab/EpHod | ~/tools/EpHod | ephod
| P2Rank | https://github.com/rdk/p2rank | ~/tools/p2rank_2.5.1 | p2rank

Please note that: 
- For SignalP6, the slow-sequential models were used. If you use another model, 
change the --mode parameter in scripts/signalp.sh
- All scripts use conda activate <env_name> - environments must have exactly these names
- Adjust paths in scripts if you install tools to different locations


## Usage

### Input preparation
**Required inputs**

1. Proteomes to be screened by the pipeline 

format should be .faa or .fasta.

The proteomes directory should be located in data/proteomes

```bash
cd ~/cel_screening
mkdir -p data/proteomes/
cd data/proteomes/

# keep all proteomic data here

cd ~/cel_screening
```

2. GH families

The HMMER profiles will be constructed from specific GH families. Those families should be
specified by the user. A GH_families.txt file shall be created with the desired GH families 
(one family per line). Subfamily level can be used. The GH_families.txt should be located
 inside data/ directory. 

<mark>Example of the GH_families.txt</mark>

GH5_1

GH7

GH9

```bash 
cd ~/cel_screening/data
touch GH_families.txt

#add families or sub families using nano / vim or echo
```

3. dbCAN HMM database (for catalytic domain annotation)

it should be located in dbcan/

```bash
conda activate thesis
cd ~/cel_screening
mkdir dbcan
wget -O dbcan/dbCAN.hmm https://pro.unl.edu/dbCAN2/download/run_dbCAN_database_
total/dbCAN.hmm
hmmpress dbCAN.hmm
cd ..
```

4. Reference structures for 3D structure alignment

Reference structures are ideally resolved PDB structures.
They should be located inside subdirectories of their GH_family
in data/reference_structures in .pdb format

<mark>Example</mark>:

data/reference_structures/GH7/1CEL.pdb

```bash
mkdir data/reference_structures/
# create subdirectories based on your GH families
```

5. ligand file for docking

In this pipeline, a universal ligand is docked against all potential candidates. In case of more ligands 
or specific ligand for a candidate, the docking script scripts/run_docking.sh should be updated. 

In case of one ligand, protonated ligand file in .sdf format should be located in results/Autodock_vina

```bash
cd ~/cel_screening
mkdir -p results/Autodock_vina

# place protonated ligand.sdf here

cd ../../
```


## Running the pipeline

```bash
conda activate thesis

# first dry run 
snakemake --cores 8 --dry-run

#full run 
snakemake --cores 8
```

## Output

Results will be saved in:
- <mark>output_files/</mark> - CSV and TSV summaries
- <mark>results/</mark> - This directory contain initial output file from all programs
- <mark>results/Autodock_vina/</mark> - Here you can find receptor for docking and docking results those files can be used for visualization


## Customization

### Modifying fitlters 
 Edit <mark>scripts/filter_top_hits.py</mark> to change"
- TM-score threshold (default: 0.75)
- RMSD threshold (default: 2.06) 
- Candide residue length (default: < 800)

 Edit ==scripts/run_docking.sh== to change: 
- <mark>EXHAUST=32</mark> (exhaustiveness)



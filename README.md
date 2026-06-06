# Cel_Screening

Pipeline for screening thermostable cellulases from GH family sequences.

## Table of Contents

1. [Environment Setup](#environment-setup)
2. [Required External Data](#required-external-data)
3. [Input Data Preparation](#input-data-preparation)
4. [Ligand Preparation](#ligand-preparation)
5. [Running the Pipeline](#running-the-pipeline)
6. [Output Files](#output-files)

---

## Environment Setup

### Create all required environments

```bash
# Main pipeline
conda env create -f environments/thesis.yml

# Docking
conda env create -f environments/vina.yml

# SignalP6
conda env create -f environments/signalp6.yml

# ThermoProt
conda env create -f environments/thermoprot.yml

# EpHod (optional)
conda env create -f environments/ephod.yml

---


conda activate thesis

---

mkdir -p dbcan
wget -O dbcan/dbCAN.hmm https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V11.txt
cd dbcan
hmmpress dbCAN.hmm
cd ..

---

GH5_5
GH6
GH7
GH9
GH12

---


conda activate thesis
snakemake --cores 8


---



## Then verify it worked

```bash
head -20 README.md


---

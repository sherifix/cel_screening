# Cel_Screening

Pipeline for screening thermostable cellulases from GH family sequences.

## Required external data

### dbCAN HMM database

Download dbCAN HMM file from:
https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V11.txt

Place it in `dbcan/dbCAN.hmm`:

```bash
mkdir -p dbcan
wget -O dbcan/dbCAN.hmm https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V11.txt
hmmpress dbcan/dbCAN.hmm

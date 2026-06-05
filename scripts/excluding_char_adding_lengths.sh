#!/bin/bash


python helper_script/download_all_characterized_cazy_accessions.py

python helper_script/remove_characterized_accession.py 

echo "CAZy data download and filtering complete!"
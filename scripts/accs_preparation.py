import pandas as pd
import os
import sys
import glob


#extract genbank accession and output them in .txt
def prepare_accessions(gh_families_dir="./data/raw", output_dir="./data/raw"):

    #find all .csv files in the directory
    csv_files = glob.glob(os.path.join(gh_families_dir, "GH*.csv"))

    if not csv_files:
        print(f"No CSV files found in {gh_families_dir}")
        return

    for csv_file in csv_files:
        #extract GH family from filename
        basename = os.path.basename(csv_file)
        gh_fam = basename.replace(".csv", "")

        # Check if file has content
        if os.path.getsize(csv_file) == 0:
            print(f"Skipping {gh_fam}: empty file")
            continue

        try:
            df = pd.read_csv(csv_file)
        except pd.errors.EmptyDataError:
            print(f"Skipping {gh_fam}: CSV has no data")
            continue
        except Exception as e:
            print(f"Skipping {gh_fam}: Error reading CSV - {e}")
            continue


        if df["GenBank_accession"].empty or df["GenBank_accession"].isna().all():
            print(f"Skipping {gh_fam}: no accessions")
            continue

        # save accessions to .txt
        output_file = os.path.join(output_dir, f"{gh_fam}.txt")
        df["GenBank_accession"].dropna().to_csv(
            output_file,
            index=False,
            header=False
        )


if __name__ == "__main__":
    prepare_accessions()





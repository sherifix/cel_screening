import requests
from bs4 import BeautifulSoup
import pandas as pd
import os
import sys
import argparse


# scraping characterized cazy webpage
#return dataframe representing tables of those webpages

def scrape_cazy_characterized(url, target_ecs=None):
    if target_ecs is None:
        target_ecs = {
            "3.2.1.4",
            "3.2.1.21",
            "3.2.1.91",
            "3.2.1.176"
        }

    filename = url.split("/")[-1]
    gh_family = filename.replace("_characterized.html", "")

    html = requests.get(url).text
    soup = BeautifulSoup(html, "html.parser")

    tables = soup.find_all("table")
    target_table = next(
        (t for t in tables if t.find("a", href=lambda h: h and "uniprot.org" in h)),
        None
    )
    if target_table is None:
        print(f"ATTENTION! {gh_family}: no characterized table — skipping")
        return pd.DataFrame()


    # Extract rows
    data_rows = []
    for row in target_table.find_all("tr"):
        cells = row.find_all("td")
        if len(cells) < 7:
            continue

        # EC number
        ec_link = cells[1].find("a")
        if not ec_link:
            continue
        ec_number = ec_link.get_text(strip=True)
        if ec_number not in target_ecs:
            continue

        # Organism filter
        organism_text = cells[3].get_text(strip=True).lower()
        if any(x in organism_text for x in ["uncultured", "unidentified", "unknown", "synthetic"]):
            continue
        else:
            organism_name = organism_text


        # Bold GenBank accession preferred (current)
        bold_acc = cells[4].find("b")
        genbank_acc = bold_acc.get_text(strip=True) if bold_acc else None

        # UniProt accession
        uniprot_link = cells[5].find("a")
        uniprot_acc = uniprot_link.get_text(strip=True) if uniprot_link else None

        # PDB accession bold preferred (current)
        pdb_cell = cells[6]
        pdb_bold = pdb_cell.find("b")
        if pdb_bold:
            pdb_ids = [pdb_bold.get_text(strip=True)]
        else:
            pdb_links = pdb_cell.find_all("a")
            pdb_ids = [a.get_text(strip=True) for a in pdb_links]

        # Append row
        data_rows.append({
            "Organism": organism_name,
            "GH_family": gh_family,
            "EC_number": ec_number,
            "GenBank_accession": genbank_acc,
            "UniProt_accession": uniprot_acc,
            "PDB": ", ".join(pdb_ids) if pdb_ids else None
        })

    # Convert to DataFrame
    df = pd.DataFrame(data_rows)
    df.to_csv(f"{gh_family}.csv", index = False)
    return df






if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download CAZy characterized data')
    parser.add_argument('--families', '-f',
                        help='Text file with GH families (one per line)')
    parser.add_argument('--families-list', '-l',
                        nargs='+',
                        help='List of GH families directly (e.g., -l GH5_5 GH6 GH7)')

    args = parser.parse_args()

    # Get GH families from either file or command line
    if args.families:
        with open(args.families, 'r') as f:
            GH_FAMILIES = [line.strip() for line in f if line.strip() and not line.startswith('#')]
        print(f"Loaded {len(GH_FAMILIES)} families from {args.families}")
    elif args.families_list:
        GH_FAMILIES = args.families_list
        print(f"Using {len(GH_FAMILIES)} families from command line")
    else:
        print("ERROR: Please provide --families FILE or --families-list 'GH5_5 GH6 GH7'")
        sys.exit(1)

    # Change to data/raw directory
    os.makedirs("./data/raw", exist_ok=True)
    os.chdir("./data/raw")
    for gh_fam in GH_FAMILIES:
        url = f'https://www.cazy.org/{gh_fam}_characterized.html'
        print(f'accessing {url}')
        scrape_cazy_characterized(url)
        print(f'{gh_fam} is DOWNLOADED in data/raw/')




import requests
from bs4 import BeautifulSoup
import pandas as pd

def scrape_cazy_characterized(url, target_ecs=None):
    """
    Scrape a CAZy characterized GH page.

    Parameters:
        url (str): URL of GH characterized page
        target_ecs (set): EC numbers to filter (optional)

    Returns:
        pd.DataFrame: columns = ['GH_family', 'EC_number', 'GenBank_accession', 'UniProt_accession', 'PDB']
    """
    if target_ecs is None:
        target_ecs = {
            "3.2.1.4",
            "3.2.1.21",
            "3.2.1.91",
            "3.2.1.176"
        }

    # Extract GH family from URL
    filename = url.split("/")[-1]
    gh_family = filename.replace("_characterized.html", "")

    # Fetch page
    html = requests.get(url).text
    soup = BeautifulSoup(html, "html.parser")

    # Find the table with UniProt links (robust)
    tables = soup.find_all("table")
    target_table = next(
        (t for t in tables if t.find("a", href=lambda h: h and "uniprot.org" in h)),
        None
    )
    if target_table is None:
        print(f"[INFO] {gh_family}: no characterized table — skipping")
        return pd.DataFrame()
        #raise RuntimeError("Could not find target CAZy table on the page.") #was used with initial download

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


        # Bold GenBank accession (preferred)
        bold_acc = cells[4].find("b")
        genbank_acc = bold_acc.get_text(strip=True) if bold_acc else None

        # UniProt accession
        uniprot_link = cells[5].find("a")
        uniprot_acc = uniprot_link.get_text(strip=True) if uniprot_link else None

        # PDB accession(s) - bold preferred
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



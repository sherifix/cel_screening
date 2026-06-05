import requests
from bs4 import BeautifulSoup
import re

OUTPUT_FILE = "/home/ubuntu/storage/thesis/test/data/CAZy_all_GenBank_accessions.txt"
TARGET_ECS = {"3.2.1.4", "3.2.1.21", "3.2.1.91", "3.2.1.176"}
keywords=[
        "endoglucanase",
        "exoglucanase",
        "cellobiohydrolase",
        "glucosidase",
    'cellulase'
    ]


all_accessions = set()
session = requests.Session()

for i in range(0, 195):
    url = f"https://www.cazy.org/GH{i}_characterized.html"
    print(f"[INFO] Accessing {url}")

    try:
        response = session.get(url, timeout=30)
        response.raise_for_status()
    except Exception as e:
        print(f"[WARN] Could not fetch {url}: {e}")
        continue

    soup = BeautifulSoup(response.text, "html.parser")
    tables = soup.find_all("table")
    target_table = next(
        (t for t in tables if t.find("a", href=lambda h: h and "uniprot.org" in h)),
        None
    )
    if target_table is None:
        print(f"[INFO] GH{i}: no characterized table found, skipping")
        continue

    for row in target_table.find_all("tr"):
        cells = row.find_all("td")
        if len(cells) < 7:
            continue

        protein_name = cells[0].get_text(strip = True).lower()

        
        # EC number check
        ec_link = cells[1].find("a")
        if not ec_link:
            continue
        ec_number = ec_link.get_text(strip=True)
        #if ec_number not in TARGET_ECS:
        #    continue

        if ec_number not in TARGET_ECS and not any(word in protein_name for word in keywords):
            continue
            

        # --- GenBank accessions (anchors + plain text) ---
        genbank_cell = cells[4]
        for text in genbank_cell.stripped_strings:
            text = text.strip()
            # pattern: GenBank protein accession (no '_')
            if "." in text and not "_" in text and re.match(r"[A-Z]{3,}\d+\.\d+", text):
                all_accessions.add(text)

# Write all unique accessions to one txt
with open(OUTPUT_FILE, "w") as f:
    for acc in sorted(all_accessions):
        f.write(acc + "\n")

print(f"[SUCCESS] Saved {len(all_accessions)} unique GenBank accessions to {OUTPUT_FILE}")

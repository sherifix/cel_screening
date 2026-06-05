from usable_scripts import cazy_scraper 
import os


GH_FAMILIES = [
    'GH5_5',
    'GH5_52',
    'GH5_2',
    'GH5_4',
    'GH5_9',
    'GH5_25',
    'GH5_26',
    'GH5_37',
    'GH5_38',
    'GH5_1',
    'GH45',
    'GH6',
    'GH7',
    'GH9',
    'GH44',
    'GH12',
    'GH48'
]



os.chdir("./data/raw")

for gh_fam in GH_FAMILIES:
    url = f'https://www.cazy.org/{gh_fam}_characterized.html'
    print(f'accessing {url}')
    cazy_scraper.scrape_cazy_characterized(url)
    print('Successful')
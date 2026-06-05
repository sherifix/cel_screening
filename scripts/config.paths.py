import os


ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# define all paths relative to ROOT
PATHS = {
    'data_dir': os.path.join(ROOT_DIR, 'data'),
    'raw_data': os.path.join(ROOT_DIR, 'data', 'raw'),
    'gh_families': os.path.join(ROOT_DIR, 'data', 'GH_families'),
    'trimmed': os.path.join(ROOT_DIR, 'data', 'trimmed'),
    'results': os.path.join(ROOT_DIR, 'results'),
    'logs': os.path.join(ROOT_DIR, 'logs'),
    'temp': os.path.join(ROOT_DIR, 'temp'),
}

def get_path(key):
    return PATHS.get(key, key)

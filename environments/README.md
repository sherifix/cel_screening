# Environment Files

## Required environments

| Environment | File | Installation |
|-------------|------|--------------|
| Main pipeline | `thesis.yml` | `conda env create -f thesis.yml` |
| Docking | `vina.yml` | `conda env create -f vina.yml` |

## External tools (install separately)

| Tool | Installation |
|------|--------------|
| SignalP6 | https://github.com/fmfi-compbio/signalp6 |
| ThermoProt | `pip install thermoprot` |
| EpHod | Custom installation (see documentation) |
| P2Rank | https://github.com/rdk/p2rank |

## Create environments

```bash
conda env create -f thesis.yml
conda env create -f vina.yml
conda activate thesis

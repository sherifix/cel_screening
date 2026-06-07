# Environment Files

## Required environments

| Environment | File | Installation |
|-------------|------|--------------|
| Main pipeline | `thesis.yml` | `conda env create -f thesis.yml` |
| Docking | `vina.yml` | `conda env create -f vina.yml` |

## External tools (install separately)

| Tool | Installation |
|------|--------------|
| SignalP6 | https://github.com/fteufel/signalp-6.0 |
| ThermoProt | https://github.com/jafetgado/ThermoProt |
| EpHod | https://github.com/beckham-lab/EpHod |
| P2Rank | https://github.com/rdk/p2rank |

## Create environments

```bash
conda env create -f thesis.yml
conda env create -f vina.yml
conda activate thesis

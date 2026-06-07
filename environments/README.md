# Environment Files

## Required environments

| Environment | File | Installation |
|-------------|------|--------------|
| Main pipeline | `thesis.yml` | `conda env create -f thesis.yml` |
| Docking | `vina.yml` | `conda env create -f vina.yml` |

## External tools (install separately)

The tools must be installed into the mentioned location below. Note that for SignalP
the slow-sequential models was used. If you use another model, please change --mode parameter in
scripts/signalp.sh to what you actually downloaded (fast / slow) 

| Tool | Installation |  Expected location
|------|--------------|---------------------|
| SignalP6 | https://github.com/fteufel/signalp-6.0 | ~/tools/signalp/
| ThermoProt | https://github.com/jafetgado/ThermoProt | ~/tools/ThermoProt
| EpHod | https://github.com/beckham-lab/EpHod | ~/tools/EpHod
| P2Rank | https://github.com/rdk/p2rank | ~/tools/p2rank_2.5.1

## Create environments

```bash
conda env create -f thesis.yml
conda env create -f vina.yml
conda activate thesis

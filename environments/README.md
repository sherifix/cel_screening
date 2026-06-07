# Environment Files

## Required environments

Thesis and Vina environments can be created using their respective .yml files
then use vina environment to install AutoDock Vina 

```bash
conda env create -f thesis.yml
conda env create -f vina.yml
conda activate thesis

# install vina 
conda activate vina
pip install vina
conda deactivate
```


## External tools (install separately)

These tools must be installed manually. The pipeline expects them at the locations
below with the exact environment names.

| Tool | Installation |  Expected location  | Environment name |
|------|--------------|---------------------|------------------|
| SignalP6 | https://github.com/fteufel/signalp-6.0 | ~/tools/signalp/ | signalp6
| ThermoProt | https://github.com/jafetgado/ThermoProt | ~/tools/ThermoProt | thermoprot
| EpHod | https://github.com/beckham-lab/EpHod | ~/tools/EpHod | ephod
| P2Rank | https://github.com/rdk/p2rank | ~/tools/p2rank_2.5.1 | p2rank


## Important notes

- For SignalP6, the slow-sequential models were used. If you use another model,
 change the --mode parameter in scripts/signalp.sh
- All scripts use conda activate <env_name> - environments must have exactly these names
- Adjust paths in scripts if you install tools to different locations


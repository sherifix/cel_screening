# Environment Package Lists

These files document the exact package versions used in this pipeline.

| Environment | File | Purpose |
|-------------|------|---------|
| thesis | thesis.txt | Main environment (HMMER, CD-HIT, MAFFT, etc.) |
| signalp6 | signalp6.txt | SignalP 6.0 for secretion prediction |
| thermoprot | thermoprot.txt | Thermostability prediction |
| vina | vina.txt | Molecular docking |
| EpHod | ephod.txt | Additional tool |

## Recreating

To recreate an environment:

conda create -n thesis --file environments/thesis.txt

Note: Due to conda dependency resolution, manual installation may be needed.

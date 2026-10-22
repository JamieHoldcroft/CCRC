# Climate Change Research Centre

Comparing CABLE and JULES (and the other PLUMBER2 models) against flux tower observations at the Australian PLUMBER2 sites.

## Setup

```
conda env create -f environment.yml
conda activate ccrc
```

Run every command from the repo root.

## Folders

| Folder | What is in it |
|---|---|
| `config/` | `models.yaml` (how each model is corrected in the build), `detection.yaml` (thresholds for step 6) |
| `data/` | raw tower files (`flux/`, `met/`), raw model files (`models/<MODEL>/<SITE>.nc`), built files (`harmonised/`). Not in git. |
| `scripts/` | building the harmonised files, checks, and failure detection |
| `hydrology/` | ET signatures, CSV tables and figures |
| `outputs/` | logs, checks and detection results from `scripts/`. Not in git. |

## scripts/: run in order

```
python scripts/step1_rename_files.py        # rename raw model files to <SITE>.nc (only for new model files)
python scripts/step2_inventory.py           # list time axis, units, variables of every model file
python scripts/step3_build_harmonised.py --all   # data/harmonised/<SITE>.nc
python scripts/step4_validate.py            # check the built files, score the models
```

`loaders.py`, `harmonise.py` and `validate.py` are not run on their own; steps 3 and 4 import them.
`view_cable_jules.py` opens CABLE and JULES vs tower plots for a look (`--site`, `--zoom`, `--rain`, `--ML`).

## hydrology/: run after step 3

```
python hydrology/step1_collect_results.py   # data -> CSVs in hydrology/results/
python hydrology/step2_make_report.py       # CSVs -> PNGs in hydrology/figures/
```

`et_signatures.py` holds the calculations and is imported by the others.
`plot_gap_filling.py SITE YYYY-MM-DD` zooms in on measured vs gap-filled Qle.
`gap_filling_checks.py` prints how many days, months and years the gap-filling rules remove per site.
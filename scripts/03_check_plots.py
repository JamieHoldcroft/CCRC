"""Draw the standard check plots for harmonised site files.

Run from the repo root:
    python scripts/03_check_plots.py AU-Cum AU-How   # named sites
    python scripts/03_check_plots.py --all           # every file in data/harmonised

Writes outputs/checks/<SITE>/diurnal.png, monthly.png and zoom.png.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # find ccrc/ from scripts/
from ccrc.harmonise import load_config
from ccrc.plots import check_plots

DATA = Path("data/harmonised")
OUT = Path("outputs/checks")


def main(args):
    if not args:
        sys.exit(__doc__)
    if args == ["--all"]:
        args = sorted(p.stem for p in DATA.glob("*.nc"))
    colours = {name: cfg["colour"] for name, cfg in load_config().items()}
    for site in args:
        check_plots(DATA / f"{site}.nc", colours, OUT / site)
        print(f"{site} -> {OUT / site}")


if __name__ == "__main__":
    main(sys.argv[1:])

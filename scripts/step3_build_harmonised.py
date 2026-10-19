"""Build data/harmonised/<SITE>.nc for one or more sites.

Run from the repo root:
    python scripts/step3_build_harmonised.py AU-Cum AU-How     # named sites
    python scripts/step3_build_harmonised.py --all             # every site in data/flux

For each site also writes, to outputs/build/:
    <SITE>_corrections.csv   every change applied to each model (time shift,
                             units, renames, rebuilt variables, missing ones)
    alignment.csv            forcing-echo lag per model x site (should be 0)
"""

import sys
from pathlib import Path

import pandas as pd

from harmonise import build_site, load_config, save

OUT = Path("data/harmonised")
LOGS = Path("outputs/build")


def main(args):
    if not args:
        sys.exit(__doc__)
    if args == ["--all"]:
        args = sorted({p.name.split("_")[0] for p in Path("data/flux").glob("*_Flux.nc")})
    cfg = load_config()
    LOGS.mkdir(parents=True, exist_ok=True)

    align_path = LOGS / "alignment.csv"
    align = pd.read_csv(align_path) if align_path.exists() else pd.DataFrame()
    for site in args:
        print(f"{site}:")
        ds, log, site_align = build_site(site, cfg)
        save(ds, OUT / f"{site}.nc")
        log.to_csv(LOGS / f"{site}_corrections.csv", index=False)
        if not align.empty:
            align = align[align.site != site]
        align = pd.concat([align, site_align], ignore_index=True)
        align.to_csv(align_path, index=False)

        bad = site_align[site_align.lag_steps.fillna(0) != 0]
        print(f"  -> {OUT / f'{site}.nc'} ({len(ds.source) - 1} models)")
        for _, r in bad.iterrows():
            print(f"  !! {r.model}: {r.echo_var} ({r.echo_strength} check) best matches the obs at lag {r.lag_steps:g} steps")
        for _, r in log[log.correction.str.startswith("WARNING")].iterrows():
            print(f"  !! {r.source} {r.variable}: {r.correction}")


if __name__ == "__main__":
    main(sys.argv[1:])

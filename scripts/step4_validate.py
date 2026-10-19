"""Check every harmonised site file for processing errors, and score the models.

Run from the repo root:
    python scripts/step4_validate.py              # all files in data/harmonised
    python scripts/step4_validate.py AU-Cum       # named sites

Writes to outputs/validation/:
    checks.csv       one row per site x source x variable x check (pass/warn/fail)
    metrics.csv      r, bias, RMSE, NME per site x source x variable, and whether
                     the model beats a linear regression on sunlight (1lin)
    obs_quality.csv  per site: share of each flux measured, energy closure, ...
and prints a summary of every failed check.
"""

import sys
from pathlib import Path

import pandas as pd

from validate import validate

DATA = Path("data/harmonised")
OUT = Path("outputs/validation")


def main(args):
    sites = args or sorted(p.stem for p in DATA.glob("*.nc"))
    align_path = Path("outputs/build/alignment.csv")
    alignment = pd.read_csv(align_path) if align_path.exists() else None
    checks, metrics, obs_q = validate([DATA / f"{s}.nc" for s in sites], alignment)

    OUT.mkdir(parents=True, exist_ok=True)
    checks.to_csv(OUT / "checks.csv", index=False)
    metrics.to_csv(OUT / "metrics.csv", index=False)
    obs_q.to_csv(OUT / "obs_quality.csv", index=False)

    pd.set_option("display.width", 200, "display.max_rows", 500, "display.max_colwidth", 70)
    counts = checks.groupby("status").size()
    print(f"\n{len(checks)} checks: " + ", ".join(f"{k} {v}" for k, v in counts.items()))
    print("\nWarnings + fails per check:")
    bad = checks[checks.status.isin(["warn", "fail"])]
    print(bad.groupby(["check", "status"]).size().unstack(fill_value=0))
    print("\nFAILED checks:")
    print(checks[checks.status == "fail"][["site", "source", "variable", "check", "value", "note"]]
          .to_string(index=False))
    print("\nShare of site x flux where the model beats 1lin (higher is better):")
    print(metrics.pivot_table(index="source", columns="variable", values="beats_1lin",
                              aggfunc="mean").round(2))
    print(f"\nWrote {OUT}/checks.csv, metrics.csv, obs_quality.csv")


if __name__ == "__main__":
    main(sys.argv[1:])

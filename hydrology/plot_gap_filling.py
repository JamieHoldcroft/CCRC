"""Where is the tower data gap-filled? One figure per site.

Run from the repo root:
    python hydrology/plot_gap_filling.py                 # all sites
    python hydrology/plot_gap_filling.py AU-Cum AU-How   # named sites only
    python hydrology/plot_gap_filling.py --Qle           # add daily mean Qle underneath

Each figure has one horizontal strip per tower variable. Every time step is
checked, then summarised per day, because a figure is ~2,000 pixels wide and a
site has up to 280,000 half-hours. The colour of each day is the share of its
steps that were gap-filled: grey = all measured, red = all gap-filled, pink =
a mix (e.g. Qle gap-filled at night only). "Measured" uses the same rule as the
analysis (et_signatures.is_measured: qc flag below 20). The strip label gives
the share of all steps that were gap-filled.

With --Qle, a panel underneath shows the tower's Qle averaged per day
(all steps, gap-filled included), so you can see what was happening when
the gaps occurred.

Saves to hydrology/figures/gap_filled/<SITE>.png
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import et_signatures as ets

# ============================================================================
# Settings
# ============================================================================
HARMONISED = ets.HARMONISED
OUT = HERE / "figures" / "gap_filled"
# Top to bottom. Qle, Precip, Rnet, Tair and Psurf feed the ET, rain and PET
# analysis; Qh is used for the energy-balance gap (et_eb.csv).
VARIABLES = ["Qle", "Qh", "Precip", "Rnet", "Tair", "Psurf"]
COLOURS = LinearSegmentedColormap.from_list("gaps", ["lightgrey", "tab:red"])   # 0 = all measured, 1 = all gap-filled


def plot_site(path, with_qle):
    site = path.stem
    ds = xr.open_dataset(path)
    times = ds.time.to_index()
    start, end = times[0], times[-1] + (times[1] - times[0])

    n_panels = 2 if with_qle else 1
    heights = [len(VARIABLES), 3] if with_qle else [len(VARIABLES)]
    fig, axes = plt.subplots(n_panels, 1, figsize=(14, 2.5 + 1.0 * sum(heights) / 2),
                             sharex=True, gridspec_kw={"height_ratios": heights}, squeeze=False)
    axes = axes[:, 0]

    # --- Top panel: one strip per variable ----------------------------------
    ax = axes[0]
    labels = []
    for row, var in enumerate(VARIABLES):
        y = len(VARIABLES) - 1 - row                    # first variable at the top
        if f"{var}_qc" not in ds:
            labels.append(f"{var} (no qc flag)")
            continue
        gapfilled = ~ets.is_measured(ds, var)                           # True/False per step
        share_per_day = gapfilled.astype(float).resample("1D").mean()   # 0-1 per day
        days = share_per_day.index
        # One row of coloured cells, one cell per day, drawn between y-0.35 and y+0.35
        ax.imshow(share_per_day.values[np.newaxis, :], aspect="auto", cmap=COLOURS,
                  vmin=0, vmax=1, interpolation="nearest",
                  extent=(mdates.date2num(days[0]), mdates.date2num(days[-1] + pd.Timedelta(days=1)),
                          y - 0.35, y + 0.35))
        labels.append(f"{var}  ({100 * gapfilled.mean():.0f}% gap-filled)")
    ax.set_yticks(range(len(VARIABLES)), labels[::-1])
    ax.set_ylim(-0.6, len(VARIABLES) - 0.4)
    ax.xaxis_date()
    ax.set_title(f"{site}: share of each day gap-filled (grey = measured, red = gap-filled)", loc="left")

    # --- Optional bottom panel: daily mean tower Qle ------------------------
    if with_qle:
        ax = axes[1]
        qle_daily = ds.Qle.sel(source="obs").to_series().resample("1D").mean()
        ax.plot(qle_daily.index, qle_daily.values, color="k", lw=0.6)
        ax.set_ylabel("tower Qle\n(daily mean, W/m2)")
        ax.grid(alpha=0.3)

    axes[-1].set_xlim(mdates.date2num(start), mdates.date2num(end))
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{site}.png", dpi=150)
    plt.close(fig)
    print(f"wrote {OUT / site}.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sites", nargs="*", help="site codes, e.g. AU-Cum (default: all)")
    parser.add_argument("--Qle", action="store_true",
                        help="add the tower's daily mean Qle as a line underneath")
    args = parser.parse_args()

    sites = args.sites or [p.stem for p in sorted(HARMONISED.glob("*.nc"))]
    for site in sites:
        plot_site(HARMONISED / f"{site}.nc", args.Qle)


if __name__ == "__main__":
    main()

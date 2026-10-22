"""Zoom in on the tower's Qle for a few days: measured vs gap-filled.

Run from the repo root:
    python hydrology/plot_gap_filling.py AU-Tum 2016-01-15             # that one day
    python hydrology/plot_gap_filling.py AU-Tum 2016-01-15 --days 3    # 3 days from that date

Black line = measured Qle, red line and dots = gap-filled Qle.
"Measured" uses the same rule as the analysis (et_signatures.is_measured:
qc flag below 20).

Saves to hydrology/figures/gap_filled/<SITE>_<DATE>.png (figures_fluxnet/ for DATASET = "fluxnet")
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import et_signatures as ets

OUT = HERE / f"figures{ets.SUFFIX}" / "gap_filled"


def plot_qle(site, date, days):
    ds = xr.open_dataset(ets.HARMONISED / f"{site}.nc")

    # the chosen days only
    start = pd.Timestamp(date)
    end = start + pd.Timedelta(days=days)
    qle = ds.Qle.sel(source="obs").to_series()[start:end - pd.Timedelta(minutes=1)]
    if qle.empty:
        first, last = ds.time.to_index()[[0, -1]]
        print(f"{site} has no data on {date}. Its record runs {first:%Y-%m-%d} to {last:%Y-%m-%d}.")
        return

    measured = ets.is_measured(ds, "Qle")[qle.index]     # True/False for the same steps

    # The red line: gap-filled steps, plus the step either side of each, so the red
    # line joins onto the black line instead of leaving a hole
    gapfilled = ~measured
    near_gap = gapfilled | gapfilled.shift(1, fill_value=False) | gapfilled.shift(-1, fill_value=False)
    qle_red = qle.where(near_gap)

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(qle.index, qle.values, color="k", linewidth=1, label="measured")            # every step in black
    ax.plot(qle_red.index, qle_red.values, color="tab:red", linewidth=1.5, label="gap-filled")   # red drawn on top
    ax.plot(qle.index[gapfilled], qle[gapfilled].values, "o", color="tab:red", markersize=2.5)  # a dot on each gap-filled step
    ax.axhline(0, color="0.7", linewidth=0.6)
    ax.set_ylabel("tower Qle (W/m2)")
    ax.set_title(f"{site}: tower Qle, {date} ({days} day{'s' if days > 1 else ''}), "
                 f"{100 * gapfilled.mean():.0f}% gap-filled")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.autofmt_xdate()
    fig.tight_layout()

    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{site}_{date}.png", dpi=150)
    plt.close(fig)
    print(f"wrote {OUT / f'{site}_{date}.png'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("site", help="site code, e.g. AU-Tum")
    parser.add_argument("date", help="first day, e.g. 2016-01-15")
    parser.add_argument("--days", type=int, default=3, help="how many days to show (default 1)")
    args = parser.parse_args()
    plot_qle(args.site, args.date, args.days)


if __name__ == "__main__":
    main()

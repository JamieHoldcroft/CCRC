"""Draw every detected failure period, with 10 days either side for context.

Run from the repo root (after scripts/08_benchmark_loss.py):
    python scripts/09_plot_periods.py

Reads outputs/detection/benchmark_loss_periods.csv and writes one PNG per period
to outputs/detection/plots/, named <site>_<variable>_<model>_<start>.png.
The figures use the same drawing code as scripts/07_cable_jules.py (obs black,
CABLE blue, JULES red, benchmarks dashed grey), half-hourly, with the detected
period shaded yellow.
"""

import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # save files, don't open windows
import matplotlib.pyplot as plt
import pandas as pd

PERIODS = Path("outputs/detection/benchmark_loss_periods.csv")
OUT = Path("outputs/detection/plots")
CONTEXT_DAYS = 10

# Reuse read_site and plot_site from scripts/07_cable_jules.py
spec = importlib.util.spec_from_file_location("viewer", "scripts/07_cable_jules.py")
viewer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(viewer)

periods = pd.read_csv(PERIODS, parse_dates=["start", "end"])
OUT.mkdir(parents=True, exist_ok=True)

for site, site_periods in periods.groupby("site"):
    data = viewer.read_site(site, ["Qle", "Qh"], with_ml=True)   # read each site once
    for _, p in site_periods.iterrows():
        start = p.start - pd.Timedelta(days=CONTEXT_DAYS)
        end = p.end + pd.Timedelta(days=CONTEXT_DAYS + 1)
        window = data.sel(time=slice(start, end))
        title = (f"{p.start.date()} to {p.end.date()}: {p.model} {p.variable} worse than every "
                 f"benchmark by {p.mean_loss:.0f} W/m2 on average (shaded)")
        fig = viewer.plot_site(window, site, ["Qle", "Qh"], None, title)
        for ax in fig.axes[:2]:
            ax.axvspan(p.start, p.end + pd.Timedelta(days=1), color="yellow", alpha=0.25, zorder=0)
        name = f"{site}_{p.variable}_{p.model}_{p.start.date()}.png"
        fig.savefig(OUT / name, dpi=100)
        plt.close(fig)
        print(f"wrote {OUT / name}")
    del data

print(f"\n{len(periods)} figures in {OUT}")

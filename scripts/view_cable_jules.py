"""Look at CABLE and JULES against the tower, site by site. Opens plot windows;
nothing is saved to disk.

Run from the repo root:
    python scripts/view_cable_jules.py                    # every site, whole record (weekly means)
    python scripts/view_cable_jules.py --site AU-How      # one site, whole record
    python scripts/view_cable_jules.py --ML               # also draw the ML benchmarks
    python scripts/view_cable_jules.py --zoom AU-Tum 2013-12-01 2014-01-02
                                                        # one window, half-hourly
Options you can add to any of the above:
    --vars Qle Qh       which fluxes to show (default: Qle Qh Rnet NEE GPP)
    --freq 7D           average to this step first (1D = daily, 7D = weekly, 1MS = monthly)
    --rain              add tower rainfall as grey bars hanging from the top of each
                        panel (7-day means, in mm/day, on an upside-down right axis).
                        The flux axis is made twice as tall so the lines sit in the
                        lower half and the bars in the upper half, like PLUMBER plots.

Whole-record plots default to weekly means, because 10+ years of half-hourly data
is a solid block. Zoom plots default to half-hourly.

With several sites, close each window to see the next one. Each figure is
closed and its data released before the next site is read, so memory stays flat.

Lines:
    obs     black. Only time steps the tower actually measured; gaps = gap-filled.
    CABLE   blue
    JULES   red
    --ML    the empirical benchmarks as thin dashed grey lines (they only report
            Qle, Qh and NEE, so they are missing from the other panels)

When averaging (whole record, or --freq), every line is averaged over the same
steps the tower measured, so the averages are comparable.
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

HARMONISED = Path("data/harmonised")
SKIP_SITES = ["AU-Ync"]  # its tower files are a copy of AU-ASM

UNITS = {"Qle": "W/m2", "Qh": "W/m2", "Rnet": "W/m2", "Qg": "W/m2",
         "NEE": "umol/m2/s", "GPP": "umol/m2/s"}
DEFAULT_VARS = ["Qle", "Qh", "Rnet", "NEE", "GPP"]   # Qg left out: tower plates are 8 cm deep

MODELS = {"CABLE": "#1f77b4", "JULES": "#d62728"}
BENCHMARKS = ["1lin", "3km27", "6km729", "6km729lag", "RF", "LSTM"]
BENCHMARK_GREYS = {b: str(0.75 - 0.1 * i) for i, b in enumerate(BENCHMARKS)}  # LSTM darkest

RAIN_FREQ = "7D"   # --rain: one bar per 7 days


def read_site(site, variables, with_ml, with_rain=False):
    """Open one site and keep only what we draw: obs, CABLE, JULES (+ benchmarks),
    and the tower rainfall if with_rain is True."""
    path = HARMONISED / f"{site}.nc"
    if not path.exists():
        sys.exit(f"No file {path}. Sites available: "
                 + ", ".join(p.stem for p in sorted(HARMONISED.glob('*.nc'))))
    with xr.open_dataset(path) as ds:
        sources = ["obs"] + list(MODELS)
        if with_ml:
            sources += BENCHMARKS
        sources = [s for s in sources if s in ds.source.values]
        keep = variables + [f"measured_{v}" for v in variables]
        if with_rain:
            keep.append("Precip_mm")                 # obs rainfall, mm per time step
        data = ds[keep].sel(source=sources).load()   # load only these, then close the file
    return data


def clamp_window(data, start, end):
    """Keep the zoom window inside the site's record. Too-early or too-late dates
    are moved to the first or last day, with a message; an empty window stops."""
    first = pd.Timestamp(data.time.values[0])
    last = pd.Timestamp(data.time.values[-1])
    start, end = pd.Timestamp(start), pd.Timestamp(end)
    if start >= end:
        sys.exit(f"Start {start.date()} is not before end {end.date()}.")
    if end < first or start > last:
        sys.exit(f"{data.attrs.get('site_code', '')} has data from {first.date()} "
                 f"to {last.date()}; {start.date()}..{end.date()} is outside it.")
    if start < first:
        print(f"Start moved from {start.date()} to the first day of data, {first.date()}")
        start = first
    if end > last:
        print(f"End moved from {end.date()} to the last day of data, {last.date()}")
        end = last
    return start, end


def weekly_rain(data):
    """Tower rainfall as 7-day means in mm/day: add up each day, then average
    the days in each 7-day block."""
    rain = data["Precip_mm"]
    if "source" in rain.dims:            # in case it is stored per source
        rain = rain.sel(source="obs")
    daily = rain.resample(time="1D").sum()          # mm per day
    return daily.resample(time=RAIN_FREQ).mean()    # mean mm/day in each 7-day block


def add_rain(ax, rain):
    """Grey rain bars hanging from the top of a panel, on an upside-down right axis."""
    # 1. Make the flux axis twice as tall: keep the bottom, double the range.
    bottom, top = ax.get_ylim()
    ax.set_ylim(bottom, bottom + 2 * (top - bottom))

    # 2. A second y axis on the right that shares the same time axis.
    rain_ax = ax.twinx()

    # 3. One bar per 7-day block, starting at the block's first day.
    rain_ax.bar(rain.time.values, rain.values, width=pd.Timedelta(RAIN_FREQ),
                align="edge", color="0.8", zorder=0)

    # 4. Turn the right axis upside down (0 at the top). Going to twice the
    #    wettest week means the longest bar reaches halfway down the panel.
    wettest = float(rain.max()) if float(rain.max()) > 0 else 1.0
    rain_ax.set_ylim(2 * wettest, 0)
    rain_ax.set_ylabel("Rain (mm/day,\n7-day mean)", color="0.4")

    # 5. Draw the flux lines in front of the bars.
    ax.set_zorder(rain_ax.get_zorder() + 1)
    ax.patch.set_visible(False)      # see-through, so the bars behind still show


def plot_site(data, site, variables, freq, title_extra):
    """One panel per flux: obs, CABLE, JULES (and benchmarks if loaded),
    plus rain bars if the data holds Precip_mm (--rain)."""
    fig, axes = plt.subplots(len(variables), 1, figsize=(15, 2.6 * len(variables)),
                             sharex=True, squeeze=False)
    rain = weekly_rain(data) if "Precip_mm" in data else None   # only with --rain
    for ax, var in zip(axes[:, 0], variables):
        measured = data[f"measured_{var}"].astype(bool)
        flux = data[var]
        if freq:
            # Average every source over the steps the tower measured, so the
            # averages cover the same hours.
            flux = flux.where(measured)
        else:
            
            # Half-hourly: hide the tower's gap-filled steps, draw models in full.
            flux = xr.where(flux.source == "obs", flux.where(measured), flux)

        for source in flux.source.values:
            line = flux.sel(source=source)
            if not np.isfinite(line.values).any():
                continue          # this source doesn't report this flux
            if freq:
                line = line.resample(time=freq).mean()
            if source == "obs":
                style = dict(color="k", lw=0.9, label="obs (measured)", zorder=5)
            elif source in MODELS:
                style = dict(color=MODELS[source], lw=0.8, label=source, zorder=4)
            else:
                style = dict(color=BENCHMARK_GREYS[source], lw=0.6, ls="--", label=source, zorder=3)
            ax.plot(line.time.values, line.values, **style)

        ax.axhline(0, color="0.85", lw=0.5, zorder=0)
        ax.set_ylabel(f"{var} ({UNITS[var]})")
        ax.margins(x=0.005)
        if rain is not None:
            add_rain(ax, rain)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center right", fontsize=8, frameon=False)
    step = f"{freq} means" if freq else "half-hourly"
    fig.suptitle(f"{site}: CABLE and JULES vs tower, {title_extra} ({step})")
    fig.tight_layout(rect=(0, 0, 0.90, 0.97))
    return fig


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--site", help="one site, e.g. AU-How (default: every site)")
    p.add_argument("--zoom", nargs=3, metavar=("SITE", "START", "END"),
                   help="one site and window, e.g. AU-Tum 2013-12-01 2014-01-02")
    p.add_argument("--ML", action="store_true", help="also draw the empirical benchmarks")
    p.add_argument("--vars", nargs="+", default=DEFAULT_VARS, choices=list(UNITS))
    p.add_argument("--freq", default=None, help="e.g. 1D, 7D, 1MS")
    p.add_argument("--rain", action="store_true",
                   help="add 7-day mean rainfall as bars hanging from the top")
    args = p.parse_args()

    if args.zoom:
        site, start, end = args.zoom
        data = read_site(site, args.vars, args.ML, args.rain)
        start, end = clamp_window(data, start, end)
        data = data.sel(time=slice(start, end))
        fig = plot_site(data, site, args.vars, args.freq, f"{start.date()} to {end.date()}")
        plt.show()
        plt.close(fig)
        return

    if args.site:
        sites = [args.site]
    else:
        sites = [p.stem for p in sorted(HARMONISED.glob("*.nc")) if p.stem not in SKIP_SITES]
    freq = args.freq or "7D"
    for site in sites:
        print(f"{site} (close the window to continue)")
        data = read_site(site, args.vars, args.ML, args.rain)
        fig = plot_site(data, site, args.vars, freq, "whole record")
        plt.show()
        plt.close(fig)
        del data   # release this site's data before reading the next


if __name__ == "__main__":
    main()

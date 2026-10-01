"""Thin line graphs of the observations and every model, one panel per flux.

    from ccrc.lineplots import line_plot
    line_plot("data/harmonised/AU-Cum.nc", "2013-12-11", "2013-12-21",
              out="outputs/lines/AU-Cum_zoom.png")

Each flux is a panel: x = time, y = flux, one thin line per source.
    obs   black (measured steps only; gap-filled steps are left blank)
    CABLE and JULES in their fixed colours (config/models.yaml), drawn on top
    every other model in its own colour, thinner and slightly transparent

Lines are 0.5-0.9 pt wide so 12 of them stay readable. For long windows, pass
freq="1D" (daily) or "7D" (weekly) to average first; half-hourly data over years
is a solid block however thin the lines are. When averaging, every source is
averaged over the same measured-obs steps.

Nothing here modifies the harmonised files.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # write files, never open windows
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

FLUXES = [("Qle", "W/m2"), ("Qh", "W/m2"), ("Rnet", "W/m2"),
          ("Qg", "W/m2"), ("NEE", "umol/m2/s"), ("GPP", "umol/m2/s")]
FOCUS = {"CABLE": "#1f77b4", "JULES": "#d62728"}
OTHER_COLOURS = plt.get_cmap("tab20")  # 20 distinct colours, one per other model

OBS_LW, FOCUS_LW, OTHER_LW = 0.9, 0.8, 0.5


def _colours(sources):
    """Fixed colours for CABLE/JULES, a distinct tab20 colour for every other model
    (config colours are shared within a family, so they can't tell ORCHIDEE2 from
    ORCHIDEE3)."""
    others = [s for s in sources if s != "obs" and s not in FOCUS]
    colours = dict(FOCUS)
    colours.update({s: OTHER_COLOURS(i % 20) for i, s in enumerate(others)})
    return colours


def line_plot(path, start=None, end=None, out=None, variables=None, freq=None,
              size=(16, 3.0)):
    """Draw one line graph per flux for `path` between `start` and `end`.

    start, end   anything xarray's .sel(time=slice(...)) accepts, e.g. "2013-12-11"
    variables    subset of FLUXES by name, e.g. ["Qle", "GPP"] (default: all six)
    freq         optional pandas frequency to average to before drawing ("1D", "7D")
    size         (width, height per panel) in inches
    """
    ds = xr.open_dataset(path).load()
    site = ds.attrs["site_code"]
    win = ds.sel(time=slice(start, end))
    if win.time.size == 0:
        raise ValueError(f"{site}: no data between {start} and {end} "
                         f"(file covers {ds.attrs['start']} to {ds.attrs['end']})")

    panels = [(v, u) for v, u in FLUXES if variables is None or v in variables]
    sources = list(win.source.values)
    colours = _colours(sources)
    order = [s for s in sources if s not in FOCUS and s != "obs"] + \
            [s for s in FOCUS if s in sources] + ["obs"]  # obs and focus on top

    fig, axes = plt.subplots(len(panels), 1, figsize=(size[0], size[1] * len(panels)),
                             sharex=True, squeeze=False)
    for ax, (var, units) in zip(axes[:, 0], panels):
        data = win[var]
        measured = win[f"measured_{var}"].astype(bool)
        # Half-hourly: obs drawn only where measured, models in full.
        # Averaged (freq): EVERY source is masked to the measured steps first, so
        # each mean covers the same steps. (Otherwise the obs mean would be mostly
        # daytime, because night fluxes are gap-filled more often, and the obs
        # would look biased against the models.)
        if freq:
            data = data.where(measured)
        obs = data.sel(source="obs").where(measured)
        for s in order:
            y = obs if s == "obs" else data.sel(source=s)
            if not np.isfinite(y.values).any():
                continue  # model does not report this flux
            if freq:
                y = y.resample(time=freq).mean()
            lw = OBS_LW if s == "obs" else FOCUS_LW if s in FOCUS else OTHER_LW
            ax.plot(y.time.values, y.values, lw=lw, color="k" if s == "obs" else colours[s],
                    alpha=1.0 if s in FOCUS or s == "obs" else 0.8, label=s)
        ax.axhline(0, color="0.85", lw=0.5, zorder=0)
        ax.set_ylabel(f"{var} ({units})")
        ax.margins(x=0.005)
    # one legend, outside the data, built from the first panel
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center right", fontsize=8, frameon=False,
               handlelength=2.5)
    span = f"{np.datetime_as_string(win.time.values[0], unit='D')} to " \
           f"{np.datetime_as_string(win.time.values[-1], unit='D')}"
    fig.suptitle(f"{site}: observations and models ({span}"
                 + (f", {freq} means" if freq else ", half-hourly") + ")")
    fig.tight_layout(rect=(0, 0, 0.93, 0.98))  # room for the legend on the right

    out = Path(out) if out else Path("outputs/lines") / f"{site}_{span.replace(' ', '')}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out

"""Standard check plots for one harmonised site file (data/harmonised/<SITE>.nc).

Style shared by every figure (and the future event cards):
    obs      black, measured steps only (gap-filled obs are not shown)
    CABLE    and JULES in their fixed colours from config/models.yaml
    others   grey band = min..max of the other land models

In the averaged plots (diurnal, monthly) every source is masked to the steps
where the obs were measured, so a model is always averaged over exactly the same
steps as the obs it is compared with. The zoom plot draws the models in full.
"""

import matplotlib

matplotlib.use("Agg")  # write files, never open windows
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

FLUXES = [("Qle", "W/m2"), ("Qh", "W/m2"), ("Rnet", "W/m2"),
          ("Qg", "W/m2"), ("NEE", "umol/m2/s"), ("GPP", "umol/m2/s")]
FOCUS = ["CABLE", "JULES"]
# Southern Hemisphere seasons
SEASONS = {"DJF (summer)": [12, 1, 2], "MAM (autumn)": [3, 4, 5],
           "JJA (winter)": [6, 7, 8], "SON (spring)": [9, 10, 11]}
BAND = dict(color="0.6", alpha=0.35, lw=0, label="other models (range)")


def measured(ds, var):
    """Every source's `var`, kept only where the obs were measured (qc == 0)."""
    mask = ds[f"measured_{var}"].astype(bool) if f"measured_{var}" in ds else True
    return ds[var].where(mask)


def other_models(ds):
    return [s for s, r in zip(ds.source.values, ds.role.values)
            if s not in FOCUS and r == "lsm"]


def _draw(ax, x, data, ds, colours):
    """obs line, CABLE/JULES lines and the grey band of the other models."""
    others = other_models(ds)
    if others:
        block = data.sel(source=others)
        if np.isfinite(block.values).any():
            ax.fill_between(x, block.min("source"), block.max("source"), **BAND)
    for s in FOCUS:
        if s in data.source and np.isfinite(data.sel(source=s).values).any():
            ax.plot(x, data.sel(source=s), color=colours[s], lw=1.2, label=s)
    ax.plot(x, data.sel(source="obs"), color="k", lw=1.5, label="obs (measured)")


def _title(ds, what):
    a = ds.attrs
    veg = a.get("IGBP_veg_long", a.get("IGBP_veg_short", ""))
    return f"{a['site_code']} ({a.get('site_name', '')}, {veg}): {what}"


def diurnal(ds, colours, path):
    """Mean daily cycle of each flux, one column per season."""
    hour = ds.time.dt.hour + ds.time.dt.minute / 60
    month = ds.time.dt.month
    fig, axes = plt.subplots(len(FLUXES), 4, figsize=(15, 14), sharex=True, sharey="row")
    for i, (var, units) in enumerate(FLUXES):
        data = measured(ds, var)
        for j, (season, months) in enumerate(SEASONS.items()):
            ax = axes[i, j]
            sub = data.where(month.isin(months))
            cycle = sub.groupby(hour).mean()
            _draw(ax, cycle[cycle.dims[-1]].values, cycle, ds, colours)
            ax.axhline(0, color="0.8", lw=0.5, zorder=0)
            if i == 0:
                ax.set_title(season)
            if j == 0:
                ax.set_ylabel(f"{var} ({units})")
            if i == len(FLUXES) - 1:
                ax.set_xlabel("hour of day (local standard time)")
    axes[0, -1].legend(fontsize=8, loc="upper right")
    fig.suptitle(_title(ds, "mean daily cycle by season, measured obs steps only"))
    fig.tight_layout(rect=(0, 0, 1, 0.97))  # leave room for the title
    fig.savefig(path, dpi=110)
    plt.close(fig)


def monthly(ds, colours, path, min_frac=0.2):
    """Monthly means over the whole record. Months with fewer than `min_frac`
    of steps measured are left blank, so a few odd measurements can't fake a mean."""
    fig, axes = plt.subplots(len(FLUXES), 1, figsize=(15, 14), sharex=True)
    for ax, (var, units) in zip(axes, FLUXES):
        data = measured(ds, var)
        means = data.resample(time="1MS").mean()
        frac = data.sel(source="obs").notnull().resample(time="1MS").mean()
        means = means.where(frac >= min_frac)
        _draw(ax, means.time.values, means, ds, colours)
        ax.axhline(0, color="0.8", lw=0.5, zorder=0)
        ax.set_ylabel(f"{var} ({units})")
    axes[0].legend(fontsize=8, loc="upper right", ncol=4)
    fig.suptitle(_title(ds, f"monthly means (months with >= {min_frac:.0%} measured)"))
    fig.tight_layout(rect=(0, 0, 1, 0.97))  # leave room for the title
    fig.savefig(path, dpi=110)
    plt.close(fig)


def pick_zoom(ds, days=10):
    """Start of the summer (Dec-Feb) window with the most measured Qle steps."""
    m = ds.measured_Qle.astype(bool).to_series()
    summer = m[m.index.month.isin([12, 1, 2])]
    if summer.empty:
        summer = m
    counts = summer.resample("1D").sum().rolling(days).sum()
    end = counts.idxmax()
    return end - np.timedelta64(days - 1, "D")


def zoom(ds, colours, path, days=10):
    """Half-hourly data for one summer window: sunlight and rain, then each flux.
    Gap-filled obs are drawn as a thin grey line so gaps are visible."""
    start = pick_zoom(ds, days)
    win = ds.sel(time=slice(start, start + np.timedelta64(days, "D")))
    fig, axes = plt.subplots(len(FLUXES) + 1, 1, figsize=(15, 16), sharex=True)
    ax = axes[0]
    ax.plot(win.time, win.SWdown, color="orange", lw=1, label="SWdown (W/m2)")
    ax.set_ylabel("SWdown (W/m2)")
    rain = ax.twinx()
    rain.bar(win.time, win.Precip_mm, width=np.timedelta64(30, "m"), color="tab:blue",
             label="rain (mm/step)")
    rain.set_ylabel("rain (mm per step)")
    for ax, (var, units) in zip(axes[1:], FLUXES):
        # Models drawn in full here (masking them would chop the lines into pieces);
        # only the obs are split into measured (black) and gap-filled (dotted).
        _draw(ax, win.time.values, win[var].where(win.source != "obs", measured(win, var)),
              win, colours)
        ax.plot(win.time, win[var].sel(source="obs"), color="0.5", lw=0.6, ls=":",
                label="obs (gap-filled)", zorder=0)
        ax.axhline(0, color="0.8", lw=0.5, zorder=0)
        ax.set_ylabel(f"{var} ({units})")
    axes[1].legend(fontsize=8, loc="upper right", ncol=5)
    fig.suptitle(_title(ds, f"half-hourly, {days} summer days with the most measured Qle"))
    fig.tight_layout(rect=(0, 0, 1, 0.97))  # leave room for the title
    fig.savefig(path, dpi=110)
    plt.close(fig)


def check_plots(path, colours, out_dir):
    """Write diurnal.png, monthly.png and zoom.png for one site."""
    ds = xr.open_dataset(path).load()
    out_dir.mkdir(parents=True, exist_ok=True)
    diurnal(ds, colours, out_dir / "diurnal.png")
    monthly(ds, colours, out_dir / "monthly.png")
    zoom(ds, colours, out_dir / "zoom.png")

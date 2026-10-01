"""Load flux-tower observations and land-model output into one aligned Dataset.

Layout of the combined Dataset returned by `combine`:
    dims:   time, source
    fluxes: Qle, Qh, Rnet, Qg, NEE, GPP   -> (source, time), source = ["obs", "JULES", ...]
    met:    Tair, VPD, Precip, SWdown, ... -> (time,)  shared by every source
    extras: *_qc flags, measured_<var> masks, Qle_cor, Qh_cor, ... -> (time,) obs only

Conventions checked against AU-Cum_2013-2018_OzFlux (PLUMBER2):
    - time is local standard time (SWdown peaks at 11:00-12:00), 30-min steps
    - qc flags: 0 = measured, 20/30/50 etc. = gap-filled
    - NEE is positive for release to the atmosphere (negative at midday)
"""

import warnings
from pathlib import Path

import numpy as np
import xarray as xr

# Fluxes that exist in both obs and model output. These get the `source` dimension.
FLUX_VARS = ["Qle", "Qh", "Rnet", "Qg", "NEE", "GPP"]

# Alternative names seen in PLUMBER2 / ALMA files -> the name used here.
RENAME = {
    "Rainf": "Precip",
    "PSurf": "Psurf",
    "Lwup_qc": "LWup_qc",  # AU-Cum Flux file has this typo
}


def _clean(ds):
    """Drop the 1x1 x/y grid, decode byte-string variables, standardise names."""
    # reset_coords first: in model files latitude/longitude are coordinates on the
    # x/y grid, and squeeze(drop=True) would throw them away along with x/y.
    ds = ds.reset_coords().squeeze(drop=True)

    for name in list(ds.data_vars):
        if ds[name].dtype.kind == "S":  # e.g. IGBP_veg_short is stored as bytes
            ds.attrs[name] = ds[name].values.item().decode().strip()
            ds = ds.drop_vars(name)
    return ds.rename({k: v for k, v in RENAME.items() if k in ds.variables})


def vpd_from_humidity(ds):
    """VPD (kPa) from Tair, Qair and Psurf: the humidity the models were forced with.

    The PLUMBER2 VPD variable is labelled hPa but holds kPa at most sites, and at
    AU-ASM, AU-Cow, AU-TTE, AU-Ync and AU-Rob it does not match the other humidity
    variables at all. Recomputing it removes both problems.
    """
    t = ds.Tair - 273.15
    es = 0.6108 * np.exp(17.27 * t / (t + 237.3))                      # kPa, Tetens
    e = ds.Qair * (ds.Psurf / 1000) / (0.622 + 0.378 * ds.Qair)       # kPa
    return (es - e).clip(min=0)


def load_site(met_path, flux_path):
    """Load one PLUMBER2 site (Met + Flux files) into a single Dataset.

    Adds:
        Tair_C     air temperature in degC
        VPD_kPa    vapour pressure deficit in kPa, recomputed (see vpd_from_humidity)
        Precip_mm  rainfall per time step in mm
        measured_<var>  True where <var>_qc == 0 (measured, not gap-filled)
    """
    met = _clean(xr.open_dataset(met_path))
    flux = _clean(xr.open_dataset(flux_path))

    # Static site info appears in both files; keep the Met copy.
    flux = flux.drop_vars([v for v in flux.data_vars if v in met.data_vars])
    ds = xr.merge([met, flux], combine_attrs="override")

    if not (ds.time.values == met.time.values).all():
        raise ValueError("Met and Flux files have different time axes")

    dt = float(np.median(np.diff(ds.time.values)) / np.timedelta64(1, "s"))
    ds.attrs["timestep_s"] = dt

    ds["Tair_C"] = ds.Tair - 273.15
    ds["Tair_C"].attrs = {"units": "degC", "long_name": "Near surface air temperature"}
    ds["VPD_kPa"] = vpd_from_humidity(ds)
    ds["VPD_kPa"].attrs = {"units": "kPa", "long_name": "Vapour pressure deficit",
                           "method": "computed from Tair, Qair, Psurf (Tetens)"}
    ds["Precip_mm"] = ds.Precip * dt  # kg/m2/s * s = mm per step
    ds["Precip_mm"].attrs = {"units": "mm per time step", "long_name": "Precipitation"}

    for var in FLUX_VARS:
        if f"{var}_qc" in ds:
            ds[f"measured_{var}"] = ds[f"{var}_qc"] == 0

    return ds


def _standardise(ds, name, path):
    """Keep only FLUX_VARS (missing ones as NaN) and label the Dataset with its source.

    Site latitude/longitude move to attrs so `combine` can check the model was
    run for the same site as the obs.
    """
    missing = [v for v in FLUX_VARS if v not in ds]
    if missing:
        warnings.warn(f"{name}: no {missing} in {path}; filled with NaN")
    for var in missing:
        ds[var] = xr.full_like(ds.time, np.nan, dtype=float)
    attrs = {"source": name, "path": str(path)}
    for coord in ["latitude", "longitude"]:
        if coord in ds.variables:
            attrs[coord] = float(ds[coord])
    # Drop every coordinate except time, so nothing clashes with the obs when merging.
    ds = ds[FLUX_VARS].reset_coords(drop=True)
    ds.attrs = attrs
    return ds


def load_model(path, name):
    """Load one land-model output file that already follows the obs conventions.

    ASSUMPTION: PLUMBER2-style ALMA netCDF, the same layout as the site Flux file:
    `time` in seconds since a stated start, singleton x/y dims, time stamps at
    the start of each step, same time zone as the obs, variables Qle, Qh, Rnet,
    Qg (W/m2) and NEE, GPP (umol/m2/s). Models that break any of this need their
    own loader, like `load_jules`. Run `check_conventions` on the combined data
    to catch violations.

    Returns a Dataset with only FLUX_VARS; any the model does not report are NaN.
    """
    ds = _clean(xr.open_dataset(path))
    return _standardise(ds, name, path)


# kg C -> g C (1e3) -> mol C (/12.011) -> umol (1e6). One mol of C fixed is one mol of CO2.
KGC_TO_UMOLCO2 = 1e3 / 12.011 * 1e6


def load_jules(path, name="JULES"):
    """Load JULES output for one site.

    Differences from the obs files, all checked on AU-Cum 2013-2018:
      - time is stored as float32 seconds. At ~1.9e8 s float32 can only step by
        16 s, so 15% of stamps are off by 8 s and would never match the obs.
        Rounded here to the nearest minute.
      - each stamp marks the END of its 30-min step (first stamp 00:30), while
        the obs stamp the START (first stamp 00:00). Using the start bound from
        `time_bounds` lines them up: JULES Rainf then equals obs Precip exactly,
        and is off by one step at any other shift. Without this, `combine`
        would find only one missing step and the whole record would be
        silently 30 min out.
      - GPP is in kg C m-2 s-1; converted to umol CO2 m-2 s-1 to match the obs.
      - there is no NEE (JULES reports GPP, NPP and plant respiration but not
        soil respiration), so NEE is NaN.
    Qle, Qh, Rnet and Qg already match the obs units and signs (Qg positive into
    the ground).
    """
    ds = xr.open_dataset(path, decode_times=False)
    start = ds.time_bounds.values[:, 0].astype("float64")
    start = np.round(start / 60) * 60
    time_attrs = {k: v for k, v in ds.time.attrs.items() if k != "bounds"}
    ds = ds.drop_vars("time_bounds").assign_coords(time=("time", start, time_attrs))
    ds = _clean(xr.decode_cf(ds))

    ds["GPP"] = ds.GPP * KGC_TO_UMOLCO2
    ds["GPP"].attrs = {"units": "umol/m2/s", "long_name": "Gross primary productivity"}
    ds["NEE"] = xr.full_like(ds.time, np.nan, dtype=float)  # not an output of these runs
    return _standardise(ds, name, path)


def combine(obs, models):
    """Put obs and model fluxes side by side along a new `source` dimension.
    Models are aligned onto the obs time axis.
    """
    fluxes = [obs[FLUX_VARS].assign_attrs(source="obs")]
    for m in models:
        for coord in ["latitude", "longitude"]:
            if coord in m.attrs and coord in obs and abs(m.attrs[coord] - float(obs[coord])) > 0.01:
                warnings.warn(
                    f"{m.attrs['source']}: {coord} {m.attrs[coord]:.2f} but obs has "
                    f"{float(obs[coord]):.2f}. Is this the right site?"
                )
        aligned = m.reindex(time=obs.time)
        frac_missing = float(aligned.Qle.isnull().mean())
        if frac_missing > 0.01:
            warnings.warn(
                f"{m.attrs['source']}: {frac_missing:.0%} of obs timesteps have no "
                "model value. Check time zone and time-stamp convention."
            )
        fluxes.append(aligned)

    names = [f.attrs["source"] for f in fluxes]
    stacked = xr.concat(fluxes, dim=xr.DataArray(names, dims="source", name="source"))

    out = xr.merge([obs.drop_vars(FLUX_VARS), stacked], combine_attrs="override")
    out.attrs = obs.attrs
    return out

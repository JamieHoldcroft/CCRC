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

# kg C -> g C (1e3) -> mol C (/12.011) -> umol (1e6). One mol of C fixed is one mol of CO2.
KGC_TO_UMOLCO2 = 1e3 / 12.011 * 1e6
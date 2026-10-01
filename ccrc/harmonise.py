"""Build one harmonised NetCDF per site: obs + every model on the obs time axis.

    data/harmonised/<SITE>.nc
        dims:    source (obs, CABLE, JULES, ...), time (obs time axis)
        coords:  role(source), family(source)  from config/models.yaml
        fluxes:  Qle, Qh, Rnet, Qg (W/m2), NEE, GPP (umol CO2/m2/s)  -> (source, time)
        obs:     met forcing, qc flags, measured_* masks, extra obs  -> (time,)

Conventions in the output, whatever the raw file did:
    - time stamps mark the START of each step, in the obs' local standard time
    - W/m2 for energy, umol CO2 m-2 s-1 for carbon
    - NEE positive = CO2 released to the air; Qg positive = into the ground
    - a flux a model does not report (and cannot be rebuilt exactly) is NaN

Every change made to a model's data is written to the correction log (returned
by `build_site`, saved as a CSV and in the file's `corrections` attribute).
"""

import re
import subprocess
import warnings
from datetime import datetime, timezone
from pathlib import Path

import cftime
import numpy as np
import pandas as pd
import xarray as xr
import yaml

from ccrc.loaders import FLUX_VARS, KGC_TO_UMOLCO2, load_site

ENERGY_VARS = ["Qle", "Qh", "Rnet", "Qg"]
CARBON_VARS = ["NEE", "GPP"]

# Extra model variables read only to rebuild a missing flux or check alignment.
HELPER_VARS = ["SWnet", "LWnet", "LWup", "SWdown", "rsds", "Rainf", "Precip"]

# After conversion, any |value| above this is a glitch or undeclared fill value
# (STEMMUS-SCOPE AU-Ync has one NEE/GPP step 1e6x too large): set to NaN and logged.
GLITCH_LIMIT = {"Qle": 5000, "Qh": 5000, "Rnet": 5000, "Qg": 5000, "NEE": 500, "GPP": 500}

# Normalised unit string -> factor to multiply by.
ENERGY_UNITS = {"w/m2": 1.0}
CARBON_UNITS = {
    "umol/m2/s": 1.0,
    "kg/m2/s": KGC_TO_UMOLCO2,          # kg C
    "g/m2/s": KGC_TO_UMOLCO2 / 1e3,     # g C
    "gco2/m2/s": 1e6 / 44.01,           # g CO2 (only via a units override in the config)
}

# Site attributes kept from the obs files; everything else (emails, processing
# history, software versions) is dropped. The raw files keep the full record.
KEEP_ATTRS = ["site_code", "site_name", "country", "site_description",
              "vegetation_description", "soil_type", "time_zone",
              "license_type", "acknowledgement"]
STATIC_VARS = ["latitude", "longitude", "elevation", "reference_height", "canopy_height"]

CITATION = ("Abramowitz et al. (2024) On the predictability of turbulent fluxes from "
            "land: PLUMBER2 MIP experimental description and preliminary results. "
            "Biogeosciences 21, 5517-5538.")


def load_config(path="config/models.yaml"):
    with open(path) as fh:
        return yaml.safe_load(fh)["models"]


def normalise_units(units):
    """'W m-2', 'W/m^2', 'gC/m^2/s', 'kg m-2 s-1' -> 'w/m2', 'g/m2/s', 'kg/m2/s'."""
    u = units.strip().lower().replace("^", "")
    u = re.sub(r"\s*m-2", "/m2", u)
    u = re.sub(r"\s*s-1", "/s", u)
    u = u.replace(" ", "").replace("gc/", "g/")
    return u


def unit_factor(var, units):
    table = ENERGY_UNITS if var in ENERGY_VARS else CARBON_UNITS
    key = normalise_units(units)
    if key not in table:
        raise ValueError(f"{var}: unrecognised units '{units}'. Add them to ccrc/harmonise.py")
    return table[key]


def model_times(ds, stamp, step_s, extra_steps=0):
    """Decode a model's time axis and move every stamp to the start of its step."""
    tname = "time" if "time" in ds.variables else "time_counter"
    t = ds[tname]
    # A time-zone suffix such as "+09:00" is only a label: the models were run on
    # the obs' local times, so it is dropped (checked by the forcing echo).
    units = re.sub(r"\s*[+-]\d{2}:\d{2}$", "", t.attrs["units"])
    raw = t.values.astype("float64")
    valid = np.isfinite(raw)  # BEPS pads some files with empty (fill-value) stamps
    dates = cftime.num2date(raw[valid], units, t.attrs.get("calendar", "standard"),
                            only_use_cftime_datetimes=False)
    good = pd.to_datetime([str(d) for d in dates], format="ISO8601")
    times = pd.Series(pd.NaT, index=range(raw.size), dtype="datetime64[ns]")
    times[valid] = good
    times = pd.DatetimeIndex(times)
    offset = {"start": 0.0, "middle": -0.5, "end": -1.0}[stamp] + extra_steps
    times = times + pd.to_timedelta(offset * step_s, unit="s")
    return times.round("min")  # float32 axes (CLM5, JULES) wobble by seconds; NaT stays NaT


def open_model(path, wanted, pft=None, log=None, name=""):
    """Open a model file, keep only `wanted` variables, drop 1x1 grid dims.

    pft: {"dim": ..., "weight": ...} for models that report some fluxes per
    plant functional type (ORCHIDEE). Those are combined into one grid-cell value
    as sum(value * weight) over the PFT dimension.
    """
    ds = xr.open_dataset(path, decode_times=False)
    if "time_counter" in ds.dims:
        ds = ds.rename({"time_counter": "time"})
    keep = [v for v in wanted if v in ds.data_vars]
    out = ds[keep + ["time"]].reset_coords(drop=True)
    if pft:
        weight = ds[pft["weight"]] if pft.get("weight") else 1
        how = f"weighted by {pft['weight']}" if pft.get("weight") else "unweighted"
        for v in keep:
            if pft["dim"] in out[v].dims:
                attrs = out[v].attrs
                out[v] = (out[v] * weight).sum(pft["dim"], min_count=1)
                out[v].attrs = attrs
                if log is not None:
                    log.append((name, v, f"per-PFT values summed ({how})"))
    # NoahMP etc. keep unused level coordinates (e.g. lev=4): drop those, squeeze 1x1 grids
    big = [d for d in out.dims if d != "time" and out.sizes[d] > 1]
    out = out.drop_dims(big)
    return out.squeeze([d for d in out.dims if d != "time"], drop=True), ds


def harmonise_model(name, cfg, path, obs, log):
    """One model at one site -> Dataset of FLUX_VARS on the obs time axis."""
    rename = cfg.get("rename", {})
    sign = cfg.get("sign", {})
    inverse = {std: raw for raw, std in rename.items()}
    raw_name = {v: inverse.get(v, v) for v in FLUX_VARS + HELPER_VARS}

    data, raw_ds = open_model(path, list(raw_name.values()), cfg.get("pft"), log, name)
    site = obs.attrs["site_code"]
    step_s = obs.attrs["timestep_s"]
    extra = cfg.get("shift_steps", {}).get(site, 0) + cfg.get("shift_minutes", 0) * 60 / step_s
    times = model_times(raw_ds, cfg["stamp"], step_s, extra)
    offset_steps = {"start": 0.0, "middle": -0.5, "end": -1.0}[cfg["stamp"]] + extra
    if offset_steps:
        log.append((name, "time", f"stamps moved {offset_steps * step_s / 60:+.0f} min "
                    f"({cfg['stamp']} of step -> start"
                    + (f", plus {extra:g} steps from shift_steps/shift_minutes" if extra else "") + ")"))
    if times.isna().any():
        log.append((name, "time", f"dropped {int(times.isna().sum())} steps with no time stamp"))
    data = data.assign_coords(time=times).isel(time=np.flatnonzero(~times.isna()))
    if data.indexes["time"].has_duplicates:
        raise ValueError(f"{name}: duplicate time stamps after shifting")

    out = xr.Dataset(coords={"time": obs.time})
    for var in FLUX_VARS + HELPER_VARS:
        rv = raw_name[var]
        if rv not in data:
            continue
        da = data[rv].astype("float64").reindex(time=obs.time)
        if var in cfg.get("drop_vars", []):
            continue  # becomes NaN below; logged as not reported
        if var in FLUX_VARS:
            units = cfg.get("units", {}).get(var, data[rv].attrs.get("units", ""))
            factor = unit_factor(var, units)
            da = da * factor * sign.get(var, 1)
            n_bad = int((abs(da) > GLITCH_LIMIT[var]).sum())
            if n_bad:
                da = da.where(abs(da) <= GLITCH_LIMIT[var])
                log.append((name, var, f"{n_bad} glitch values (|x| > {GLITCH_LIMIT[var]}) set to NaN"))
            notes = []
            if rv != var:
                notes.append(f"renamed from {rv}")
            if factor != 1:
                notes.append(f"units {units} x {factor:.6g}"
                             + (" (units overridden in config)" if var in cfg.get("units", {}) else ""))
            if sign.get(var, 1) == -1:
                notes.append("sign flipped")
            if notes:
                log.append((name, var, "; ".join(notes)))
        out[var] = da

    missing_after = float(out[[v for v in FLUX_VARS if v in out][0]].isnull().mean()) if any(
        v in out for v in FLUX_VARS) else 1.0
    if missing_after > 0.01:
        log.append((name, "time", f"WARNING {missing_after:.1%} of obs steps have no model value"))

    # Rebuild net radiation where it can be done exactly from the model's own terms.
    if "Rnet" not in out and {"SWnet", "LWnet"} <= set(out):
        out["Rnet"] = out.SWnet + out.LWnet
        log.append((name, "Rnet", "rebuilt as SWnet + LWnet"))
    elif "Rnet" not in out and {"SWnet", "LWup"} <= set(out):
        # LWdown is the forcing, identical for model and obs.
        out["Rnet"] = out.SWnet + obs.LWdown - out.LWup
        log.append((name, "Rnet", "rebuilt as SWnet + LWdown(forcing) - LWup"))

    for var in FLUX_VARS:
        if var not in out:
            out[var] = xr.full_like(obs.Tair, np.nan, dtype="float64")
            why = "dropped in config/models.yaml" if var in cfg.get("drop_vars", []) else "not reported"
            log.append((name, var, f"{why} -> NaN"))

    out.attrs = {"source": name, "path": str(path)}
    return out


def forcing_echo_lag(model, obs, max_lag=4, n_steps=17520):
    """Lag (in steps) at which the model's copy of the forcing best matches the obs.

    0 means the time axes line up. Uses SWdown or rain if the model writes them
    out (they are the forcing, so they should match exactly); otherwise SWnet or
    Rnet, which follow SWdown with no lag.
    """
    for mvar, ovar, strength in [("SWdown", "SWdown", "exact"), ("rsds", "SWdown", "exact"),
                                 ("Rainf", "Precip", "exact"), ("Precip", "Precip", "exact"),
                                 ("SWnet", "SWdown", "strong"), ("Rnet", "SWdown", "strong"),
                                 ("Qle", "SWdown", "weak")]:
        if mvar in model and np.isfinite(model[mvar].values).any():
            m = model[mvar].values[:n_steps]
            o = obs[ovar].values[:n_steps]
            best, best_r = None, -np.inf
            for lag in range(-max_lag, max_lag + 1):
                a, b = (m[lag:], o[:len(o) - lag]) if lag >= 0 else (m[:lag], o[-lag:])
                ok = np.isfinite(a) & np.isfinite(b)
                if ok.sum() < 100:
                    continue
                r = np.corrcoef(a[ok], b[ok])[0, 1]
                if r > best_r:
                    best, best_r = lag, r
            return {"echo_var": mvar, "echo_strength": strength, "lag_steps": best,
                    "r": round(float(best_r), 4)}
    return {"echo_var": None, "echo_strength": None, "lag_steps": None, "r": None}


def git_commit():
    try:
        sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True,
                             text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], capture_output=True,
                               text=True).stdout.strip()
        return sha + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def build_site(site, models_cfg, data_dir=Path("data")):
    """Return (harmonised Dataset, correction log DataFrame, alignment DataFrame)."""
    met = next((data_dir / "met").glob(f"{site}_*_Met.nc"))
    flux = next((data_dir / "flux").glob(f"{site}_*_Flux.nc"))
    obs = load_site(met, flux)
    raw_attrs = dict(obs.attrs)

    # load_site recomputes VPD_kPa; keep the file's value for reference.
    ratio = float((obs.VPD / obs.VPD_kPa).where(obs.VPD_kPa > 0.3).median())
    obs = obs.rename({"VPD": "VPD_file"})
    obs["VPD_file"].attrs["note"] = ("as stored in the PLUMBER2 Met file; units label "
                                     "unreliable, use VPD_kPa")

    log = [("obs", "VPD_kPa", f"recomputed from Tair/Qair/Psurf; file VPD / recomputed "
                              f"= {ratio:.3f} (median)")]
    sources, align = [obs[FLUX_VARS]], []
    for name, cfg in models_cfg.items():
        if cfg["role"] not in ("focus", "lsm", "sensitivity", "benchmark", "baseline"):
            continue
        path = data_dir / "models" / name / f"{site}.nc"
        if site in cfg.get("exclude_sites", []):
            log.append((name, "-", "excluded at this site in config/models.yaml"))
            continue
        if not path.exists():
            log.append((name, "-", "no file for this site"))
            continue
        model = harmonise_model(name, cfg, path, obs, log)
        align.append({"site": site, "model": name, **forcing_echo_lag(model, obs)})
        sources.append(model[FLUX_VARS])
        print(f"  {name}")

    names = ["obs"] + [s.attrs["source"] for s in sources[1:]]
    stacked = xr.concat([s.reset_coords(drop=True) for s in sources],
                        dim=pd.Index(names, name="source"), combine_attrs="drop")
    roles = ["obs"] + [models_cfg[n]["role"] for n in names[1:]]
    families = ["obs"] + [models_cfg[n]["family"] for n in names[1:]]
    stacked = stacked.assign_coords(role=("source", roles), family=("source", families))
    for var in ENERGY_VARS:
        stacked[var].attrs = {**obs[var].attrs, "units": "W/m2"}
    for var in CARBON_VARS:
        stacked[var].attrs = {**obs[var].attrs, "units": "umol/m2/s"}
    for var in FLUX_VARS:
        stacked[var].attrs.pop("_FillValue", None)

    # Sign check: every model flux should rise and fall with the obs.
    for var in FLUX_VARS:
        with warnings.catch_warnings():  # all-NaN columns (e.g. BEPS Qh) give NaN, not an error
            warnings.simplefilter("ignore", RuntimeWarning)
            r = xr.corr(stacked[var], stacked[var].sel(source="obs"), dim="time")
        for name, value in zip(names[1:], r.values[1:]):
            if np.isfinite(value) and value < 0:
                log.append((name, var, f"WARNING correlation with obs is {value:.2f}: "
                                       "sign flipped? set `sign` in config/models.yaml"))

    rest = obs.drop_vars(FLUX_VARS + [v for v in STATIC_VARS if v in obs])
    out = xr.merge([rest, stacked], combine_attrs="drop_conflicts")

    log_df = pd.DataFrame(log, columns=["source", "variable", "correction"])
    attrs = {k: raw_attrs[k] for k in KEEP_ATTRS if k in raw_attrs}
    for k in STATIC_VARS:
        if k in obs:
            attrs[k] = float(obs[k])
    for k in ["IGBP_veg_short", "IGBP_veg_long"]:
        if k in raw_attrs:
            attrs[k] = raw_attrs[k]
    attrs.update({
        "timestep_s": raw_attrs["timestep_s"],
        "time_convention": "stamp = start of step, local standard time (as the obs files)",
        "start": str(pd.Timestamp(obs.time.values[0])),
        "end": str(pd.Timestamp(obs.time.values[-1])),
        "dataset": "PLUMBER2 v1.0 (modelevaluation.org)",
        "citation": CITATION,
        "obs_files": f"{met.name}, {flux.name}",
        "model_files": ", ".join(f"models/{n}/{site}.nc" for n in names[1:]),
        "built": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "built_by": "scripts/02_build_site.py",
        "git_commit": git_commit(),
        "corrections": "\n".join(f"{s} {v}: {c}" for s, v, c in log),
    })
    out.attrs = attrs
    return out, log_df, pd.DataFrame(align)


def save(ds, path):
    """Write as compressed float32 (booleans as int8)."""
    enc = {}
    for name, da in ds.data_vars.items():
        if da.dtype == bool:
            enc[name] = {"dtype": "int8", "zlib": True}
        elif da.dtype.kind == "f":
            enc[name] = {"dtype": "float32", "zlib": True, "complevel": 4, "_FillValue": np.nan}
    path.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(path, encoding=enc)

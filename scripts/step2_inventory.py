"""Inventory every model output file before building the harmonised per-site files.

Run from the repo root:
    python scripts/step2_inventory.py                # full run (reads every file, slow)
    python scripts/step2_inventory.py --report-only  # re-print the summary from files.csv

Reads only metadata, the time axis and the key flux variables of each file in
data/models/<MODEL>/<SITE>.nc. Nothing is modified. Writes three tables to
outputs/inventory/:

    files.csv      one row per model x site: time axis (units, calendar, step,
                   first/last stamp, offset from the obs), lat/lon vs the obs
    variables.csv  one row per model x site x variable: name, dims, units
    key_vars.csv   one row per model: for each flux we compare (Qle, Qh, ...),
                   the name the model uses, its units, and the number of sites
                   that have it. "-" means the model never reports it.

and prints a short list of the things that need a decision (non-half-hourly
output, time offsets, unit differences, missing sites).
"""

import sys
from pathlib import Path

import cftime
import netCDF4
import numpy as np
import pandas as pd

MODELS_DIR = Path("data/models")
FLUX_DIR = Path("data/flux")
OUT_DIR = Path("outputs/inventory")

# Variables we want from every model, and the other names models use for them.
# CMIP-style names (hfls, gpp, ...) appear in CLM5 and LPJ-GUESS.
# NOTE: nep = -NEE (opposite sign), so it is listed but needs flipping later.
KEY_VARS = {
    "Qle": ["Qle", "hfls"],
    "Qh": ["Qh", "hfss"],
    "Rnet": ["Rnet"],
    "Qg": ["Qg", "hfds"],
    "NEE": ["NEE", "nee", "nep"],
    "GPP": ["GPP", "gpp"],
    "Evap": ["Evap", "et"],
}
TIME_NAMES = ["time", "time_counter"]
LAT_NAMES = ["latitude", "lat"]
LON_NAMES = ["longitude", "lon"]


def find(ds, names):
    """First of `names` that exists in the file, else None."""
    return next((n for n in names if n in ds.variables), None)


def scalar(ds, names):
    """Lat or lon as a plain float (they are stored as 0-d, 1-d or 2-d arrays)."""
    name = find(ds, names)
    if name is None:
        return np.nan
    return float(np.ravel(ds[name][:])[0])


def read_time(ds):
    """Decode the time axis. Returns a dict of facts, never raises."""
    name = find(ds, TIME_NAMES)
    if name is None:
        return {"time_name": None, "time_error": "no time variable"}
    t = ds[name]
    raw = np.asarray(t[:], dtype="float64")
    units = getattr(t, "units", "")
    calendar = getattr(t, "calendar", "standard")
    info = {
        "time_name": name,
        "time_units": units,
        "calendar": calendar,
        "time_dtype": str(t.dtype),
        "n_time": raw.size,
        "has_bounds": hasattr(t, "bounds") or f"{name}_bounds" in ds.variables
                      or f"{name}_bnds" in ds.variables,
    }
    try:
        stamps = cftime.num2date(raw[[0, -1]], units, calendar,
                                 only_use_cftime_datetimes=False)
        info["first_stamp"] = str(stamps[0])
        info["last_stamp"] = str(stamps[-1])
    except Exception as e:  # odd units such as "hours since ... +09:00"
        info["time_error"] = f"{type(e).__name__}: {e}"
    if raw.size > 1:
        # step in seconds, from the units' base (seconds / hours / days)
        per = {"seconds": 1, "minutes": 60, "hours": 3600, "days": 86400}
        factor = per.get(units.split(" ")[0].lower(), np.nan)
        steps = np.diff(raw) * factor
        info["step_s_median"] = float(np.median(steps))
        info["step_s_min"] = float(steps.min())
        info["step_s_max"] = float(steps.max())
    return info


def obs_reference(site):
    """First time stamp, step and lat/lon of the observed Flux file for a site."""
    files = sorted(FLUX_DIR.glob(f"{site}_*_Flux.nc"))
    if not files:
        return None
    with netCDF4.Dataset(files[0]) as ds:
        t = ds["time"]
        raw = np.asarray(t[:2], dtype="float64")
        first = cftime.num2date(raw[0], t.units, getattr(t, "calendar", "standard"),
                                only_use_cftime_datetimes=False)
        return {
            "obs_first": pd.Timestamp(str(first)),
            "obs_step_s": float(raw[1] - raw[0]),  # obs units are "seconds since ..."
            "obs_n_time": t.size,
            "obs_lat": scalar(ds, LAT_NAMES),
            "obs_lon": scalar(ds, LON_NAMES),
        }


def nan_fraction(var):
    """Share of values that are missing (NaN or masked as fill values)."""
    data = np.ma.masked_invalid(np.ma.asarray(var[:], dtype="float64"))
    return float(np.ma.getmaskarray(data).mean())


def inspect(path, model, site, ref):
    """Collect one files.csv row and the variables.csv rows for one file."""
    with netCDF4.Dataset(path) as ds:
        row = {"model": model, "site": site, "file": str(path),
               "dims": " ".join(f"{k}={len(v)}" for k, v in ds.dimensions.items())}
        row.update(read_time(ds))
        row["lat"] = scalar(ds, LAT_NAMES)
        row["lon"] = scalar(ds, LON_NAMES)

        if ref is not None:
            row["obs_n_time"] = ref["obs_n_time"]
            row["obs_step_s"] = ref["obs_step_s"]
            row["lat_diff"] = round(row["lat"] - ref["obs_lat"], 3)
            row["lon_diff"] = round(row["lon"] - ref["obs_lon"], 3)
            if "first_stamp" in row:
                # positive = model's first stamp is later than the obs' first stamp
                delta = pd.Timestamp(row["first_stamp"]) - ref["obs_first"]
                row["first_offset_min"] = delta.total_seconds() / 60

        for key, names in KEY_VARS.items():
            name = find(ds, names)
            row[f"{key}_name"] = name
            row[f"{key}_units"] = getattr(ds[name], "units", "") if name else None
            row[f"{key}_nan_frac"] = round(nan_fraction(ds[name]), 4) if name else None

        var_rows = [{"model": model, "site": site, "variable": name,
                     "dims": " ".join(v.dimensions),
                     "units": getattr(v, "units", ""),
                     "long_name": getattr(v, "long_name", "")}
                    for name, v in ds.variables.items()]
    return row, var_rows


def summarise_key_vars(files):
    """One row per model: name/units used for each key variable, and site count."""
    rows = []
    for model, g in files.groupby("model"):
        row = {"model": model, "n_sites": len(g)}
        for key in KEY_VARS:
            have = g[g[f"{key}_name"].notna()]
            if have.empty:
                row[key] = "-"
                continue
            names = "/".join(sorted(have[f"{key}_name"].unique()))
            units = "/".join(sorted(have[f"{key}_units"].astype(str).unique()))
            row[key] = f"{names} [{units}] ({len(have)}/{len(g)} sites)"
        rows.append(row)
    return pd.DataFrame(rows)


def report(files, sites):
    """Print the things that need a decision before building the harmonised files."""
    print("\n=== Needs attention ===")
    for model, g in files.groupby("model"):
        notes = []
        missing = sorted(set(sites) - set(g.site))
        if missing:
            notes.append(f"no file for {', '.join(missing)}")
        # CLM5 stores time as float32, so its steps wobble by ~1 s: allow 5 s
        off_step = g[(g.step_s_median - g.obs_step_s).abs() > 5]
        if not off_step.empty:
            steps = sorted(off_step.step_s_median.unique())
            notes.append(f"time step {steps} s differs from the obs at "
                         f"{len(off_step)} sites")
        if "time_error" in g and g.time_error.notna().any():
            notes.append(f"time not decodable: {g.time_error.dropna().iloc[0]}")
        if "first_offset_min" in g:
            # as a fraction of the step: 0 = start, 0.5 = middle, 1 = end of step
            frac = (g.first_offset_min * 60 / g.obs_step_s).round(3)
            shown = sorted(frac.dropna().unique())
            if shown and shown != [0.0]:
                notes.append(f"first stamp vs obs, in steps: {shown}")
        if (g.n_time != g.obs_n_time).any():
            n_bad = int((g.n_time != g.obs_n_time).sum())
            notes.append(f"length differs from obs at {n_bad} sites")
        far = g[(g.lat_diff.abs() > 0.01) | (g.lon_diff.abs() > 0.01)]
        if not far.empty:
            notes.append(f"lat/lon differs from obs at {', '.join(far.site)}")
        for key in KEY_VARS:
            units = g[f"{key}_units"].dropna().unique()
            if len(units) > 1:
                notes.append(f"{key} units differ between sites: {list(units)}")
            names = g[f"{key}_name"].dropna().unique()
            if len(names) and names[0] != key:
                notes.append(f"{key} is called {names[0]}")
            fully_nan = g[g[f"{key}_nan_frac"] == 1].site.tolist()
            if fully_nan:
                notes.append(f"{key} all missing at {', '.join(fully_nan)}")
        print(f"\n{model}:")
        for n in notes or ["nothing unusual"]:
            print(f"  - {n}")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sites = sorted({p.name.split("_")[0] for p in FLUX_DIR.glob("*_Flux.nc")})
    if "--report-only" in sys.argv:
        report(pd.read_csv(OUT_DIR / "files.csv"), sites)
        return
    refs = {s: obs_reference(s) for s in sites}

    file_rows, var_rows = [], []
    for model_dir in sorted(p for p in MODELS_DIR.iterdir() if p.is_dir()):
        for path in sorted(model_dir.glob("*.nc")):
            site = path.stem
            row, vrows = inspect(path, model_dir.name, site, refs.get(site))
            file_rows.append(row)
            var_rows.extend(vrows)
        print(f"read {model_dir.name}")

    files = pd.DataFrame(file_rows)
    variables = pd.DataFrame(var_rows)
    key_vars = summarise_key_vars(files)

    files.to_csv(OUT_DIR / "files.csv", index=False)
    variables.to_csv(OUT_DIR / "variables.csv", index=False)
    key_vars.to_csv(OUT_DIR / "key_vars.csv", index=False)
    print(f"\nWrote {len(files)} files, {len(variables)} variable rows to {OUT_DIR}/")

    report(files, sites)


if __name__ == "__main__":
    main()

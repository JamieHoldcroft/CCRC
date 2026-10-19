"""Data checks and skill metrics for harmonised site files.

Two different questions, kept apart on purpose:

1. Is the data CORRECT? (checks.csv)  Things that would mean a processing error
   rather than a bad model: wrong timing, units, signs, fill values, flat lines,
   energy or water not adding up, impossible values.
   Each check gives status pass / warn / fail with the value and threshold.

2. How GOOD is each model? (metrics.csv)  r, bias, RMSE and NME against the
   measured obs, plus the NME of a linear regression on sunlight ("1lin") as a
   yardstick. A model can pass every check in (1) and still score badly in (2):
   that is model error, which is what we are looking for.

Plus obs_quality.csv: how much of each obs flux is measured, and how well the
tower's energy balance closes.
"""

import numpy as np
import pandas as pd
import xarray as xr

FLUXES = ["Qle", "Qh", "Rnet", "Qg", "NEE", "GPP"]
# Plausible half-hourly bounds for a flux tower site (generous on purpose).
BOUNDS = {"Qle": (-150, 900), "Qh": (-400, 900), "Rnet": (-250, 1150),
          "Qg": (-400, 500), "NEE": (-70, 50), "GPP": (-5, 80)}
LATENT = 2.45e6  # J/kg, to turn Qle into evaporation


def status(value, warn, fail, higher_is_bad=True):
    if value is None or not np.isfinite(value):
        return "n/a"
    if higher_is_bad:
        return "fail" if value > fail else "warn" if value > warn else "pass"
    return "fail" if value < fail else "warn" if value < warn else "pass"


def longest_flat_fraction(x, min_run=48):
    """Share of non-NaN steps inside runs of >= min_run identical values (a stuck
    output). Zeros are ignored: GPP is legitimately 0 every night."""
    ok = np.isfinite(x)
    if ok.sum() == 0:
        return np.nan
    v = np.where(ok, x, np.nan)
    same = np.r_[False, v[1:] == v[:-1]] & ok & (v != 0)
    # run lengths of consecutive `same`
    edges = np.diff(np.r_[0, same.astype(int), 0])
    starts, ends = np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)
    flat = sum(e - s + 1 for s, e in zip(starts, ends) if e - s + 1 >= min_run)
    return flat / ok.sum()


def peak_hour(series, hours):
    """Hour of day of the mean diurnal maximum."""
    df = pd.DataFrame({"x": series, "h": hours}).dropna()
    if df.empty:
        return np.nan
    return float(df.groupby("h").x.mean().idxmax())


def nme(m, o):
    """Normalised mean error (PLUMBER): sum|m-o| / sum|o-mean(o)|. 0 = perfect,
    1 = no better than always predicting the obs mean."""
    ok = np.isfinite(m) & np.isfinite(o)
    if ok.sum() < 100:
        return np.nan
    return np.abs(m[ok] - o[ok]).sum() / np.abs(o[ok] - o[ok].mean()).sum()


def one_lin(sw, o):
    """In-sample linear regression of the flux on SWdown: PLUMBER's '1lin'
    benchmark, but fitted on the same data it is scored on (so it is a bit
    easier to beat a real out-of-sample 1lin)."""
    ok = np.isfinite(sw) & np.isfinite(o)
    if ok.sum() < 100:
        return np.full_like(o, np.nan)
    b, a = np.polyfit(sw[ok], o[ok], 1)
    return a + b * sw


def check_site(ds, alignment=None):
    """Return (checks, metrics, obs_quality) DataFrames for one site."""
    site = ds.attrs["site_code"]
    step_s = ds.attrs["timestep_s"]
    sw = ds.SWdown.values
    day, night = sw > 50, sw < 1
    hours = (ds.time.dt.hour + ds.time.dt.minute / 60).values
    months = pd.DatetimeIndex(ds.time.values).to_period("M")
    precip_mm = ds.Precip_mm.values
    rows, metrics = [], []

    def add(source, var, check, value, warn, fail, higher_is_bad=True, note=""):
        rows.append({"site": site, "source": source, "variable": var, "check": check,
                     "value": None if value is None else round(float(value), 4),
                     "warn_at": warn, "fail_at": fail,
                     "status": status(value, warn, fail, higher_is_bad), "note": note})

    obs = {v: ds[v].sel(source="obs").values for v in FLUXES}
    meas = {v: ds[f"measured_{v}"].values.astype(bool) for v in FLUXES}

    for source in ds.source.values:
        if source == "obs":
            continue
        m = {v: ds[v].sel(source=source).values.astype("float64") for v in FLUXES}
        have = [v for v in FLUXES if np.isfinite(m[v]).any()]

        # --- timing --------------------------------------------------------
        if alignment is not None:
            a = alignment[(alignment.site == site) & (alignment.model == source)]
            if len(a) and pd.notna(a.lag_steps.iloc[0]):
                add(source, a.echo_var.iloc[0], "forcing_echo_lag_steps",
                    abs(a.lag_steps.iloc[0]), 0, 0, note=f"r={a.r.iloc[0]}")
        for v in ["Qh", "Rnet"]:
            if v in have:
                diff = abs(peak_hour(m[v], hours) - peak_hour(np.where(meas[v], obs[v], np.nan), hours))
                add(source, v, "diurnal_peak_hour_diff", diff, 1.0, 2.5)

        for v in have:
            x = m[v]
            # --- completeness and broken values ----------------------------
            add(source, v, "missing_fraction", np.isnan(x).mean(), 0.01, 0.10)
            add(source, v, "fill_value_like", float((np.abs(x) > 1e4).any()), 0, 0,
                note="any |value| > 1e4")
            add(source, v, "flatline_fraction", longest_flat_fraction(x), 0.005, 0.05,
                note="share of steps in >=1 day runs of an identical non-zero value")
            lo, hi = BOUNDS[v]
            add(source, v, "out_of_bounds_fraction", np.nanmean((x < lo) | (x > hi)),
                0.001, 0.01, note=f"bounds {lo}..{hi}")

            # --- sign and size vs obs --------------------------------------
            ok = meas[v] & np.isfinite(x)
            if ok.sum() > 100:
                r = np.corrcoef(x[ok], obs[v][ok])[0, 1]
                add(source, v, "corr_with_obs", r, 0.3, 0.0, higher_is_bad=False)
            dm = np.nanmean(np.abs(x[day]))
            do = np.nanmean(np.abs(np.where(meas[v], obs[v], np.nan)[day]))
            if do > 0 and v != "Qg":
                ratio = dm / do
                add(source, v, "daytime_size_ratio_offscale", max(ratio, 1 / ratio) if ratio > 0 else np.inf,
                    2.0, 3.0, note=f"model/obs daytime mean |flux| = {ratio:.2f}")

            # --- spin-up: first year residual vs the rest ------------------
            res = pd.Series(np.where(meas[v], x - obs[v], np.nan)).groupby(months).mean()
            if len(res.dropna()) > 24:
                first, rest = res.iloc[:12].mean(), res.iloc[12:]
                z = abs(first - rest.mean()) / rest.std()
                add(source, v, "first_year_drift_z", z, 1.5, 3.0,
                    note="|first-year mean residual - later mean| / later monthly std")

        # --- physics inside the model -------------------------------------
        if "GPP" in have:
            add(source, "GPP", "night_GPP_mean", np.nanmean(m["GPP"][night]), 0.5, 2.0,
                note="photosynthesis in the dark (umol/m2/s)")
            add(source, "GPP", "negative_GPP_fraction", np.nanmean(m["GPP"] < -0.5), 0.01, 0.05)
        if "NEE" in have:
            add(source, "NEE", "daytime_NEE_mean", np.nanmean(m["NEE"][day]), -0.5, 1.0,
                note="should be negative (uptake) on average in daylight")
        if {"Rnet", "Qle", "Qh", "Qg"} <= set(have):
            resid = m["Rnet"] - m["Qle"] - m["Qh"] - m["Qg"]
            add(source, "Rnet", "energy_balance_residual_frac",
                abs(np.nanmean(resid)) / np.nanmean(np.abs(m["Rnet"])), 0.05, 0.15,
                note="|mean(Rnet-Qle-Qh-Qg)| / mean|Rnet|; models should close")
        if "Qle" in have:
            et = np.nansum(m["Qle"]) / LATENT * step_s   # mm over the record
            p = np.nansum(precip_mm)
            if p > 0:
                add(source, "Qle", "ET_over_rain", et / p, 1.2, 2.0,
                    note="whole-record evaporation / rainfall")
        if "Rnet" in have:
            ok = meas["Rnet"] & np.isfinite(m["Rnet"])
            bias = np.mean(m["Rnet"][ok] - obs["Rnet"][ok]) / np.mean(np.abs(obs["Rnet"][ok]))
            add(source, "Rnet", "Rnet_relative_bias", abs(bias), 0.15, 0.30,
                note=f"bias {bias:+.0%}; same forcing, so a big gap is suspicious")

        # --- skill (not pass/fail) ----------------------------------------
        for v in have:
            o = np.where(meas[v], obs[v], np.nan)
            ok = np.isfinite(o) & np.isfinite(m[v])
            if ok.sum() < 100:
                continue
            bench = one_lin(sw, o)
            metrics.append({
                "site": site, "source": source, "variable": v, "n": int(ok.sum()),
                "r": np.corrcoef(m[v][ok], o[ok])[0, 1],
                "bias": np.mean(m[v][ok] - o[ok]),
                "rmse": np.sqrt(np.mean((m[v][ok] - o[ok]) ** 2)),
                "nme": nme(m[v], o),
                "nme_1lin_insample": nme(bench, o),
            })

    # --- obs quality ----------------------------------------------------------
    q = {"site": site}
    for v in FLUXES:
        q[f"measured_{v}"] = round(float(meas[v].mean()), 3)
    q["LWdown_gapfilled"] = round(float((ds.LWdown_qc.values != 0).mean()), 3)
    all_meas = meas["Qle"] & meas["Qh"] & meas["Rnet"] & meas["Qg"] & day
    if all_meas.sum() > 100:
        q["energy_closure"] = round(float(
            (obs["Qle"][all_meas] + obs["Qh"][all_meas]).sum()
            / (obs["Rnet"][all_meas] - obs["Qg"][all_meas]).sum()), 3)
    q["VPD_file_ratio"] = ds.attrs.get("corrections", "").split("recomputed = ")[-1].split(" ")[0] \
        if "recomputed = " in ds.attrs.get("corrections", "") else None

    checks = pd.DataFrame(rows)
    metrics = pd.DataFrame(metrics)
    if not metrics.empty:
        metrics["beats_1lin"] = metrics.nme < metrics.nme_1lin_insample
        metrics = metrics.round(3)
    return checks, metrics, pd.DataFrame([q])


def validate(paths, alignment=None):
    out = [check_site(xr.open_dataset(p).load(), alignment) for p in paths]
    return tuple(pd.concat([o[i] for o in out], ignore_index=True) for i in range(3))

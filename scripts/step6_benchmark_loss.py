"""Find periods where a model does worse than every empirical benchmark.

Idea (Cranko Page et al. 2026, "LSM loss"): the benchmarks only see the weather.
If a land model is worse than all of them, the information was in the weather
but the model did not use it -- a sign of a missing or wrong process.

Run from the repo root:
    python scripts/step6_benchmark_loss.py                 # all sites
    python scripts/step6_benchmark_loss.py AU-How AU-Tum   # named sites

Settings: config/detection.yaml (benchmark_loss section).
Output:   outputs/detection/benchmark_loss_periods.csv, one row per period.

Steps, for each site, model and variable:
 1. Daily error. For each day, the mean absolute error (MAE) of the model and
    of each benchmark against the tower, using only steps the tower measured.
    Days with too few measured steps are dropped.
 2. Daily loss = model MAE - best benchmark MAE that day.
    Loss > 0 means the model was worse than every benchmark that day.
 3. Smooth the loss with a centred rolling mean (default 7 days), so one odd
    day does not make a period and one good day does not break one.
 4. Flag days whose smoothed loss is in the model's own worst 5% at this site
    (95th percentile) AND above a minimum size.
 5. Join flagged days into periods (bridging short gaps) and keep periods that
    are long enough.
 6. Keep a period only if its mean loss is several times the model's usual loss
    at that site (min_loss_ratio), so chronic bias alone can't create one.
 7. Keep a period only if the tower measured enough of it (min_period_measured).
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import yaml

HARMONISED = Path("data/harmonised")
OUT = Path("outputs/detection")
SKIP_SITES = ["AU-Ync"]  # copy of AU-ASM (its tower files are a copy of AU-ASM)

with open("config/detection.yaml") as fh:
    CFG = yaml.safe_load(fh)["benchmark_loss"]


def daily_mae(model, obs, measured):
    """Mean absolute error per day, over the steps the tower measured.
    model, obs: pandas Series on the half-hourly time axis."""
    error = (model - obs).abs().where(measured)
    return error.resample("1D").mean()


def find_runs(flags, max_gap, min_length):
    """Turn a True/False daily series into a list of (start, end) dates.
    Gaps of up to max_gap days are bridged; runs shorter than min_length are dropped."""
    days = flags[flags].index
    if len(days) == 0:
        return []
    runs = []
    start = previous = days[0]
    for day in days[1:]:
        if (day - previous).days > max_gap + 1:   # gap too long: close the run
            runs.append((start, previous))
            start = day
        previous = day
    runs.append((start, previous))
    return [(s, e) for s, e in runs if (e - s).days + 1 >= min_length]


def detect_site(path):
    site = path.stem
    ds = xr.open_dataset(path)
    rows = []

    for var in CFG["variables"]:
        obs = ds[var].sel(source="obs").to_series()
        measured = ds[f"measured_{var}"].to_series().astype(bool)

        # Days with enough measured steps
        measured_share = measured.resample("1D").mean()
        good_day = measured_share >= CFG["min_measured_fraction"]

        # Best benchmark per day
        bench_mae = {}
        for b in CFG["benchmarks"]:
            if b in ds.source.values and np.isfinite(ds[var].sel(source=b).values).any():
                bench_mae[b] = daily_mae(ds[var].sel(source=b).to_series(), obs, measured)
        if not bench_mae:
            continue
        bench_mae = pd.DataFrame(bench_mae)
        bench_mae = bench_mae.dropna(how="all")       # days with no benchmark values
        best_bench_mae = bench_mae.min(axis=1)
        best_bench_name = bench_mae.idxmin(axis=1)

        for model in CFG["models"]:
            if model not in ds.source.values:
                continue
            model_series = ds[var].sel(source=model).to_series()
            if not np.isfinite(model_series.values).any():
                continue                              # model doesn't report this variable

            model_mae = daily_mae(model_series, obs, measured)
            loss = (model_mae - best_bench_mae).where(good_day)

            # Smooth, then compare with this model's own distribution at this site
            smooth = loss.rolling(CFG["smooth_days"], center=True,
                                  min_periods=CFG["smooth_days"] // 2 + 1).mean()
            threshold = max(np.nanpercentile(smooth, CFG["percentile"]), CFG["min_loss"])
            flagged = smooth > threshold

            whole_record_loss = loss.mean()
            for start, end in find_runs(flagged, CFG["max_gap_days"], CFG["min_days"]):
                # Episode test: much worse than this model's usual loss at this site,
                # so a chronic bias alone can't make a period.
                if loss[start:end].mean() < CFG["min_loss_ratio"] * max(whole_record_loss, 1):
                    continue
                period = slice(start, end + pd.Timedelta(hours=23, minutes=59))
                m = measured[period]
                # Data-quality test: enough of the period measured by the tower
                if m.mean() < CFG["min_period_measured"]:
                    continue
                o = obs[period].where(m)
                md = model_series[period].where(m)
                rows.append({
                    "site": site,
                    "variable": var,
                    "model": model,
                    "start": start.date(),
                    "end": end.date(),
                    "n_days": (end - start).days + 1,
                    "mean_loss": round(loss[start:end].mean(), 2),
                    "peak_smoothed_loss": round(smooth[start:end].max(), 2),
                    "model_mae": round(model_mae[start:end].mean(), 2),
                    "best_benchmark_mae": round(best_bench_mae[start:end].mean(), 2),
                    "best_benchmark_usually": best_bench_name[start:end].mode().iloc[0],
                    "obs_mean": round(o.mean(), 2),
                    "model_mean": round(md.mean(), 2),
                    "model_bias": round((md - o).mean(), 2),
                    "measured_share": round(m.mean(), 2),
                    "site_threshold": round(threshold, 2),
                    # whole-record context, to tell chronic bias from an episode
                    "model_bias_whole_record": round((model_series - obs).where(measured).mean(), 2),
                    "model_loss_whole_record": round(loss.mean(), 2),
                })
    return rows


def main(args):
    sites = args or [p.stem for p in sorted(HARMONISED.glob("*.nc")) if p.stem not in SKIP_SITES]
    rows = []
    for site in sites:
        print(f"{site}")
        rows += detect_site(HARMONISED / f"{site}.nc")

    periods = pd.DataFrame(rows)
    # Record how the periods were found
    for key in ["min_measured_fraction", "smooth_days", "percentile", "min_loss",
                "max_gap_days", "min_days", "min_loss_ratio", "min_period_measured"]:
        periods[f"cfg_{key}"] = CFG[key]

    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "benchmark_loss_periods.csv"
    periods.to_csv(out, index=False)

    print(f"\n{len(periods)} periods -> {out}")
    if len(periods):
        print("\nPeriods per model and variable:")
        print(periods.groupby(["model", "variable"]).size().to_string())
        print("\nTen biggest by mean loss:")
        cols = ["site", "variable", "model", "start", "end", "n_days", "mean_loss",
                "model_bias", "model_bias_whole_record"]
        print(periods.sort_values("mean_loss", ascending=False)[cols].head(10).to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1:])

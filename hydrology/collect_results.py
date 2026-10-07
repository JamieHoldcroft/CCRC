"""Collect everything the paper needs into hydrology/results/.

Run from the repo root, or from anywhere:
    python hydrology/collect_results.py

Uses the functions in hydrology/et_signatures.py, so the signatures are exactly
the same as when you run that script. AU-Ync is skipped (its tower files are a
copy of AU-ASM).

Writes to hydrology/results/:
    et_sites.csv        one row per site x source: the Table 2 signatures, the
                        dry/wet-season ratio, and site metadata
    et_monthly_clim.csv mean monthly ET, rain and PET per site and source
    et_annual.csv       ET, PET and water stress (WS) for every complete year,
                        per site and source (for the WS time series)
    et_monthly.csv      the same for every complete month (monthly WS time series)
    et_eb.csv           tower energy gap and energy-balance-corrected tower ET
    drydown_curves.csv  ET against days since rain, per site and source
    drydown_sites.csv   dry-down summary per site and source (half-life etc.)
    event_composites.csv ET around big rain events that end a dry spell
    event_sites.csv     event summary per site and source (peak response etc.)
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))
import et_signatures as ets  # noqa: E402  (the Table 2 code)

HARMONISED = REPO / "data" / "harmonised"
RESULTS = HERE / "results"
SKIP_SITES = ["AU-Ync"]          # its tower files are a copy of AU-ASM

# Dry-down settings: a "rain day" has at least this much rain (mm/day)
DRYDOWN_RAIN_MM = 1.0
DRYDOWN_MAX_DAYS = 150           # follow dry spells up to this many days
DRYDOWN_MIN_SAMPLES = 8          # need this many days at a given day-since-rain

# Event composite settings: a big rain day that ends a dry spell
EVENT_RAIN_MM = 10.0             # rain on the event day (mm)
EVENT_DRY_BEFORE_DAYS = 7        # previous days must each have < 1 mm
EVENT_BEFORE, EVENT_AFTER = 5, 20  # days shown before and after the event


def days_since_rain(rain_daily):
    """For each day, how many days since the last day with >= DRYDOWN_RAIN_MM.
    0 on a rain day; NaN before the first rain day in the record."""
    out = np.full(len(rain_daily), np.nan)
    count = np.nan
    for i, r in enumerate(rain_daily.values):
        if np.isfinite(r) and r >= DRYDOWN_RAIN_MM:
            count = 0
        elif np.isfinite(count):
            count += 1
        out[i] = count
    return pd.Series(out, index=rain_daily.index)


def drydown(et_daily, since):
    """Median ET at each day since rain (1..DRYDOWN_MAX_DAYS), normalised by the
    median over days 1-3. Returns the curve and summary numbers."""
    df = pd.DataFrame({"et": et_daily, "d": since}).dropna()
    df = df[(df.d >= 1) & (df.d <= DRYDOWN_MAX_DAYS)]
    counts = df.groupby("d").et.count()
    med = df.groupby("d").et.median()[counts >= DRYDOWN_MIN_SAMPLES]
    if len(med) < 10 or med.loc[med.index <= 3].empty:
        return None, {}
    start = med.loc[med.index <= 3].mean()
    norm = med / start
    smooth = norm.rolling(7, center=True, min_periods=3).median()
    below = smooth[smooth <= 0.5]
    half_life = float(below.index[0]) if len(below) else np.nan
    # e-folding time from a straight-line fit of log(ET) over days 1-60
    fit = med[(med.index <= 60) & (med > 0)]
    tau = np.nan
    if len(fit) >= 15:
        slope = np.polyfit(fit.index.values, np.log(fit.values), 1)[0]
        tau = -1 / slope if slope < 0 else np.inf
    summary = {"et_days1to3": start, "half_life_days": half_life, "efold_days": tau,
               "max_day_available": int(med.index.max()),
               "ratio_day30": float(smooth.get(30, np.nan)),
               "ratio_day60": float(smooth.get(60, np.nan))}
    return pd.DataFrame({"day": norm.index, "et_median": med.values, "et_norm": norm.values,
                         "et_norm_smooth": smooth.values}), summary


def find_events(rain_daily):
    """Days with >= EVENT_RAIN_MM rain after EVENT_DRY_BEFORE_DAYS days with < 1 mm."""
    r = rain_daily.values
    events = []
    for i in range(EVENT_DRY_BEFORE_DAYS, len(r) - EVENT_AFTER):
        if r[i] >= EVENT_RAIN_MM and np.all(r[i - EVENT_DRY_BEFORE_DAYS:i] < 1.0):
            events.append(i)
    return events


def event_composite(et_daily, events):
    """Mean change in daily ET relative to the 5 days before each event."""
    lags = np.arange(-EVENT_BEFORE, EVENT_AFTER + 1)
    rows = []
    v = et_daily.values
    for i in events:
        if i - EVENT_BEFORE < 0:
            continue
        window = v[i - EVENT_BEFORE: i + EVENT_AFTER + 1]
        base = np.nanmean(v[i - EVENT_BEFORE: i])
        if not np.isfinite(base):
            continue
        rows.append(window - base)
    if not rows:
        return None
    return pd.DataFrame({"lag": lags, "d_et": np.nanmean(np.array(rows), axis=0),
                         "n_events": len(rows)})


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    site_rows, clim_rows, eb_rows, annual_rows, monthly_rows = [], [], [], [], []
    dd_curves, dd_rows, ev_curves, ev_rows = [], [], [], []

    for path in sorted(HARMONISED.glob("*.nc")):
        site = path.stem
        if site in SKIP_SITES:
            continue
        print(site)
        ds = xr.open_dataset(path)

        rain_daily = ds.Precip_mm.to_series().resample("1D").sum(min_count=1)
        pet_daily = ets.daily_pet(ds)
        pet_monthly = ets.totals(pet_daily, "MS")
        pet_annual = ets.totals(pet_daily, "YS")
        rain_monthly = ets.totals(rain_daily, "MS")
        rain_clim = rain_monthly.groupby(rain_monthly.index.month).mean()
        pet_clim = pet_monthly.groupby(pet_monthly.index.month).mean()
        dry_months = rain_clim.nsmallest(3).index      # 3 calendar months with least rain
        wet_months = rain_clim.nlargest(3).index       # 3 with most
        since = days_since_rain(rain_daily)
        events = find_events(rain_daily)

        meta = {
            "site": site, "site_name": ds.attrs.get("site_name"),
            "veg": ds.attrs.get("IGBP_veg_short"), "veg_long": ds.attrs.get("IGBP_veg_long"),
            "lat": ds.attrs.get("latitude"), "lon": ds.attrs.get("longitude"),
            "elevation": ds.attrs.get("elevation"), "canopy_height": ds.attrs.get("canopy_height"),
            "rain_mm_yr": ets.totals(rain_daily, "YS").median(),
            "pet_mm_yr": ets.totals(pet_daily, "YS").median(),
            "tair_mean": float(ds.Tair_C.mean()), "vpd_mean": float(ds.VPD_kPa.mean()),
            "lai_mean": float(ds.LAI.mean()) if "LAI" in ds else np.nan,
            "frac_qle_measured": float(ds.measured_Qle.mean()),
            "dry_months": " ".join(str(m) for m in sorted(dry_months)),
            "wet_months": " ".join(str(m) for m in sorted(wet_months)),
            "n_dry_events": len(events),
        }

        for source in ets.SOURCES:
            aet_daily = ets.daily_et(ds, source)
            aet_monthly = ets.totals(aet_daily, "MS")
            aet_annual = ets.totals(aet_daily, "YS")
            aet_clim = aet_monthly.groupby(aet_monthly.index.month).mean()
            r, n_events = ets.event_responsiveness(aet_daily, rain_daily)
            site_rows.append({**meta, "source": source,
                              "years": len(aet_annual), "n_events": n_events,
                              "ET_median": ets.long_term_median(aet_annual),
                              "CV_annual": ets.cv_annual(aet_annual),
                              "TSP": ets.tsp(aet_monthly),
                              "CV_monthly": ets.cv_monthly(aet_monthly),
                              "WS": ets.water_stress(aet_monthly, pet_monthly),
                              "R": r,
                              "ET_driest3": aet_clim[dry_months].mean(),
                              "ET_wettest3": aet_clim[wet_months].mean(),
                              "dry_wet_ratio": aet_clim[dry_months].mean() / aet_clim[wet_months].mean()})
            # Water stress year by year: the same formula as signature 6 (WS),
            # but using one year's ET and PET instead of the whole record
            for year_start, aet_year in aet_annual.items():
                if year_start not in pet_annual.index:
                    continue                    # PET incomplete that year
                pet_year = pet_annual[year_start]
                annual_rows.append({"site": site, "source": source, "year": year_start.year,
                                    "ET": aet_year, "PET": pet_year,
                                    "WS": (pet_year - aet_year) / pet_year})
            # ... and month by month
            for month_start, aet_month in aet_monthly.items():
                if month_start not in pet_monthly.index:
                    continue                    # PET incomplete that month
                pet_month = pet_monthly[month_start]
                monthly_rows.append({"site": site, "source": source, "month": month_start,
                                     "ET": aet_month, "PET": pet_month,
                                     "WS": (pet_month - aet_month) / pet_month})
            for m in range(1, 13):
                clim_rows.append({"site": site, "source": source, "month": m,
                                  "ET": aet_clim.get(m, np.nan), "rain": rain_clim.get(m, np.nan),
                                  "PET": pet_clim.get(m, np.nan)})

            curve, summary = drydown(aet_daily, since)
            if curve is not None:
                dd_curves.append(curve.assign(site=site, source=source))
                dd_rows.append({"site": site, "source": source, **summary})

            comp = event_composite(aet_daily, events)
            if comp is not None:
                ev_curves.append(comp.assign(site=site, source=source))
                after = comp[comp.lag.between(0, 5)]
                peak = after.d_et.max()
                peak_lag = int(after.loc[after.d_et.idxmax(), "lag"])
                late = comp[comp.lag.between(10, 20)].d_et.mean()
                ev_rows.append({"site": site, "source": source, "n_events": int(comp.n_events.iloc[0]),
                                "peak_d_et": peak, "peak_lag": peak_lag, "late_d_et": late,
                                "persistence": late / peak if peak > 0 else np.nan})

        # Tower energy gap and energy-balance-corrected tower ET
        measured = (ds.measured_Qle.astype(bool) & ds.measured_Qh.astype(bool)).values
        raw = (ds.Qle.sel(source="obs") + ds.Qh.sel(source="obs")).values[measured]
        corrected = (ds.Qle_cor + ds.Qh_cor).values[measured]
        step = ds.attrs["timestep_s"]
        lam = ets.latent_heat(ds.Tair_C.to_series())
        et_cor = (ds.Qle_cor.to_series() / lam * step).resample("1D").sum(min_count=1)
        et_cor = et_cor.where(et_cor >= 0)
        eb_rows.append({"site": site,
                        "eb_gap_pct": 100 * (np.nanmean(corrected) / np.nanmean(raw) - 1),
                        "ET_obs_cor": ets.totals(et_cor, "YS").median()})

    pd.DataFrame(site_rows).to_csv(RESULTS / "et_sites.csv", index=False)
    pd.DataFrame(clim_rows).to_csv(RESULTS / "et_monthly_clim.csv", index=False)
    pd.DataFrame(eb_rows).to_csv(RESULTS / "et_eb.csv", index=False)
    pd.DataFrame(annual_rows).to_csv(RESULTS / "et_annual.csv", index=False)
    pd.DataFrame(monthly_rows).to_csv(RESULTS / "et_monthly.csv", index=False)
    pd.concat(dd_curves).to_csv(RESULTS / "drydown_curves.csv", index=False)
    pd.DataFrame(dd_rows).to_csv(RESULTS / "drydown_sites.csv", index=False)
    pd.concat(ev_curves).to_csv(RESULTS / "event_composites.csv", index=False)
    pd.DataFrame(ev_rows).to_csv(RESULTS / "event_sites.csv", index=False)
    print(f"wrote results to {RESULTS}")


if __name__ == "__main__":
    main()

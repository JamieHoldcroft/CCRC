"""Flag drought, heatwave and high-VPD periods at a flux-tower site.

Every flag is computed on daily values, then copied back to each half-hourly
timestep of that day, so the result is a boolean variable along `time`:

    ds.Qle.where(ds.drought)   # latent heat during drought only
    ds.where(ds.normal)        # every variable outside all events

Definitions (all thresholds are function arguments):
    drought   rolling 90-day rainfall below the 10th percentile for that time of
              year (all 90-day totals within +/-15 days of the same calendar
              day, pooled over every year), kept only in runs of >= 14 days
    heatwave  Excess Heat Factor (Nairn & Fawcett 2015, used by the Bureau):
              a 3-day mean Tair above the record's 95th percentile of daily
              mean Tair. Every day of that 3-day window is flagged. The
              `heatwave_ehf` variable adds severity, which grows when the
              3 days are also hotter than the 30 days before them.
    high_vpd  daily daytime-mean VPD (SWdown > 50 W/m2) above its 90th
              percentile over the record
    normal    none of the above

Flags depend only on the met forcing, so they are identical for obs and models.
"""

import numpy as np
import pandas as pd
import xarray as xr


def _daily(da, how):
    """Daily series (pandas) from a half-hourly DataArray."""
    s = da.to_series()
    return getattr(s.resample("1D"), how)()


def _runs_at_least(flags, min_days):
    """Keep only runs of True that last at least `min_days` consecutive days."""
    run_id = (flags != flags.shift()).cumsum()
    run_len = flags.groupby(run_id).transform("size")
    return flags & (run_len >= min_days)


def _to_timesteps(daily_flags, time):
    """Copy a daily boolean series onto every timestep of `time`."""
    days = pd.DatetimeIndex(time.values).floor("D")
    values = daily_flags.reindex(days).fillna(False).astype(bool).values
    return xr.DataArray(values, coords={"time": time}, dims="time")


def _seasonal_percentile(series, percentile, half_window_days):
    """Percentile of `series` for each calendar day, pooling all years.

    For calendar day d, uses every value whose day-of-year is within
    +/-half_window_days of d (wrapping around the new year).
    Returns a Series aligned with `series`.
    """
    doy = series.index.dayofyear.values
    values = series.values
    by_doy = np.full(367, np.nan)
    for d in range(1, 367):
        dist = np.abs(doy - d)
        dist = np.minimum(dist, 366 - dist)
        by_doy[d] = np.nanpercentile(values[dist <= half_window_days], percentile)
    return pd.Series(by_doy[doy], index=series.index)


def drought_flag(ds, window_days=90, percentile=10, min_days=14, half_window_days=15):
    """Rainfall-based drought: unusually low rainfall *for the time of year*.

    Comparing against the same season stops the normal dry season being
    flagged every year.
    """
    rain = _daily(ds.Precip_mm, "sum")
    total = rain.rolling(window_days, min_periods=window_days).sum()
    threshold = _seasonal_percentile(total, percentile, half_window_days)
    flags = _runs_at_least(total < threshold, min_days)
    out = _to_timesteps(flags, ds.time)
    out.attrs = {
        "long_name": "Drought (low rolling rainfall for the season)",
        "definition": f"{window_days}-day rainfall < seasonal p{percentile} "
                      f"(+/-{half_window_days} days, {threshold.min():.0f} to "
                      f"{threshold.max():.0f} mm), runs >= {min_days} days",
    }
    return out


def excess_heat_factor(ds, percentile=95, accl_days=30):
    """Daily Excess Heat Factor (Nairn & Fawcett 2015), indexed by the first
    day of each 3-day window.

    EHI_sig  = T3 - T_p95    how hot vs the whole record (fixed threshold)
    EHI_accl = T3 - T30      how hot vs the 30 days before (acclimatisation)
    EHF      = EHI_sig * max(1, EHI_accl)

    T3 is the mean daily Tair over days i, i+1, i+2. Positive EHF = heatwave.
    Because max(1, EHI_accl) >= 1, EHF > 0 exactly when EHI_sig > 0; the
    30-day term only scales severity. Nairn & Fawcett use 9am-9am Tmax/Tmin
    averages; here the midnight-midnight mean of half-hourly Tair stands in.
    """
    tmean = _daily(ds.Tair_C, "mean")
    t95 = np.nanpercentile(tmean, percentile)
    t3 = tmean[::-1].rolling(3, min_periods=3).mean()[::-1]  # days i..i+2
    t30 = tmean.shift(1).rolling(accl_days, min_periods=accl_days).mean()  # days i-30..i-1
    ehi_sig = t3 - t95
    ehi_accl = t3 - t30
    return ehi_sig * np.maximum(1, ehi_accl), t95


def heatwave_flag(ds, percentile=95, accl_days=30):
    """Heatwave days from the Excess Heat Factor: all 3 days of any window
    with EHF > 0 are flagged. Returns (flag, severity) on the half-hourly axis.
    """
    ehf, t95 = excess_heat_factor(ds, percentile, accl_days)
    starts = (ehf > 0).astype(int)
    # a day is in a heatwave if a positive window started on it or 1-2 days before
    flags = starts.rolling(3, min_periods=1).max().astype(bool)
    severity = ehf.clip(lower=0).rolling(3, min_periods=1).max().where(flags, 0.0)

    out = _to_timesteps(flags, ds.time)
    out.attrs = {
        "long_name": "Heatwave (Excess Heat Factor > 0)",
        "definition": f"3-day mean Tair > p{percentile} of daily mean Tair "
                      f"({t95:.1f} degC), EHF with {accl_days}-day acclimatisation",
    }
    days = pd.DatetimeIndex(ds.time.values).floor("D")
    sev = xr.DataArray(severity.reindex(days).fillna(0.0).values,
                       coords={"time": ds.time}, dims="time")
    sev.attrs = {"long_name": "Excess Heat Factor (heatwave severity)",
                 "units": "degC2"}
    return out, sev


def high_vpd_flag(ds, percentile=90, sw_threshold=50):
    """High VPD: days whose daytime-mean VPD is in the top of the record."""
    vpd_day = _daily(ds.VPD_kPa.where(ds.SWdown > sw_threshold), "mean")
    threshold = np.nanpercentile(vpd_day, percentile)
    out = _to_timesteps(vpd_day > threshold, ds.time)
    out.attrs = {
        "long_name": "High VPD",
        "definition": f"daytime-mean VPD > p{percentile} ({threshold:.2f} kPa)",
    }
    return out


def add_event_flags(ds, drought_kw=None, heatwave_kw=None, high_vpd_kw=None):
    """Return a copy of `ds` with drought, heatwave, high_vpd and normal added."""
    ds = ds.copy()
    ds["drought"] = drought_flag(ds, **(drought_kw or {}))
    ds["heatwave"], ds["heatwave_ehf"] = heatwave_flag(ds, **(heatwave_kw or {}))
    ds["high_vpd"] = high_vpd_flag(ds, **(high_vpd_kw or {}))
    ds["normal"] = ~(ds.drought | ds.heatwave | ds.high_vpd)
    ds["normal"].attrs = {"long_name": "No event flagged"}
    return ds


def summarise_events(ds):
    """Days flagged and number of separate events, per flag."""
    rows = {}
    for name in ["drought", "heatwave", "high_vpd", "normal"]:
        daily = ds[name].to_series().resample("1D").first()
        starts = daily & ~daily.shift(fill_value=False)
        rows[name] = {"days": int(daily.sum()), "events": int(starts.sum()),
                      "fraction": round(float(daily.mean()), 3)}
    return pd.DataFrame(rows).T

"""Evapotranspiration (ET) signatures for the tower (obs), CABLE and JULES."""

from pathlib import Path
import numpy as np
import pandas as pd

# ============================================================================
# Settings
# ============================================================================
HARMONISED = Path(__file__).resolve().parents[1] / "data" / "harmonised"   # CCRC/data/harmonised
SOURCES = ["obs", "CABLE", "JULES"]

PT_ALPHA = 1.26                  # Priestley-Taylor coefficient
EVENT_RAIN_MM = 10.0             # rain day threshold for R (mm/day)
EVENT_WINDOW_DAYS = 10           # look for max ET over days i ... i+10
MIN_COMPLETE = 0.001               # share of days needed for a month or year to count
MIN_MEASURED_PER_DAY = 0.001       # a day counts if >= x% of its steps were measured (qc < 20)


def latent_heat(tair_c):
    """Latent heat of vaporisation (J/kg) at air temperature tair_c (degC)."""
    return (2.501 - 0.00237 * tair_c) * 1e6


# ============================================================================
# Step 0: keeping measured tower data only
# ============================================================================
def is_measured(ds, var):
    """True for each step where the tower measured `var`, False where it was
    gap-filled. OzFlux qc flags: 0 = measured, 1-19 = QA checks and corrections
    to measured data (AU-DaP rain is flagged 10), 20 and up = gap-filling.
    var is e.g. "Qle", "Precip", "Rnet"."""
    return (ds[f"{var}_qc"] < 20).to_series()


def measured_daily(per_step, measured, how="sum"):
    """Daily values from measured steps only.
    1. Steps where measured is False (gap-filled) are dropped.
    2. Days where fewer than MIN_MEASURED_PER_DAY of the steps were measured
       are dropped (NaN).
    how="sum" adds up the kept steps of each day (ET, rain, radiation);
    how="mean" averages them (temperature, pressure)."""
    kept = per_step.where(measured)

    """   
    # Count gap-filled runs of 4 steps or fewer (2 h at half-hourly sites) as measured
    run_id = (measured != measured.shift()).cumsum()            # give each run of T's or F's its own number
    run_length = measured.groupby(run_id).transform("size")     # how long the run each step is in
    measured = measured | (run_length <= 4)                     # short runs become True
    
    """
    if how == "sum":
        daily = kept.resample("1D").sum(min_count=1)
    else:
        daily = kept.resample("1D").mean()
    measured_share = measured.astype(float).resample("1D").mean()   # 0-1 per day
    return daily.where(measured_share >= MIN_MEASURED_PER_DAY)


# ============================================================================
# Step 1: daily ET (mm/day) for one source
# This is the ONLY place ET is made. collect_results.py calls it too, so any
# change here (e.g which tower data are used) reaches every output
# ============================================================================


def daily_et(ds, source):
    step_seconds = ds.attrs["timestep_s"]
    qle = ds.Qle.sel(source=source).to_series()          # W/m2, half-hourly (or hourly)
    lam = latent_heat(ds.Tair_C.to_series())             # J/kg, for each step's temperature
    et_per_step = qle / lam * step_seconds               # mm per step (1 kg/m2 of water = 1 mm)
    # The TOWER's Qle mask for every source: models use the same steps as obs
    et_daily = measured_daily(et_per_step, is_measured(ds, "Qle"))
    et_daily = et_daily.where(et_daily >= 0)             # drop negative days (paper excludes negative AET)
    return et_daily


# ============================================================================
# Step 2: daily rain (mm/day), measured steps only
# ============================================================================
def daily_rain(ds):
    return measured_daily(ds.Precip_mm.to_series(), is_measured(ds, "Precip"), "sum")


# ============================================================================
# Step 2b: daily PET (mm/day), Priestley-Taylor, from the tower weather
# ============================================================================
def daily_pet(ds):
    step_seconds = ds.attrs["timestep_s"]
    rnet = ds.Rnet.sel(source="obs").to_series()          # W/m2, tower
    tair = ds.Tair_C.to_series()                          # degC
    psurf_kpa = ds.Psurf.to_series() / 1000               # Pa -> kPa

    # The day's net radiation in J/m2 and the day's mean temperature and
    # pressure, each from its own measured steps only. If any of the three is
    # missing that day, PET is missing that day.
    rnet_day = measured_daily(rnet * step_seconds, is_measured(ds, "Rnet"), "sum").clip(lower=0)
    t_day = measured_daily(tair, is_measured(ds, "Tair"), "mean")
    p_day = measured_daily(psurf_kpa, is_measured(ds, "Psurf"), "mean")

    # Priestley-Taylor: PET = alpha * delta / (delta + gamma) * Rnet / lambda
    es = 0.6108 * np.exp(17.27 * t_day / (t_day + 237.3))       # kPa, saturation vapour pressure
    delta = 4098 * es / (t_day + 237.3) ** 2                    # kPa/degC, slope of es
    gamma = 0.000665 * p_day                                    # kPa/degC, psychrometric constant
    pet = PT_ALPHA * delta / (delta + gamma) * rnet_day / latent_heat(t_day)   # mm/day
    return pet


# ============================================================================
# Step 2c: put ET (every source) and PET on exactly the same days
# ============================================================================
def analysis_days(ds):
    """Daily ET for obs, CABLE and JULES, and daily PET, all on the same days.
    A day is kept only if every source's ET AND PET have a value that day
    (a source can lose a day on its own, e.g. a negative daily total).
    Returns ({source: daily ET}, daily PET)."""
    et = {source: daily_et(ds, source) for source in SOURCES}
    pet = daily_pet(ds)
    everything = pd.concat(list(et.values()) + [pet], axis=1)    # one column each, same dates
    keep = everything.notna().all(axis=1)                       # True = all have a value
    return {source: e.where(keep) for source, e in et.items()}, pet.where(keep)


# ============================================================================
# Step 3: monthly and annual totals, keeping only complete months / years
# ============================================================================
def totals(daily, freq):
    """Monthly ('MS') or annual ('YS') totals; drop incomplete months / years.
    Total = mean of the days present x days in the period, so a month with
    missing days is not counted as smaller than a full one."""
    total = daily.resample(freq).mean() * daily.resample(freq).size()   # mean day x days in period
    days_present = daily.notna().resample(freq).sum()
    days_in_period = daily.resample(freq).size()
    complete = days_present >= MIN_COMPLETE * days_in_period
    return total.where(complete).dropna()


# ============================================================================
# Step 4: the signatures (numbering as in Table 2 of the paper)
# ============================================================================
def long_term_median(annual):
    return annual.median()


def cv_annual(annual):
    return annual.std() / annual.mean()


def tsp(monthly):
    if monthly.empty:                    # no complete month at all
        return np.nan
    median_by_month = monthly.groupby(monthly.index.month).median()
    return int(median_by_month.idxmax())


def cv_monthly(monthly):
    return monthly.std() / monthly.mean()


def water_stress(monthly_aet, monthly_pet):
    pet_mean = monthly_pet.mean()
    aet_mean = monthly_aet.mean()
    return (pet_mean - aet_mean) / pet_mean


def event_responsiveness(daily_aet, daily_rain):
    event_days = daily_rain[daily_rain > EVENT_RAIN_MM].index
    rain_values = []
    max_aet_values = []
    for day in event_days:
        window = daily_aet[day: day + pd.Timedelta(days=EVENT_WINDOW_DAYS)]
        if window.notna().sum() == 0:
            continue
        rain_values.append(daily_rain[day])
        max_aet_values.append(window.max())
    if len(rain_values) < 5:
        return np.nan, len(rain_values)
    rain_anomaly = np.array(rain_values) - np.mean(rain_values)
    aet_anomaly = np.array(max_aet_values) - np.mean(max_aet_values)
    r = np.corrcoef(rain_anomaly, aet_anomaly)[0, 1]
    return r, len(rain_values)
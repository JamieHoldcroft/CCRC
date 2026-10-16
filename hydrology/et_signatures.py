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
MIN_COMPLETE = 0.7               # share of days needed for a month or year to count
MIN_MEASURED_PER_DAY = 0.8       # a day counts if >= x% of its steps were measured (qc < 20)


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
    return (ds[f"{var}_qc"] == 0 ).to_series()


DAY_SWDOWN = 5             # W/m2: a step is daytime if SWdown >= this
LONG_DAY_GAP_HOURS = 3     # a day is dropped if it has a daytime gap longer than this

def measured_daily(per_step, measured, daytime, how="sum"):
    """Daily values from ALL steps (gap-filled included), but a day is dropped
    (NaN) if it has a run of gap-filled DAYTIME steps longer than
    LONG_DAY_GAP_HOURS. Gaps at night are ignored (a night step ends a run).
    daytime: True/False per step, True when SWdown >= DAY_SWDOWN.
    how="sum" adds up each day's steps (ET, rain, radiation);
    how="mean" averages them (temperature, pressure)."""
    gapfilled = ~measured                                          # True = gap-filled step
    day_gap = gapfilled & daytime                                  # True = gap-filled AND in daylight

    # Find runs of daytime gap-filled steps and how long each run is (in hours)
    run_id = (day_gap != day_gap.shift()).cumsum()                 # number each run of T's or F's
    run_length = day_gap.groupby(run_id).transform("size")         # length of the run each step is in (steps)
    step_hours = (per_step.index[1] - per_step.index[0]).total_seconds() / 3600   # 0.5 or 1
    run_hours = run_length * step_hours                            # length of the run in hours
    in_long_day_gap = day_gap & (run_hours > LONG_DAY_GAP_HOURS)   # True on every step of a long daytime gap

    # A day is bad if any of its steps is part of a long daytime gap
    bad_day = in_long_day_gap.resample("1D").max()                 # True/False per day

    if how == "sum":
        daily = per_step.resample("1D").sum(min_count=1)
    else:
        daily = per_step.resample("1D").mean()
    return daily.where(~bad_day)                                    # NaN on bad days

# def measured_daily(per_step, measured, how="sum"):
#     """Daily values from measured steps only.
#     1. Steps where measured is False (gap-filled) are dropped.
#     2. Days where fewer than MIN_MEASURED_PER_DAY of the steps were measured
#        are dropped (NaN).
#     how="sum" adds up the kept steps of each day (ET, rain, radiation);
#     how="mean" averages them (temperature, pressure)."""
       
#     # Count gap-filled runs of 6 steps or fewer (3 h at half-hourly sites) as measured
#     run_id = (measured != measured.shift()).cumsum()            # give each run of T's or F's its own number
#     run_length = measured.groupby(run_id).transform("size")     # how long the run each step is in
#     measured = measured | (run_length <= 6)                     # short runs become True

#     kept = per_step.where(measured)

#     if how == "sum":
#         daily = kept.resample("1D").sum(min_count=1)
#     else:
#         daily = kept.resample("1D").mean()
#     measured_share = measured.astype(float).resample("1D").mean()   # 0-1 per day
#     return daily.where(measured_share >= MIN_MEASURED_PER_DAY)


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
    daytime = ds.SWdown.to_series() >= DAY_SWDOWN         # True in daylight
    et_daily = measured_daily(et_per_step, is_measured(ds, "Qle"), daytime)
    et_daily = et_daily.where(et_daily >= 0)             # drop negative days (paper excludes negative AET)
    return et_daily


# ============================================================================
# Step 2: daily rain (mm/day), measured steps only
# ============================================================================
def daily_rain(ds):
    daytime = ds.SWdown.to_series() >= DAY_SWDOWN         # True in daylight
    return measured_daily(ds.Precip_mm.to_series(), is_measured(ds, "Precip"), daytime, "sum")


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
    daytime = ds.SWdown.to_series() >= DAY_SWDOWN         # True in daylight
    rnet_day = measured_daily(rnet * step_seconds, is_measured(ds, "Rnet"), daytime, "sum")
    t_day = measured_daily(tair, is_measured(ds, "Tair"), daytime, "mean")
    p_day = measured_daily(psurf_kpa, is_measured(ds, "Psurf"), daytime, "mean")

    # Priestley-Taylor: PET = alpha * delta / (delta + gamma) * Rnet / lambda
    es = 0.6108 * np.exp(17.27 * t_day / (t_day + 237.3))       # kPa, saturation vapour pressure
    delta = 4098 * es / (t_day + 237.3) ** 2                    # kPa/degC, slope of es
    gamma = 0.000665 * p_day                                    # kPa/degC, psychrometric constant
    pet = PT_ALPHA * delta / (delta + gamma) * rnet_day / latent_heat(t_day)   # mm/day
    return pet


# FAO-56 Penman-Monteith reference ET (mm/day). To use it, comment out the
# daily_pet above and un-comment this one. Same filtering as above; it also
# needs daily mean VPD and wind.
# References (Allen, Pereira, Raes & Smith 1998, FAO Irrigation and Drainage
# Paper 56, https://www.fao.org/4/x0490e/x0490e00.htm):
#   eq. 6  (p. 24) the daily Penman-Monteith equation used below
#   eq. 8  (p. 32) gamma = 0.000665 * P
#   eq. 11 (p. 36) es(T); eq. 13 (p. 37) delta
#   eq. 42 (p. 54) G is about 0 for daily steps, so it is left out
#   eq. 47 (p. 56) wind at height z converted to wind at 2 m
# 0.408 = 1 / 2.45, turning MJ/m2/day of energy into mm/day of water.
#
# def daily_pet(ds):
#     step_seconds = ds.attrs["timestep_s"]
#     rnet = ds.Rnet.sel(source="obs").to_series()          # W/m2, tower
#     tair = ds.Tair_C.to_series()                          # degC
#     psurf_kpa = ds.Psurf.to_series() / 1000               # Pa -> kPa
#     vpd = ds.VPD_kPa.to_series()                          # kPa, vapour pressure deficit (es - ea)
#     wind = ds.Wind.to_series()                            # m/s, at the tower's measurement height
#
#     daytime = ds.SWdown.to_series() >= DAY_SWDOWN         # True in daylight
#     rnet_day = measured_daily(rnet * step_seconds / 1e6, is_measured(ds, "Rnet"), daytime, "sum")   # MJ/m2/day
#     t_day = measured_daily(tair, is_measured(ds, "Tair"), daytime, "mean")
#     p_day = measured_daily(psurf_kpa, is_measured(ds, "Psurf"), daytime, "mean")
#     vpd_day = measured_daily(vpd, is_measured(ds, "VPD"), daytime, "mean")
#     wind_day = measured_daily(wind, is_measured(ds, "Wind"), daytime, "mean")
#
#     # Wind at 2 m from wind at height z (eq. 47)
#     z = float(ds.attrs["reference_height"])               # m
#     u2 = wind_day * 4.87 / np.log(67.8 * z - 5.42)
#
#     es = 0.6108 * np.exp(17.27 * t_day / (t_day + 237.3))       # kPa, saturation vapour pressure (eq. 11)
#     delta = 4098 * es / (t_day + 237.3) ** 2                    # kPa/degC, slope of es (eq. 13)
#     gamma = 0.000665 * p_day                                    # kPa/degC, psychrometric constant (eq. 8)
#
#     # eq. 6: radiation term + dry-air-and-wind term, over delta + gamma * (1 + 0.34 u2)
#     radiation_term = 0.408 * delta * rnet_day
#     aero_term = gamma * 900 / (t_day + 273) * u2 * vpd_day
#     pet = (radiation_term + aero_term) / (delta + gamma * (1 + 0.34 * u2))   # mm/day
#     return pet


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
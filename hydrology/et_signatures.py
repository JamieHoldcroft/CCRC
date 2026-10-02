"""Evapotranspiration (ET) signatures for the tower (obs), CABLE and JULES.

Metrics from Table 2 (p. 6) of:
    Gardiya Weligamage et al. (2025), "Using evapotranspiration signatures to
    assess evapotranspiration realism in rainfall-runoff models",
    EGUsphere preprint, doi:10.5194/egusphere-2025-3373.
The signatures themselves are defined in Gardiya Weligamage et al. (2025),
HESSD, doi:10.5194/hess-2024-373.

Run from the repo root (or from anywhere; paths are relative to this file):
    python hydrology/et_signatures.py            # all sites
    python hydrology/et_signatures.py AU-How     # one or more named sites

Everything below is written step by step on purpose, so each number can be
traced back by hand. Nothing is saved; results are printed.

CHOICES THAT ARE OURS, NOT THE PAPER'S (change them in the settings below):
  1. AET from latent heat: ET (mm per step) = Qle / lambda * seconds per step,
     with lambda = (2.501 - 0.00237 * Tair) * 1e6 J/kg (Tair in degC), the
     latent heat of vaporisation at that step's air temperature. The paper uses
     OzFlux ET directly, so it gives no conversion.
  2. Obs: we use the gap-filled tower Qle, as the paper does (gap-filled OzFlux
     ET, v2024_V1, p. 6 lines 122-126), so annual and monthly totals are not
     biased by missing steps. CABLE and JULES have no gaps, so all three use
     the same time steps.
  3. Negative ET: the paper drops negative AET values. Half-hourly Qle is often
     slightly negative at night (dew), so we drop negative DAILY totals instead.
  4. PET: the paper uses Morton wet-environment PET from SILO, which we don't
     have. As a SUBSTITUTE we use Priestley-Taylor PET from the tower's net
     radiation and air temperature (alpha = 1.26). PET describes the weather,
     not the model, so obs, CABLE and JULES all get the same PET.
  5. Event responsiveness (R): the paper's rainfall threshold is not in Table 2.
     We use 10 mm/day.
  6. Completeness: a month counts only if at least 90% of its days are present,
     a year only if at least 90% of its days are present.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

# Settings
HARMONISED = Path(__file__).resolve().parents[1] / "data" / "harmonised"   # CCRC/data/harmonised
SOURCES = ["obs", "CABLE", "JULES"]

PT_ALPHA = 1.26                  # Priestley-Taylor coefficient
EVENT_RAIN_MM = 10.0             # rain day threshold for R (mm/day)
EVENT_WINDOW_DAYS = 10           # look for max ET over days i ... i+10
MIN_COMPLETE = 0.9               # percentage of days needed for a month or year to count


def latent_heat(tair_c):
    return (2.501 - 0.00237 * tair_c) * 1e6


# Step 1: daily ET (mm/day) for one source
def daily_et(ds, source):
    step_seconds = ds.attrs["timestep_s"]
    qle = ds.Qle.sel(source=source).to_series()          # W/m2, half-hourly (or hourly)
    lam = latent_heat(ds.Tair_C.to_series())             # J/kg, for each step's temperature
    et_per_step = qle / lam * step_seconds               # mm per step (1 kg/m2 of water = 1 mm)
    et_daily = et_per_step.resample("1D").sum(min_count=1)
    et_daily = et_daily.where(et_daily >= 0)             # drop negative days (paper excludes negative AET)
    return et_daily


# Step 2: daily PET (mm/day), Priestley-Taylor, from the tower weather
def daily_pet(ds):
    step_seconds = ds.attrs["timestep_s"]
    rnet = ds.Rnet.sel(source="obs").to_series()          # W/m2 (gap-filled tower)
    tair = ds.Tair_C.to_series()                          # degC
    psurf_kpa = ds.Psurf.to_series() / 1000               # Pa -> kPa

    # Daily means of the weather, and the day's net radiation in J/m2
    rnet_day = (rnet * step_seconds).resample("1D").sum(min_count=1).clip(lower=0)
    t_day = tair.resample("1D").mean()
    p_day = psurf_kpa.resample("1D").mean()

    # Priestley-Taylor: PET = alpha * delta / (delta + gamma) * Rnet / lambda
    es = 0.6108 * np.exp(17.27 * t_day / (t_day + 237.3))       # kPa, saturation vapour pressure
    delta = 4098 * es / (t_day + 237.3) ** 2                    # kPa/degC, slope of es
    gamma = 0.000665 * p_day                                    # kPa/degC, psychrometric constant
    pet = PT_ALPHA * delta / (delta + gamma) * rnet_day / latent_heat(t_day)   # mm/day
    return pet


# ----------------------------------------------------------------------------
# Step 3: monthly and annual totals, keeping only complete months / years
# ----------------------------------------------------------------------------
def totals(daily, freq):
    """Sum daily values to months ('MS') or years ('YS'); drop incomplete ones."""
    total = daily.resample(freq).sum(min_count=1)
    days_present = daily.notna().resample(freq).sum()
    days_in_period = daily.resample(freq).size()
    complete = days_present >= MIN_COMPLETE * days_in_period
    return total.where(complete).dropna()


def long_term_median(annual):
    return annual.median()


def cv_annual(annual):
    return annual.std() / annual.mean()


def tsp(monthly):
    median_by_month = monthly.groupby(monthly.index.month).median()
    return int(median_by_month.idxmax())


def cv_monthly(monthly):
    """It measures how big the ups and downs of monthly ET are compared with the
    average month. 0 = every month evaporates the same amount; 0.5 = a typical
    month is about 50% above or below the average. A site with a strong wet and
    dry season has a high value and a site that evaporates steadily all year has a
    low one. Mostly it reflects the size of the seasonal cycle, plus some
    year-to-year variation (a drought year adds low months)."""
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


# Step 5: run it site by site and print a small table for each
def report_site(path):
    site = path.stem
    ds = xr.open_dataset(path)

    rain_daily = ds.Precip_mm.to_series().resample("1D").sum(min_count=1)
    pet_daily = daily_pet(ds)
    pet_monthly = totals(pet_daily, "MS")

    rows = {}
    years_used = None
    n_events = None
    for source in SOURCES:
        aet_daily = daily_et(ds, source)
        aet_monthly = totals(aet_daily, "MS")
        aet_annual = totals(aet_daily, "YS")
        years_used = len(aet_annual)
        r, n_events = event_responsiveness(aet_daily, rain_daily)

        rows[source] = {
            "1 Long-term median (mm/yr)": long_term_median(aet_annual),
            "2 CV_annual": cv_annual(aet_annual),
            "4 TSP (month 1-12)": tsp(aet_monthly),
            "5 CV_monthly": cv_monthly(aet_monthly),
            "6 WS (water stress)": water_stress(aet_monthly, pet_monthly),
            "8 R (event responsiveness)": r,
        }

    table = pd.DataFrame(rows)
    print("=" * 66)
    print(f"{site}  ({ds.attrs.get('IGBP_veg_short', '?')}, {years_used} complete years, "
          f"{n_events} rain events > {EVENT_RAIN_MM:g} mm)")
    print("-" * 66)
    print(table.round(4).to_string())
    if years_used < 3:
        print("  NOTE: fewer than 3 complete years, so 1 and 2 are unreliable here.")
    print()
    return table.assign(site=site)


def main(sites):
    if not sites:
        sites = [p.stem for p in sorted(HARMONISED.glob("*.nc"))]
    tables = []
    for site in sites:
        tables.append(report_site(HARMONISED / f"{site}.nc"))

    # Across all sites: how far CABLE and JULES are from the tower, per metric
    if len(tables) > 1:
        print("=" * 66)
        print("All sites: median of (model - obs) across sites")
        print("-" * 66)
        allsites = pd.concat(tables)
        allsites.index.name = "metric"
        gap = pd.DataFrame({
            "CABLE - obs": (allsites.CABLE - allsites.obs).groupby("metric").median(),
            "JULES - obs": (allsites.JULES - allsites.obs).groupby("metric").median(),
        })
        print(gap.round(4).to_string())
        print("-" * 66)


if __name__ == "__main__":
    main(sys.argv[1:])

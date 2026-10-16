# The hydrology pipeline, step by step as it runs

Written 12 Oct 2026 from the code on your machine at that time:

- `hydrology/et_signatures.py`
- `hydrology/collect_results.py`
- `hydrology/make_report.py`

The settings quoted are the ones in the files now, for example `MIN_COMPLETE = 0.9`. If you change a setting, the rules below change with it.

You run two commands, in this order:

```
python hydrology/collect_results.py     # data -> CSV tables in hydrology/results/
python hydrology/make_report.py         # CSV tables -> PNG figures in hydrology/figures/
```

`et_signatures.py` is never run on its own. It is a toolbox of functions that `collect_results.py` imports. `make_report.py` never touches the raw data; it only reads the CSVs.

---

## 0. The big picture in one paragraph

For every site, the pipeline does five things:

1. It opens the harmonised file: one NetCDF holding tower data plus every model, on the tower's half-hourly (or hourly) time axis.
2. It turns half-hourly Qle into daily ET, in mm/day, for the tower, CABLE and JULES.
3. It turns half-hourly rain into daily rain, and half-hourly Rnet, Tair and pressure into daily Priestley-Taylor PET.
4. It throws away days that rely on long daytime gap-filling, then keeps only days on which all three ET sources and PET have a value.
5. It adds the days up into months and years, keeping only months and years that are at least 90% complete, and computes the signatures from those.

Everything is written to CSV. `make_report.py` then reads the CSVs and draws the figures.

---

## PART A. `collect_results.py` starts

### A1. Imports and settings (lines 23-41)

- **`sys.path.insert(0, str(HERE))` and `import et_signatures as ets`** make the functions in `et_signatures.py` available as `ets.something`.
- **Importing `et_signatures.py` runs its top level once.** That sets the settings:

| setting | value | what it controls |
|---|---|---|
| `SOURCES` | `["obs", "CABLE", "JULES"]` | which sources get ET |
| `PT_ALPHA` | 1.26 | Priestley-Taylor coefficient |
| `EVENT_RAIN_MM` | 10 | a "rain event" is a day with more than 10 mm |
| `EVENT_WINDOW_DAYS` | 10 | look for peak ET from the rain day to 10 days after |
| `MIN_COMPLETE` | 0.9 | a month or year counts only if at least 90% of its days have a value |
| `MIN_MEASURED_PER_DAY` | 0.8 | **no longer used by anything** (left over from the old `measured_daily`) |
| `DAY_SWDOWN` | 5 W/m² | a half-hour is "daytime" if SWdown ≥ 5 |
| `LONG_DAY_GAP_HOURS` | 3 | a day is thrown out if it has a daytime gap-filled run longer than 3 h |

- **In `collect_results.py` itself:**
  - `HARMONISED` = `data/harmonised/`
  - `RESULTS` = `hydrology/results/`
  - `SKIP_SITES = ["AU-Ync"]`, because its tower files are a copy of AU-ASM.

### A2. `main()` begins (line 157)

1. It creates `hydrology/results/` if it doesn't exist.
2. It makes five empty lists that will collect rows: `site_rows`, `clim_rows`, `eb_rows`, `annual_rows` and `monthly_rows`.
3. **The site loop** goes through every `.nc` file in `data/harmonised/` in alphabetical order (AU-ASM, AU-Cow, AU-Cpr, ...). It skips AU-Ync and prints the site name, so you see progress.
4. **`ds = xr.open_dataset(path)`** opens that site's file. `ds` holds, on one time axis:
   - tower variables (Qle, Qh, Rnet, SWdown, Tair_C, Psurf, Precip_mm and their `_qc` flags);
   - model variables, along a `source` dimension (obs, CABLE, JULES, ...);
   - attributes such as `timestep_s` (1800 or 3600), latitude and vegetation type.

Then, for that one site, it calls `site_weather` (Part B), `site_metadata` (Part C), and for each of obs, CABLE and JULES the four row-making functions (Part D). After the source loop it also calls `energy_balance_row` (Part E). When all sites are done, it writes the CSVs (Part F).

---

## PART B. `site_weather(ds)`: all the daily data for one site

This is where all the data filtering happens. It builds a dictionary `w` that everything later reads from.

### B1. `w["aet"], pet_paired = ets.analysis_days(ds)`

`analysis_days` (et_signatures.py line 168) runs three steps.

**Step 1. Daily ET for each source.** It calls `daily_et(ds, source)` for obs, CABLE and JULES. Inside `daily_et` (line 120), for one source:

1. `step_seconds` is 1800 (half-hourly) or 3600 (hourly, at AU-Tum and AU-Otw).
2. `qle` is that source's latent heat flux at every step, in W/m².
3. `lam = latent_heat(Tair_C)` is the energy needed to evaporate 1 kg of water at that step's temperature: (2.501 − 0.00237·T) × 10⁶ J/kg. It always uses the **tower** temperature, even for the models.
4. `et_per_step = qle / lam * step_seconds`. Watts are joules per second, so W/m² ÷ J/kg × s = kg/m² of water in that step, and 1 kg/m² = 1 mm.
5. `daytime = SWdown >= 5` gives True/False at each step. It always uses the **tower** SWdown.
6. `measured_daily(et_per_step, is_measured(ds, "Qle"), daytime)` turns the steps into days. **It uses the tower's Qle qc flags for every source**, so CABLE and JULES lose exactly the same days as the tower. See B1a.
7. `et_daily.where(et_daily >= 0)` sets any day whose total ET is negative to NaN (the paper excludes negative AET). This is done per source, so CABLE could lose a day the tower keeps.

**B1a. `is_measured(ds, "Qle")`** (line 28) is True where `Qle_qc == 0` and False otherwise. Note: the docstring says qc < 20, but the code says `== 0`. With `== 0`:
- any QA-flagged or corrected step (flags 1-19) counts as gap-filled;
- AU-DaP rain, which is flagged 10 everywhere, counts as entirely gap-filled.

**B1b. `measured_daily(per_step, measured, daytime, how="sum")`** (line 39) goes through these steps:

1. `gapfilled = ~measured`: True on gap-filled steps.
2. `day_gap = gapfilled & daytime`: True only on gap-filled steps **in daylight**. Night steps are always False here.
3. `run_id` numbers each unbroken run of True or False. `run_length` is the number of steps in the run each step belongs to.
4. `step_hours` is 0.5 or 1, worked out from the gap between the first two timestamps. `run_hours = run_length * step_hours`.
5. `in_long_day_gap = day_gap & (run_hours > 3)`: True on every step of a daylight gap-filled run **longer than** 3 h. That means more than 6 half-hours, or more than 3 hours at hourly sites.
   - A night step ends a run. A gap from 03:00 to 09:00 only counts its daylight part.
   - If SWdown is missing at a step, `NaN >= 5` is False, so that step counts as night.
6. `bad_day = in_long_day_gap.resample("1D").max()` is True for a calendar day if **any** of its steps is in a long daylight gap. A run that crosses midnight would mark both days bad, but daylight runs almost never cross midnight.
7. `daily = per_step.resample("1D").sum(min_count=1)`: the day's total from **all** steps, gap-filled ones included. `min_count=1` makes a day with no numbers at all NaN instead of 0.
8. `daily.where(~bad_day)` keeps good days and sets bad days to NaN.

So for each calendar day the result is the full-day total, or NaN. **Nothing is scaled up. Short gaps and night gaps are kept with their gap-filled values.** NaN days stay in the series as rows; they are not deleted.

**Step 2. Daily PET.** It calls `pet = daily_pet(ds)` (line 143), which does this:

1. Takes the tower `rnet` (W/m²), `tair` (°C) and `psurf_kpa` (Pa ÷ 1000).
2. Computes `daytime` from SWdown again.
3. Gets three daily values, each with **its own** qc flags and the same long-daylight-gap rule:
   - `rnet_day` = the day's sum of Rnet × step seconds, in J/m² per day ("sum");
   - `t_day` = the day's mean temperature ("mean");
   - `p_day` = the day's mean pressure ("mean").
4. If any of the three is NaN, the arithmetic below gives NaN, so PET is NaN that day.
5. Priestley-Taylor, with FAO-56 formulas:
   - `es = 0.6108·exp(17.27·T/(T+237.3))` is the saturation vapour pressure (kPa).
   - `delta = 4098·es/(T+237.3)²` is how fast es rises with temperature (kPa/°C).
   - `gamma = 0.000665·P` is the psychrometric constant (kPa/°C).
   - `PET = 1.26 · delta/(delta+gamma) · rnet_day / latent_heat(t_day)`, in mm/day.
6. No ground heat flux is subtracted. A day with negative net radiation gives negative PET; it is not clipped.

**Step 3. Same days for everything.**

1. `everything = pd.concat([ET obs, ET CABLE, ET JULES, PET], axis=1)` lines the four daily series up as four columns on one date index.
2. `keep = everything.notna().all(axis=1)` is True on a day only if **all four** have a value.
3. Each series has `.where(keep)` applied, so a day missing in any one of them becomes NaN in all of them.

**Result:**
- `w["aet"]` = `{"obs": ..., "CABLE": ..., "JULES": ...}`, three daily ET series with NaN in exactly the same places.
- `pet_paired` = PET on those same days.

### B2. `w["rain_daily"] = ets.daily_rain(ds)`

This is the daily rain sum from `Precip_mm`, with the rain qc flags and the same long-daylight-gap rule.

**Rain is NOT put on the ET days.** It keeps its own days, because rain is used to describe the climate and to find rain events, not to compare with ET day by day.

### B3. `w["pet_daily"] = ets.daily_pet(ds)`

This is PET again, **unpaired**: every day PET itself can be computed, whether or not ET exists. It is only used for the site's climate (`pet_mm_yr`, and from that aridity).

### B4. Monthly and annual PET from the paired days

```
w["pet_monthly"] = ets.totals(pet_paired, "MS")
w["pet_annual"]  = ets.totals(pet_paired, "YS")
```

**How `totals(daily, freq)` works** (line 183). This is the rule used for **every** monthly and annual number in the pipeline.

- `freq` is `"MS"` (month start, one bin per calendar month) or `"YS"` (year start, one bin per calendar year).
- `daily.resample(freq).size()` is the number of rows (days) in each bin, NaN days included. For a full month that is 28-31 days, and for a full year 365 or 366.
- `daily.resample(freq).mean()` is the mean of the days that have a value (NaN days are skipped).
- `total = mean × days in period`. The month's total is "an average day × number of days". A month with 25 valid days out of 30 is not counted as 25/30 of a month: its missing days are assumed to look like its average day. **It linearly extrapolates the mean to the full period.**
- `days_present = daily.notna().resample(freq).sum()` counts the days that have a value.
- `complete = days_present >= MIN_COMPLETE × days_in_period`. With 0.9:
  - a 30-day month needs at least 27 days with a value;
  - a 31-day month needs at least 27.9, so 28;
  - February needs at least 25.2, so 26;
  - a year needs at least 328.5, so 329 days (330 in a leap year).
- `total.where(complete).dropna()` **removes** incomplete months or years from the result entirely. The result is a series with only the complete periods.

One edge: "days in period" counts the rows that exist in the data. If a site's record started on 15 March, March would have 17 rows and be judged against 17, not 31, and its total would be 17 days' worth. The harmonised records start on 1 January (AU-Tum: 2002-01-01 to 2017-12-31), so in practice this doesn't happen.

### B5. Rain and PET climatology, and the dry and wet months

1. `rain_monthly = ets.totals(w["rain_daily"], "MS")` gives the complete rain months (the same 90% rule, on rain's own days).
2. `w["rain_clim"]` takes the mean of those months by calendar month: 12 numbers, mean January rain, mean February rain, and so on. A calendar month with no complete month in any year is simply missing.
3. `w["pet_clim"]` is the same for the paired PET months.
4. `w["dry_months"]` = the 3 calendar months with the least mean rain (`nsmallest(3)`), and `w["wet_months"]` = the 3 with the most.

---

## PART C. `site_metadata(ds, site, w)`: the site description columns

This is one dictionary that becomes the first columns of every `et_sites.csv` row for this site.

- **Copied from the file's attributes:** `site_name`, `veg`, `veg_long`, `lat`, `lon`, `elevation` and `canopy_height`.
- **`rain_mm_yr`:** the median of the complete rain years (`totals(rain_daily, "YS").median()`). It is NaN if the site has no complete rain year.
- **`pet_mm_yr`:** the median of the complete **unpaired** PET years.
- **`tair_mean`, `vpd_mean` and `lai_mean`:** plain means over the whole record, all steps, gap-filled included, no filtering.
- **`frac_qle_measured`:** the share of steps where the build script's `measured_Qle` mask is True (qc 0).
- **`dry_months` / `wet_months`:** the month numbers as text, for example "6 7 8".

---

## PART D. For each source (obs, then CABLE, then JULES)

`aet_daily = w["aet"][source]` is that source's paired daily ET. The same four functions run for each source.

### D1. `signature_row(meta, source, aet_daily, w)` produces one row of `et_sites.csv`

1. `aet_monthly = totals(aet_daily, "MS")` gives the complete months.
2. `aet_annual = totals(aet_daily, "YS")` gives the complete years.
3. `aet_clim` takes the mean of the complete months by calendar month, giving 12 or fewer values.
4. **The signatures** (paper numbering, Table 2, Gardiya Weligamage et al. 2025):
   - **`years`** = the number of complete years.
   - **1 `ET_median`** = the median of the complete annual totals (mm/yr).
   - **2 `CV_annual`** = std ÷ mean of the annual totals. With 1 year the std is NaN, so this is NaN. With 2 years it is computed, but the figures only show it for sites with 3 or more years.
   - **4 `TSP`** (`tsp`): groups the complete months by calendar month, takes the median of each, and returns the month number (1-12) with the highest median. It is NaN if there is no complete month at all.
   - **5 `CV_monthly`** = std ÷ mean of **all** complete monthly totals, every month of every year in one list.
   - **6 `WS`** = (mean monthly PET − mean monthly ET) ÷ mean monthly PET.
     - It uses `w["pet_monthly"]` (paired PET) and this source's complete months.
     - Because ET and PET sit on the same days, their complete months are almost always the same months.
     - Read it as 1 − ET/PET: 0 means ET keeps up with demand, and 1 means almost no ET.
   - **8 `R`** (`event_responsiveness`, line 222), with its count `n_events`:
     1. Event days are the days where daily rain is above 10 mm. They come from `rain_daily`, so they use rain's own days.
     2. For each event day, it takes `aet_daily[day : day + 10 days]`. Pandas label slicing includes both ends, so that is **11 days**, the rain day included.
     3. If none of those 11 days has ET, the event is skipped. Otherwise it records the day's rain and the **maximum** ET in the window.
     4. If fewer than 5 events are left, R is NaN.
     5. Otherwise it subtracts the mean from both lists and takes the Pearson correlation. Correlation is unchanged by subtracting means, so this is just the Pearson r between event rain and peak ET.
   - **`ET_driest3`** = the mean of `aet_clim` over the 3 driest calendar months. `ET_wettest3` is the same for the 3 wettest. `reindex` turns a missing month into NaN instead of an error, and `.mean()` skips NaN.
   - **`dry_wet_ratio`** = `ET_driest3` ÷ `ET_wettest3`.

### D2. `water_stress_rows(..., "YS")` adds rows to `et_annual.csv`

1. `aet = totals(aet_daily, "YS")` gives the complete ET years.
2. For each complete year, if that year is also a complete **paired PET** year, it writes `site, source, year, ET, PET, WS` with WS = (PET − ET) / PET. Otherwise it skips the year.

### D3. `water_stress_rows(..., "MS")` adds rows to `et_monthly.csv`

This is exactly the same, but per complete month. The `month` column is the date of the month's first day.

### D4. `climatology_rows(...)` adds rows to `et_monthly_clim.csv`

For calendar months 1-12 it writes the source's mean ET, the site's mean rain (`rain_clim`) and its mean PET (`pet_clim`). A month with no data is NaN. That gives 12 rows per site and source.

---

## PART E. `energy_balance_row(ds, site, days)`: one row of `et_eb.csv`

`days` = `w["aet"]["obs"].notna()`, which is True on the paired analysis days.

1. **Energy gap.** It takes the half-hours where the build script marked **both** Qle and Qh as measured (`measured_Qle` & `measured_Qh`, qc 0), with no daily rule.
   - `raw` = tower Qle + Qh at those steps.
   - `corrected` = Qle_cor + Qh_cor, the fluxes scaled up by PLUMBER2 so they close the energy balance.
   - `eb_gap_pct = 100 × (mean corrected ÷ mean raw − 1)`, meaning "the corrected fluxes have this many % more energy".
2. **Corrected ET.**
   - Qle_cor is turned into mm per step exactly as in `daily_et`: same λ, same qc flags, same daytime rule.
   - It is then restricted to the analysis days that are non-negative.
   - `ET_obs_cor` = the median of its complete years.

---

## PART F. Writing the CSVs (end of `main()`)

When every site is done, the five lists are written to `hydrology/results/`:

| file | rows | made by |
|---|---|---|
| `et_sites.csv` | one per site × source (22 × 3 = 66) | D1 + C |
| `et_monthly_clim.csv` | site × source × 12 months | D4 |
| `et_eb.csv` | one per site | E |
| `et_annual.csv` | one per complete year per site × source | D2 |
| `et_monthly.csv` | one per complete month per site × source | D3 |

The annual and monthly files always have their header row, even if empty. `et_sites_lat_changed.csv` is **not** written by this script. It is a hand-made copy with 4 sites' positions nudged for the maps.

---

## PART G. `make_report.py`: tables to figures

### G1. What happens at import (top level, runs before `main()`)

1. **Settings.**
   - `OUT` = `hydrology/results/` and `FIGS` = `hydrology/figures/`, which is created if missing.
   - `USE_NUDGED_POSITIONS = True`.
   - `GROUP` puts each site in one of 4 climate groups, and `GCOL` gives each group a colour.
   - `MCOL` sets the colours: tower = black, CABLE = blue, JULES = red.
   - `METRICS` lists the 6 signatures plotted in Fig 2.
2. **Load the tables.**
   - `d` = `et_sites.csv`. If nudged positions are on, `lat`/`lon` are replaced by the values in `et_sites_lat_changed.csv`, matched by site.
   - `clim` = `et_monthly_clim.csv`, `eb` = `et_eb.csv` (indexed by site), `annual` = `et_annual.csv`, and `monthly` = `et_monthly.csv` (month read as a date).
3. **Derived numbers.**
   - `d["aridity"] = pet_mm_yr / rain_mm_yr`. It uses unpaired PET and the median complete rain year, and is NaN if either is missing.
   - `d["group"]` is the climate group.
   - `o`, `c` and `j` are the obs, CABLE and JULES rows, each indexed by site.
   - `order` = the sites sorted by aridity, wettest first. That order is used in Figs 6-8 and the WS plots.
   - `rel_err(m) = (model ET_median − tower ET_median) / tower ET_median`.
   - `tsp_c` / `tsp_j` = the model peak month minus the tower peak month, the short way round the year (−6 ... +5).
   - `grp[model]` = the median error per climate group for 5 error types: relative annual ET, WS, R, dry/wet ratio and CV_monthly.

### G2. `main()` then draws, in this order

1. **`fig1_sites_climate`**
   - Left: site positions, with dot size growing with aridity and coloured by group.
   - Right: rain vs PET per site on a log rain axis, with a PET = rain line.
2. **`fig2_model_vs_tower`**
   - Six panels: ET_median, CV_annual, CV_monthly, WS, R and dry/wet ratio.
   - Each site is one dot per model, with tower on x, model on y and a 1:1 line.
   - CV_annual only uses sites with 3 or more complete years.
   - The Pearson r is computed over sites that have both values.
   - Sites whose error is more than a quarter of the axis range are labelled.
3. **`fig4_errors_vs_aridity`**
   - Four panels: relative ET error, WS error, R error and dry/wet ratio error, each against aridity.
   - The Spearman rho and p are shown if at least 3 sites have values.
4. **`fig5_group_errors`**: horizontal bars of the median error per climate group (from `grp`), CABLE vs JULES, one panel per error type.
5. **`fig6_seasonal_cycles`**
   - 22 small panels, wettest site first.
   - Each shows the mean seasonal ET cycle for tower, CABLE and JULES from `et_monthly_clim.csv`, with rain as pale blue bars on a hidden second axis.
6. **`fig7_site_detail`**
   - Left: the TSP offsets per site.
   - Right: a heatmap of the per-site errors, divided by a fixed scale for each column so that ±1 is full colour. The printed numbers are the unscaled errors.
7. **`fig8_energy_closure`**
   - Left: the model relative ET error against the tower energy gap, with Spearman rho.
   - Right: median annual ET per site for raw tower, corrected tower, CABLE and JULES.
8. **`map_figures`**: one PNG per signature (ET_median, CV_annual, TSP, CV_monthly, WS, R). Each has four cartopy maps:
   - CABLE − obs;
   - JULES − obs, on the same colour scale as CABLE;
   - CABLE − JULES;
   - a sign map: both too low, both too high, disagree, or both match.

   Other details:
   - TSP differences go the short way round the year.
   - A map is skipped if no site has both model values.
   - A site with a NaN error shows as an unfilled dot.
9. **`ws_annual_all`**: WS for each complete year (`et_annual.csv`), one panel per site, three lines.
10. **`ws_monthly_all`**
    - WS for each complete month.
    - Missing months are filled with NaN so the line breaks instead of joining across a gap.
    - The y axis is cut at −1 to +1.
11. **`ws_per_site`**: the same yearly and monthly WS plots, one PNG per site, in `figures/ws_timeseries/`.

Every figure is saved at 200 dpi and closed, and the script prints its path.

---

## PART H. The filters, all in one place

From half-hour to signature, a number has to pass each of these filters in order:

1. **Site:** AU-Ync is skipped.
2. **Step → day (`measured_daily`):** the day is lost if it has a gap-filled **daylight** run longer than 3 h in that variable's qc flags. "Measured" means qc == 0. For ET the tower Qle flags are used for every source.
3. **Negative day:** a day whose total ET is below 0 is lost, for that source.
4. **Pairing (`analysis_days`):** the day is lost for everything unless obs ET, CABLE ET, JULES ET and PET all exist.
   - PET itself needs Rnet, Tair and Psurf to each survive step 2.
   - Rain is exempt: it keeps its own days.
5. **Day → month or year (`totals`):** the period is kept only if at least 90% of its days survived. The total is the mean day × days in the period.
6. **Signature-level requirements:**
   - CV_annual needs 2 or more years (the figures use 3 or more).
   - R needs 5 or more rain events with ET in the window.
   - TSP needs at least one complete month.
   - WS per year or month needs both ET and PET complete in that period.

## PART I. Things worth knowing

- **`MIN_MEASURED_PER_DAY` does nothing now.** It is only referenced in the commented-out old versions, and in the table title in `gap_filling_checks.py`.
- **`is_measured` uses `qc == 0`, but its docstring says qc < 20.** With `== 0`:
  - AU-DaP's rain (flag 10) is all "gap-filled", so its rain days only survive where daylight runs are 3 h or shorter. In practice that means no rain at AU-DaP.
  - Corrected or QA-flagged steps count as gaps.
- **Rain uses the daylight-gap rule too.** Rain at night is never checked.
- **Climate-description numbers** (`tair_mean`, `vpd_mean`, `lai_mean`) use all steps with no filtering. `frac_qle_measured` and the energy gap use the build script's `measured_*` masks rather than `is_measured`.
- **`pet_mm_yr` (and so aridity) uses unpaired PET days.** WS uses paired PET, so the PET in WS and the PET in aridity can differ slightly.

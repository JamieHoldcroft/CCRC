# Project briefing: CABLE vs JULES against flux-tower observations (UNSW CCRC)

This is background, not a task list. It answers the three questions you asked (what CCRC is and what the project should achieve, what comes first, where the work lives), then gives the domain knowledge you need to reason well here. Anything marked **UNCONFIRMED** is my best understanding and should be checked against the real files, or with me.

---

## 1. What is CCRC, and what should this project achieve?

### The institution
- **CCRC** is the **Climate Change Research Centre** at UNSW Sydney. It is one of Australia's main climate-science groups and a hub for land-surface modelling and model evaluation.
- My placement is a **21st Century Weather** undergraduate research scholarship. 21st Century Weather is an ARC Centre of Excellence hosted across several Australian universities.
- **Duration:** 6 weeks full-time, from **28 Sep 2026 to 13 Nov 2026**, paused for about a week (roughly 22–30 Oct) while I'm overseas at a conference.

### People
- **Dr Anna Ukkola** (CCRC; ARC DECRA Fellow) is my primary supervisor.
  - She co-leads **PLUMBER2**, the major international land-model benchmarking experiment.
  - She built the **PLUMBER2 flux-tower dataset** and the **FluxnetLSM** R package that produced the files in this repo. Her address is in the files' `Package contact` attribute.
  - Her 2016 paper showed land models badly overestimate evaporative drought.
- **Dr Mat Lipson** is my co-supervisor. He does land-surface modelling, leads Urban-PLUMBER, and is based in Canberra.

### About me (Jamie)
- Final-year UNSW student: BSc Computer Science + BSc Maths (Physics minor).
- Strong background in ML, information retrieval and applied maths, with first-authored ML papers.
- **New to land-surface science.** Explain the physics when it matters. You don't need to explain Python, statistics or ML basics.

### The project
The official title is **"How well can we simulate land surface processes across Australia?"** In practice: evaluate **CABLE** and **JULES**, two land surface models (LSMs), against flux-tower observations of the water, energy and carbon cycles, and write it up as a comparative report.

- **The scope is deliberately open.** I'm expected to design the analysis myself, with supervisor guidance.
- **Candidate focuses:** hydrology and evapotranspiration, carbon fluxes, or **extremes (drought, heat, high vapour-pressure deficit)**. Extremes are my leading interest. I prefer natural ecosystems over the urban angle.

### Why it matters
- Australia's climate model, **ACCESS**, uses **CABLE** as its land surface.
- The UK Met Office's Unified Model family, which the **Bureau of Meteorology (BoM)** uses, runs **JULES**. BoM has also been building JULES into hydrological forecasting.
- So Australia relies on both. The project feeds a wider collaboration across several universities, BoM and **ACCESS-NRI**, the national facility supporting ACCESS.
- A planned **Australian Land Model Intercomparison Project (AusLMIP)** was proposed in August 2026. It starts with CABLE and JULES-ES, motivated by the fact that most LSMs were developed for Northern Hemisphere climates.

### The gap in the literature
- CABLE and JULES appear together in big **multi-model benchmarks**:
  - PLUMBER1 (Best et al. 2015)
  - Ukkola et al. 2016
  - PLUMBER2 (Abramowitz et al. 2024)
  - Cranko Page et al. 2026
- But there are almost no dedicated **head-to-head, mechanistic** comparisons for Australia. The only direct pairwise study I know of is Dumedah & Walker (2014), which compared soil-moisture data assimilation at Yanco, NSW. It didn't look at surface fluxes or extremes.
- Anna says this kind of comparison is scarce. I'm the first in the group doing it this way.

### Ambitions
- **Minimum:** a solid, well-evidenced report.
- **Stretch:** a publishable paper.
  - "Model A has lower RMSE than model B" is not a paper; much of that can be read off PLUMBER2's tables.
  - What would be publishable is **where and why** the models diverge, especially under stress, attributed to specific processes. For example: the stomatal response to VPD, soil-water stress, soil depth, or LAI. Ideally this would be confirmed with controlled model runs.

---

## 2. Where things stand now

- **Day 2 of the placement.** There is no pipeline code yet.
- **The repo currently holds:**
  - **Two dataset files.** One is certainly the PLUMBER2 **forcing ("Met") file** for Cumberland Plain: `AU-Cum_2013-2018_OzFlux_Met.nc`. Inspect the second to confirm what it is; it may be the matching Flux file.
  - **My exploratory notebook,** from playing with the data in Google Colab.
- **No CABLE or JULES output yet.** I'm waiting to hear from my supervisors what model output I'll get and in what form.
- **UNCONFIRMED — open questions I'm taking to Anna and Mat:**
  - Will I analyse existing **PLUMBER2 model outputs** (hosted on modelevaluation.org, ~150+ sites), run CABLE/JULES myself, or both?
  - **Site-level offline runs** (PLUMBER2-style) or **gridded high-resolution Australian runs**? The project title hints at gridded runs; those would compare a grid cell with a tower.
  - Exact model versions and configurations: the JULES configuration (e.g. GL9), CABLE's stomatal scheme, prescribed vs prognostic LAI, soil layers.
  - Should I evaluate against raw or energy-balance-corrected fluxes?
  - Access to **NCI Gadi** (Australia's national supercomputer, where CCRC data probably live, under `/g/data/...`) and to modelevaluation.org.
  - The preferred focus and audience of the report (ACCESS-NRI? BoM? a paper?).

---

## 3. Where the work lives

- **This repo** is the home for the work. It doubles as my learning record, so clear docstrings and explanations are welcome.
- **My past exploration** was in Google Colab, with data on Google Drive. That's why notebook paths look like `/content/drive/MyDrive/...`. Model outputs and bigger data will likely end up on NCI Gadi (often via its Jupyter service, ARE).
- **The stack:** Python with xarray, netCDF4, numpy, pandas, scipy and matplotlib; dask if scaling to many sites or gridded data. Command-line netCDF tools (`ncdump -h`, NCO, CDO) are common in this community.
- **Existing community tools, worth knowing so we don't reinvent them:**
  - modelevaluation.org, the PLUMBER2 evaluation platform;
  - ILAMB, the International Land Model Benchmarking package;
  - the FluxnetLSM R package, which produced these files;
  - Cranko Page et al.'s released analysis code (Zenodo, doi 10.5281/zenodo.17962998).

---

## 4. Physics you need to understand the data

### Surface energy balance
**Rnet = Qh + Qle + Qg** (+ small storage terms), all in W m⁻².

| Term | Meaning |
|---|---|
| **Rnet** | Net radiation: incoming shortwave + longwave, minus reflected and emitted. |
| **Qh** | Sensible heat: heat carried away by turbulence, driven by the surface-to-air temperature difference. |
| **Qle** | Latent heat: evapotranspiration expressed as energy. |
| **Qg** | Ground heat flux, into the soil. |

**The surface temperature is the balancing knob.** If water is limited, Qle can't take its share, the surface heats up and more energy goes into Qh. So most model errors show up as a **wrong split between Qh and Qle**.

Useful derived quantities:
- **Evaporative fraction:** EF = Qle / (Qle + Qh). It's bounded between 0 and 1 and fairly stable at midday.
- **Bowen ratio:** Qh / Qle.
- **Converting latent heat to water:** ET (mm per step) = Qle / λ × Δt, with λ ≈ 2.45 × 10⁶ J kg⁻¹. So 300 W m⁻² for 30 minutes is about 0.22 mm.

### Water balance
**P = ET + runoff + drainage + ΔS.**
- ET links the water and energy balances.
- About **30 % of PLUMBER2 sites have ET > P** over their whole record (groundwater, lateral inflow, irrigation). Most LSMs structurally can't reproduce that.

### Penman–Monteith: the attribution lens
ET is controlled by four things:
1. available energy (Rn − G);
2. vapour pressure deficit (**VPD**, the air's drying power);
3. aerodynamic conductance **g_a** (wind, roughness, canopy height);
4. surface or stomatal conductance **g_s** (biology plus soil-water stress).

Any CABLE–JULES difference in ET should be traceable to one or more of these. Rearranging Penman–Monteith gives an "observed" surface conductance you can compare with each model's.

Forests are tall and rough, so they are strongly **coupled** to the atmosphere: stomata and VPD dominate their ET. Short grass is more **decoupled**: energy dominates.

### Humidity
- Saturation vapour pressure rises exponentially with temperature (Clausius–Clapeyron). The Tetens approximation: e_s(T) = 0.6108·exp(17.27T/(T+237.3)) kPa, with T in °C.
- **Specific humidity** (`Qair`) is absolute and fairly flat through the day.
- **Relative humidity** (RH) swings hugely each day, just because temperature swings.
- **VPD** = e_s − e peaks mid-afternoon. Hot afternoons produce very high VPD even with ordinary moisture.

### Stomata and photosynthesis
Both models are "third-generation" LSMs, coupling carbon uptake and water loss through stomata:
- **Photosynthesis** is the Farquhar (C3) or Collatz (C4) model: the minimum of a Rubisco limit (set by Vcmax) and a light limit.
- **Stomatal conductance** rises with photosynthesis and falls as VPD rises. The formulations differ:
  - **Leuning:** CABLE's traditional default. Hyperbolic decline with VPD.
  - **Medlyn:** optimality-based, declining with 1/√VPD. Available in both models.
  - **Jacobs:** JULES's traditional default. The ratio of internal to external CO₂ falls linearly with VPD to zero at a critical deficit D_crit, which is a hard cutoff.
- **These diverge most at high VPD.** That's exactly where Cranko Page et al. found both models fail for Qh.

### Soil-moisture stress
- A factor β between 0 and 1 (called `fwsoil` in CABLE, `beta` in JULES).
- It is 1 above a critical moisture level, 0 at the wilting point, and a ramp in between, weighted by the root distribution across soil layers.
- **What β multiplies** (Vcmax, the stomatal slope, or net photosynthesis), the soil depth and the root profile together largely determine drought behaviour.
- Ukkola et al. 2016: LSMs simulate about **3× too many evaporative-drought days, at ~8× the magnitude**. Alternative soil-hydrology versions did better.

---

## 5. The two models

### How they run in this context: "offline" site simulations
Think of an LSM as a **hand-built recurrent model**:

h_t = f(h_{t−1}, x_t; θ), y_t = g(h_t, x_t; θ)

| Symbol | Meaning |
|---|---|
| x | The met forcing: exactly the variables in the Met file, one half-hour at a time. |
| h | Internal state: soil moisture and temperature per layer, snow, canopy water, carbon pools. |
| θ | Parameters, chosen from lookup tables by vegetation class (from `IGBP_veg_short`) and soil type. **Not fitted to the site.** PLUMBER2 runs were "out of the box". |
| y | Output fluxes (Qle, Qh, NEE, …), written every step. |

- **Spin-up:** the forcing years are looped until h stops drifting, then the analysis run is done.
- **It's a single deterministic run.** Not an ensemble, not a forecast, and no data assimilation.
- **The models never see the observed fluxes.** At every step they get the true observed weather, but their soil water is their own accumulated calculation. So errors build up over dry spells, which is why drought is such a good test.
- This framing is why PLUMBER2 benchmarks LSMs against LSTMs and random forests given the same inputs. **On average, no LSM beat a linear regression on shortwave radiation** in its independent scoring. For Qh, the physics models do especially badly.

### Inside one timestep, both models do roughly this
1. Partition radiation between leaves and soil, and between sunlit and shaded leaves.
2. Compute photosynthesis and stomatal conductance, given light, temperature, CO₂ and VPD.
3. Apply soil-water stress from the model's own soil moisture.
4. Solve the energy balance for the surface or leaf temperature, which gives Qh, Qle and Qg.
5. Update soil water and heat, and carbon pools.
6. Write the outputs.

### Structural differences
Defaults are shown; both models are highly configurable, so **confirm the actual configuration** of any runs.

| | CABLE | JULES |
|---|---|---|
| Origin / main use | CSIRO; ACCESS climate and Earth-system model; now community-developed with ACCESS-NRI support | UK Met Office + UKCEH; Unified Model weather and climate, UKESM; used by BoM |
| Design pressure | Climate and ecology; built with Australian ecosystems in mind | Weather prediction and climate; numerical robustness |
| Canopy radiation | Two big leaves (sunlit/shaded), direct vs diffuse light | Configurable; recent configurations use a multilayer canopy with sunlit/shaded fractions |
| Surface energy balance | **Two-source**: separate canopy and soil temperatures; canopy turbulence from Raupach theory; roughness depends on LAI | **One "skin" temperature per tile**, solved by linearisation |
| Stomata (default → option) | Leuning → Medlyn | Jacobs (D_crit) → Medlyn (from ~vn5.6) |
| Soil (default) | 6 layers, ~4.6 m deep; options for a groundwater aquifer and the SLI soil scheme | 4 layers, 3 m deep; tiles share one soil column by default |
| Carbon / vegetation extras | CASA-CNP (carbon, nitrogen **and phosphorus**), POP tree demography | TRIFFID dynamic vegetation, RothC soil carbon, JULES-CN |
| Vegetation types | IGBP-like classes | 5 or 13 plant types + non-vegetated tiles |
| Output naming | ALMA-style (`Qle`, `Qh`, `NEE`, `GPP`, `TVeg`, `ESoil`, `ECanop`, `SoilMoist`) | Own names (e.g. `latent_heat`, `ftl_gb`, `smcl`, `gpp_gb`; `_gb` = grid-box mean, `_tile` = per tile). **UNCONFIRMED:** exact names depend on version and output settings. |

**PLUMBER2 runs of these models:**
- `CABLE`: prescribed LAI, no carbon–nitrogen cycle.
- `CABLE-POP`: carbon–nitrogen cycle + demography.
- `JULES_GL9`: its own LAI.
- `JULES_GL9_LAI`: the site's prescribed LAI. This makes the pair a natural experiment on LAI alone.
- JULES submitted **no half-hourly NEE**.

**The central methodological issue:** a CABLE–JULES difference can come from four sources.
1. **Structure:** the equations.
2. **Parameters:** e.g. Vcmax, g1, root depth, soil hydraulics.
3. **Configuration:** soil depth, LAI source, spin-up, switches.
4. **Forcing handling:** how each model ingests the same met.

"Out of the box" runs mix all four together. Separating them is what would make this work more than a leaderboard. Model outputs that split evapotranspiration into transpiration, soil evaporation and canopy interception (`TVeg`, `ESoil`, `ECanop`) are very useful for explaining differences.

---

## 6. Flux-tower observations

### How they're measured
- **Eddy covariance.** A sonic anemometer measures 3-D wind, and a gas analyser (here a LI-7500) measures CO₂ and water vapour, both at 10–20 Hz on top of a tower.
- A flux is the **covariance** between vertical wind and the quantity of interest, averaged over **30 minutes**.
- The tower "sees" an upwind **footprint** of roughly 100 m to 1 km, which shifts with wind direction.

### Processing
- Raw data go through coordinate rotation, spectral corrections, the WPL density correction and QC.
- **OzFlux** (Australia's network, run with TERN) processes data through levels **L1 (raw) → L3 (QC and corrections) → L4/L5 (gap-filling) → L6 (partitioning)**, using PyFluxPro.
- PLUMBER2 takes L6 data and reformats it with FluxnetLSM into model-ready netCDF files. Only whole years are kept.

### Known data issues
- **Energy-balance non-closure.** Qh + Qle is typically only 70–90 % of (Rnet − Qg).
  - Models close the balance by construction.
  - "Corrected" fluxes (`*_cor`) scale Qh and Qle up while preserving the Bowen ratio.
  - PLUMBER2 argues the missing energy is more likely in Qh, so it treats **raw fluxes as primary**.
- **Night-time data are unreliable.** Calm, stable nights mean poor turbulence. This is handled with a u* (friction velocity) filter; PLUMBER2 also tested a wind > 2 m s⁻¹ filter.
- **Gaps.** Fluxes are often 20–60 % gap-filled, by an ANN called SOLO in OzFlux or by MDS elsewhere. **Only measured values (`_qc == 0`) are valid for evaluating models.** Evaluating against gap-filled values means evaluating against another statistical model.
- **GPP and ecosystem respiration are not measured.** They're derived by partitioning NEE with a model (night-time respiration vs temperature, or daytime light-response curves).
- **Sign convention:** NEE = respiration − GPP, so **negative means uptake**. Turbulent fluxes are positive upward.
- **Random measurement error** is roughly 10–20 %.
- **Forcing caveats:**
  - Longwave is often mostly *synthesised* from temperature and humidity (Abramowitz et al. 2012).
  - Rain gauges under-catch.
  - Pressure may be estimated from elevation.
  - CO₂ may be a global value rather than measured.

---

## 7. The data in this repo

### PLUMBER2 file pairs
Each site has two files:
- `<SITE>_<YEARS>_<Network>_Met.nc`: the **forcing**. It's what models are driven with, and it's continuous because it's gap-filled.
- `<SITE>_<YEARS>_<Network>_Flux.nc`: the **evaluation target**. Expected contents: `Qle`, `Qh`, `NEE`, `Rnet`, `Qg`, `GPP`, each with `_qc`, plus corrected versions of Qle and Qh (e.g. `Qle_cor`). **UNCONFIRMED** until inspected.

### Verified structure of `AU-Cum_2013-2018_OzFlux_Met.nc`
- **Dimensions** are `time: 105168`, `y: 1`, `x: 1`.
  - A single site is stored as a 1×1 grid so the same model code reads site and gridded files; squeeze it away.
  - 105,168 = 6 years × 17,520 half-hours + 48 extra for the 2016 leap day.
- **Time series** (all float32, dims `(time, y, x)`), each with a `_qc` companion: `SWdown`, `LWdown`, `Tair`, `Qair`, `RH`, `VPD`, `Wind`, `Psurf`, `Precip`, `CO2air`.
- **LAI** (no `_qc`): `LAI` and `LAI_alternative`. These are two satellite products (PLUMBER2 used MODIS and Copernicus); check the attributes for which is which.
- **Static:** `latitude`, `longitude`, `elevation`, `reference_height`, `canopy_height`.
- **Text:** `IGBP_veg_short` and `IGBP_veg_long` are byte strings (`|S200`). The vegetation class is probably EBF, evergreen broadleaf forest.
- **Units.** Tair is in **K**. `Precip` is a **rate** in kg m⁻² s⁻¹ (= mm s⁻¹), so multiply by Δt for an amount. VPD is **likely hPa** (UNCONFIRMED; check by recomputing from Tair and RH). Psurf is in Pa.
- **Naming differs between PLUMBER2 files.** Some use `Rainf`/`Snowf` instead of `Precip`, and `PSurf` instead of `Psurf`.
- **QC flags.** In the `_qc` variables, **0 = measured**; anything else means gap-filled. Check the attributes.
  - The global attribute `QC_flag_descriptions` lists **OzFlux's internal L1–L6 codes, not the `_qc` codes.** They're easy to confuse.

### Site facts from the global attributes
- **AU-Cum**, Cumberland Plain, Hawkesbury Valley, NSW (−33.615, 150.724).
- Dry sclerophyll woodland dominated by *Eucalyptus moluccana* and *E. fibrosa*, with spreading mistletoe.
- Canopy ~23 m, 30 m tower, elevation ~20 m, clay soil, mean annual precipitation ~800 mm.
- Run by Western Sydney University (HIE); the tower started on 19 Oct 2012.
- **The title "…for 2012" is legacy.** The data are 2013–2018.

### Time conventions (important and easy to get wrong)
- The `time_zone = Australia/Sydney` attribute names a zone that uses daylight saving. Flux data are normally kept in **local standard time (UTC+10) all year**.
- Stamps may mark the **end** of each half-hour. A first stamp of `00:30` and a last of `2019-01-01 00:00` would suggest this. **UNCONFIRMED.**
- **A robust check:**
  1. Compute top-of-atmosphere solar radiation from solar geometry (declination, equation of time, the offset of 150.72° E from the 150° E meridian).
  2. Find the time offset that maximises the correlation with the measured `SWdown`.
  3. Read the result:
     - −15 min → end-labelled stamps;
     - +15 min → start-labelled;
     - 0 → centred;
     - summer and winter differing by 60 min → daylight saving contamination;
     - about ±600 min → UTC.
- My tutorial notebook implements this check.
- **Model output has its own convention, and possibly a `noleap` calendar.** A half-hour misalignment looks like a large model error at sunrise and sunset.

### OzFlux sites likely in play later
Diverse Australian sites in PLUMBER2:
- **AU-How** (Howard Springs, tropical savanna, NT)
- **AU-Tum** (Tumbarumba, wet sclerophyll forest)
- **AU-ASM** (Alice Springs Mulga, arid)
- **AU-Cpr** (Calperum, mallee)
- AU-Whr (Whroo), AU-Wom (Wombat), and others

Site IDs are country code + three letters.

---

## 8. The evaluation landscape

**Benchmarking philosophy (Best et al. 2015; PLUMBER2).** Evaluating against observations alone ("the error is 30 W m⁻²") has no reference point. Instead, set *a priori* expectations with **out-of-sample empirical models** fed the same inputs:
- `1lin`: linear regression on shortwave.
- `2lin`: adds air temperature.
- `3km27`: k-means clustering plus regression.
- `6km729`, and `6km729lag`, which adds lagged rain and temperature as crude memory.
- Random forest and LSTM.

**PLUMBER2's metrics:**
- **Seven "independent" metrics:** mean bias, standard-deviation difference, correlation, normalised mean error, 5th- and 95th-percentile errors, and PDF overlap.
- They are combined into **iNMV**, normalised so that `1lin` = 1 and `LSTM` = 0.

**Cranko Page et al. 2026 — where models fail:**
- **"LSM Loss"** = timesteps where an LSM is worse than *all* of the best empirical models.
- Loss ratios are binned in 8-D meteorological space.
- Failures concentrate at **edge conditions**: unusual *combinations*, e.g. very high VPD for a given temperature. They are not simply extremes of one variable.
- Removing the worst ~12–31 % of timesteps greatly improves LSM scores. Several models, including CABLE and JULES_GL9, then beat the LSTM for Qle.
- CABLE and both JULES versions share a Qh failure at high VPD for a given temperature.

**The methods space, broadly:**
- **Standard statistics:** bias, RMSE, NME, correlation, KGE/NSE, Taylor diagrams.
- **Timescale decomposition:** diurnal, seasonal, anomaly, interannual.
- **Process diagnostics:**
  - EF vs antecedent rain or soil moisture;
  - surface conductance vs VPD;
  - GPP light response;
  - dry-down rates after rain;
  - ET partitioning;
  - water balance.
- **Conditional and event analysis:** heatwave and drought composites; temperature × VPD "fingerprints".
- **Rigour:**
  - **block** bootstrap, because half-hourly data are strongly autocorrelated;
  - observation uncertainty;
  - raw vs corrected sensitivity;
  - caution about multiple comparisons when slicing data many ways.
- **Attribution by controlled runs,** if model runs are possible: matched LAI, soil depth or stomatal scheme.

---

## 9. Key references

- **Pitman 2003**, *Int. J. Climatol.*: the evolution of land surface schemes. Energy and water coupling; the "generations" of LSMs.
- **Best et al. 2015**, *J. Hydrometeorol.*: PLUMBER. LSMs were beaten by linear regression for Qh.
- **Ukkola et al. 2016**, *Environ. Res. Lett.*: LSMs overestimate evaporative drought intensity, duration and magnitude.
- **Ukkola et al. 2022**, *ESSD*: the PLUMBER2 flux-tower dataset (doi 10.25914/5fdb0902607e1).
- **Abramowitz et al. 2024**, *Biogeosciences*: the PLUMBER2 experiment and first results.
- **Cranko Page et al. 2026**, *Biogeosciences* 23, 263: meteorological conditions where LSMs fail.
- **Dumedah & Walker 2014**, *Adv. Water Resour.*: JULES vs CABLE soil-moisture assimilation, Yanco.
- **CABLE:**
  - Kowalczyk et al. 2006 (description);
  - Haverd et al. 2018 (CABLE-POP);
  - De Kauwe et al. 2015 (Medlyn stomata in CABLE).
- **JULES:**
  - Best et al. 2011 and Clark et al. 2011 (description);
  - Harper et al. 2021 (drought stress);
  - Oliver et al. 2022 (Medlyn and thermal acclimation).
- **Community:** AusLMIP proposal (ACCESS Hive forum, Aug 2026); the ACCESS-NRI land surface model team.

---

## 10. How I'd like to work

- **Verify rather than assume.** Many details above are marked UNCONFIRMED for good reason. Where the data and this briefing disagree, trust the data and tell me.
- **Keep an `OPEN_QUESTIONS.md`** of domain or data questions for me to take to Anna and Mat.
- **Never modify the raw data files.**
- **Keep reusable logic out of notebooks,** and make things work per `site_code`, because the goal is eventually many sites and both models.
- **Explain domain reasoning briefly** in docstrings or commit messages. I'm learning the field as I go.
- **Plots:** labelled units, one y-axis per chart, and consistent colours for observations, CABLE and JULES throughout.
- **The direction is open.** I'd rather discuss options than have a plan imposed. I'll decide the first steps once I've seen what's in the repo and heard back from my supervisors.

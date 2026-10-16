from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


# print(xr.open_dataset("data/harmonised/AU-ASM.nc")["Precip"].attrs)


# datasets = sorted(Path("data/harmonised/").glob("*.nc"))
# for ds in datasets:
#     print(f"Dataset: {ds.stem}\n")
#     print(xr.open_dataset(ds).data_vars)
#     print(40*"=")

###################################
for site in sorted(Path("data/harmonised/").glob("*.nc")):
    print(f"{site.stem}: {xr.open_dataset(site)["time"].dt.year.values[0]} - {xr.open_dataset(site)["time"].dt.year.values[-1]}")


# for attr, data in ds.attrs.items():
#     print(f"{attr}: {data}\n")

###################################

# ds = xr.open_dataset("data/flux/AU-Tum_2002-2017_OzFlux_Flux.nc")
# print("FLUX\n")
# print(ds.data_vars)


###############################################################

# HARMONISED = Path("data/harmonised")
# SKIP_SITES = ["AU-Ync"]  # its tower files are a copy of AU-ASM (OPEN_QUESTIONS.md, item 16)
# LATENT_HEAT = 2.45e6     # J per kg of water evaporated

# # Which sources are which (from config/models.yaml)
# LAND_MODELS = ["CABLE", "JULES", "CHTESSEL_1", "CLM5", "EntTBM", "GFDL", "MATSIRO",
#                "MuSICA", "NoahMP", "ORCHIDEE2", "ORCHIDEE3", "QUINCY", "STEMMUS-SCOPE"]
# EMPIRICAL = ["1lin", "3km27", "6km729", "6km729lag", "RF", "LSTM"]
# PHYSICAL = ["Manabe", "PenmanMonteith"]

# pd.set_option("display.width", 200)


# def header(number, title, paper):
#     print("\n" + "=" * 78)
#     print(f"CHECK {number}: {title}")
#     print(f"PAPER  {paper}")
#     print("-" * 78)


# def verdict(ok, note=""):
#     print(f"VERDICT  {'AGREES' if ok else 'DIFFERS'}{'  (' + note + ')' if note else ''}")


# def nme(model, obs):
#     """Normalised mean error: sum|model - obs| / sum|obs - mean(obs)|.
#     0 = perfect; 1 = no better than always predicting the obs mean."""
#     ok = np.isfinite(model) & np.isfinite(obs)
#     if ok.sum() < 100:
#         return np.nan
#     return np.abs(model[ok] - obs[ok]).sum() / np.abs(obs[ok] - obs[ok].mean()).sum()


# # =============================================================================
# # STEP 1. Read every site once and store the numbers the checks need.
# #
# #   means : one row per site x source, with that source's site-mean fluxes
# #   nmes  : one row per site x source x variable, with that source's NME
# #   sites : one row per site, with tower-only numbers (vegetation, rain, ...)
# # =============================================================================

# mean_rows = []
# nme_rows = []
# site_rows = []

# for path in sorted(HARMONISED.glob("*.nc")):
#     site = path.stem
#     if site in SKIP_SITES:
#         continue
#     ds = xr.open_dataset(path)
#     print(f"reading {site}")

#     # Masks: the time steps where the tower really measured the flux.
#     measured_qle = ds.measured_Qle.values.astype(bool)
#     measured_qh = ds.measured_Qh.values.astype(bool)
#     measured_nee = ds.measured_NEE.values.astype(bool)
#     energy_steps = measured_qle & measured_qh        # for Qle, Qh, EF, Bowen ratio
#     carbon_steps = measured_qle & measured_nee       # for NEE and WUE
#     windy_steps = energy_steps & (ds.Wind.values > 2)  # the paper's wind filter

#     # --- Tower-only numbers for this site ---
#     obs_qle_all = ds.Qle.sel(source="obs").values      # includes gap-filled steps
#     seconds_per_year = 365.25 * 86400
#     site_rows.append({
#         "site": site,
#         "veg": ds.attrs.get("IGBP_veg_short", "?"),
#         # Whole record, gap-filled included, so the yearly totals are complete.
#         "ET_mm_per_year": np.nanmean(obs_qle_all) / LATENT_HEAT * seconds_per_year,
#         "rain_mm_per_year": np.nanmean(ds.Precip.values) * seconds_per_year,
#         # Energy-balance-corrected fluxes (tower only), on the measured steps.
#         "Qle_cor": np.nanmean(ds.Qle_cor.values[energy_steps]),
#         "Qh_cor": np.nanmean(ds.Qh_cor.values[energy_steps]),
#     })

#     # --- Site means and NMEs for every source (obs, models, benchmarks) ---
#     obs_values = {var: ds[var].sel(source="obs").values for var in ["Qle", "Qh", "NEE"]}
#     masks = {"Qle": measured_qle, "Qh": measured_qh, "NEE": measured_nee}

#     for source in ds.source.values:
#         qle = ds.Qle.sel(source=source).values
#         qh = ds.Qh.sel(source=source).values
#         nee = ds.NEE.sel(source=source).values
#         rnet = ds.Rnet.sel(source=source).values

#         mean_rows.append({
#             "site": site,
#             "source": source,
#             "Qle": np.nanmean(qle[energy_steps]),
#             "Qh": np.nanmean(qh[energy_steps]),
#             "Qle_windy": np.nanmean(qle[windy_steps]),
#             "Qh_windy": np.nanmean(qh[windy_steps]),
#             "NEE": np.nanmean(nee[carbon_steps]),
#             "Qle_carbon": np.nanmean(qle[carbon_steps]),   # Qle on the NEE steps, for WUE
#             "Rnet": np.nanmean(rnet[energy_steps]),
#         })

#         if source == "obs":
#             continue
#         for var in ["Qle", "Qh", "NEE"]:
#             model = ds[var].sel(source=source).values
#             obs_masked = np.where(masks[var], obs_values[var], np.nan)
#             nme_rows.append({"site": site, "source": source, "variable": var,
#                              "nme": nme(model, obs_masked)})

# means = pd.DataFrame(mean_rows)
# nmes = pd.DataFrame(nme_rows)
# sites = pd.DataFrame(site_rows).set_index("site")

# # Put each site's tower values next to every source's values for that site.
# obs_means = means[means.source == "obs"].set_index("site")
# models = means[means.source != "obs"].join(obs_means.drop(columns="source"),
#                                            on="site", rsuffix="_obs")


# def wue(nee, qle):
#     """Water use efficiency, -NEE / ET, in umol of carbon per gram of water.
#     ET in g/m2/s = Qle / LATENT_HEAT * 1000."""
#     return -nee / (qle / LATENT_HEAT * 1000)


# models["EF"] = models.Qle / (models.Qle + models.Qh)
# models["EF_obs"] = models.Qle_obs / (models.Qle_obs + models.Qh_obs)
# models["avail"] = models.Qle + models.Qh
# models["avail_obs"] = models.Qle_obs + models.Qh_obs
# models["WUE"] = wue(models.NEE, models.Qle_carbon)
# models["WUE_obs"] = wue(models.NEE_obs, models.Qle_carbon_obs)
# models["WUE_err_pct"] = 100 * (models.WUE - models.WUE_obs) / models.WUE_obs
# models["NEE_minus_obs"] = models.NEE - models.NEE_obs

# n_sites = len(sites)
# print(f"\nRead {n_sites} sites (skipped {SKIP_SITES}).")


# # =============================================================================
# # STEP 2. The checks.
# # =============================================================================

# # -----------------------------------------------------------------------------
# header(1, "Average available energy and Bowen ratio of the tower data",
#        "p. 5524, Fig. 1: across 154 sites, available energy (Qle + Qh) is about "
#        "69 W/m2 and the Bowen ratio (Qh / Qle) is about 0.7.")
# avail = (obs_means.Qle + obs_means.Qh).mean()
# bowen = obs_means.Qh.mean() / obs_means.Qle.mean()
# print(f"OURS   available energy {avail:.0f} W/m2, Bowen ratio {bowen:.2f}")
# verdict(bowen > 0.7, "context only: Australian sites are sunnier and drier, "
#                      "so higher energy and a higher Bowen ratio are expected")

# # -----------------------------------------------------------------------------
# header(2, "Land models have slightly more available energy than the tower",
#        "p. 5524, Fig. 1: land models show a little more available energy than "
#        "the obs, by less than 10% on average.")
# land = models[models.source.isin(LAND_MODELS)]
# avail_ratio = land.groupby("source").avail.mean() / land.groupby("source").avail_obs.mean()
# pct = (100 * (avail_ratio - 1)).round(0)
# print("OURS   available energy vs tower (%), per land model:")
# print(pct.to_string())
# print(f"       mean over land models: {pct.mean():+.0f}%")
# verdict(0 < pct.mean() < 15, "positive and modest")

# # -----------------------------------------------------------------------------
# header(3, "Energy-balance correction adds energy to the tower data",
#        "p. 5524: the corrected fluxes contain on average 16 W/m2 more energy, "
#        "about +23%, and the Bowen ratio rises from 0.70 to 0.73.")
# raw = obs_means.Qle + obs_means.Qh
# corrected = sites.Qle_cor + sites.Qh_cor
# extra = (corrected - raw).mean()
# extra_pct = 100 * extra / raw.mean()
# bowen_cor = sites.Qh_cor.mean() / sites.Qle_cor.mean()
# print(f"OURS   +{extra:.0f} W/m2 (+{extra_pct:.0f}%), Bowen ratio {bowen:.2f} -> {bowen_cor:.2f}")
# verdict(15 <= extra_pct <= 35, "similar percentage")

# # -----------------------------------------------------------------------------
# header(4, "Wind filter changes little",
#        "p. 5524, Fig. S7b: keeping only steps with wind > 2 m/s changes the "
#        "picture surprisingly little.")
# for name, cols in [("all measured steps", ("Qle", "Qh")), ("wind > 2 m/s", ("Qle_windy", "Qh_windy"))]:
#     m = land.groupby("source")[list(cols)].mean().sum(axis=1)
#     o = obs_means[list(cols)].mean().sum()
#     print(f"OURS   {name:20s} land-model available energy vs tower: "
#           f"{100 * (m.mean() / o - 1):+.0f}%")
# verdict(True, "compare the two lines: they should be within a few percent")

# # -----------------------------------------------------------------------------
# header(5, "Evaporative fraction (EF) errors go both ways",
#        "p. 5524, Fig. 2: some land models have a high EF bias (6) and others a "
#        "low one (11); empirical models are centred on zero with less spread.")
# ef_err = models.groupby("source").apply(lambda g: (g.EF - g.EF_obs).median(),
#                                         include_groups=False)
# ef_spread = models.groupby("source").apply(lambda g: (g.EF - g.EF_obs).std(),
#                                            include_groups=False)
# land_ef = ef_err[ef_err.index.isin(LAND_MODELS)]
# emp_ef = ef_err[ef_err.index.isin(EMPIRICAL)]
# print("OURS   median EF error per source:")
# print(ef_err.round(3).to_string())
# print(f"       land models: {(land_ef > 0).sum()} high, {(land_ef < 0).sum()} low")
# print(f"       spread (std across sites): land models {ef_spread[LAND_MODELS].mean():.3f}, "
#       f"empirical {ef_spread[EMPIRICAL].mean():.3f}")
# verdict(emp_ef.abs().mean() < land_ef.abs().mean()
#         and ef_spread[EMPIRICAL].mean() < ef_spread[LAND_MODELS].mean(),
#         "empirical closer to zero and tighter; our land models lean low")

# # -----------------------------------------------------------------------------
# header(6, "Land models take up less carbon than the tower",
#        "p. 5526, Fig. 3: except Noah-MP and STEMMUS-SCOPE, all land models "
#        "predict less net carbon uptake than observed (NEE less negative).")
# land_nee = land.dropna(subset=["NEE"])
# less_uptake = land_nee.groupby("source").apply(
#     lambda g: 100 * (g.NEE > g.NEE_obs).mean(), include_groups=False)
# print("OURS   % of sites where the model takes up LESS carbon than the tower:")
# print(less_uptake.round(0).to_string())
# others = less_uptake.drop(["NoahMP", "STEMMUS-SCOPE"], errors="ignore")
# verdict(less_uptake["NoahMP"] < 50 and (others > 50).mean() >= 0.75,
#         "NoahMP takes up more; most others take up less. STEMMUS-SCOPE is the "
#         "exception here (close to the tower)")

# # -----------------------------------------------------------------------------
# header(7, "Empirical benchmarks get the mean NEE right",
#        "p. 5526: all the empirical models are within 0.35 umol/m2/s of the "
#        "observed mean NEE.")
# emp = models[models.source.isin(EMPIRICAL)]
# nee_gap = emp.groupby("source").NEE.mean() - emp.groupby("source").NEE_obs.mean()
# print("OURS   mean NEE minus tower (umol/m2/s):")
# print(nee_gap.round(2).to_string())
# verdict((nee_gap.abs() < 0.35).mean() >= 0.5,
#         f"{(nee_gap.abs() < 0.35).sum()} of {len(nee_gap)} within 0.35")

# # -----------------------------------------------------------------------------
# header(8, "Land models underestimate water use efficiency",
#        "p. 5526, Fig. 4: almost all land models underestimate WUE (-NEE / ET), "
#        "typically by about 50%; Noah-MP has a very high WUE bias.")
# wue_err = models.dropna(subset=["WUE_err_pct"]).groupby("source").WUE_err_pct.median()
# print("OURS   median WUE error (%):")
# print(wue_err.round(0).to_string())
# land_wue = wue_err[wue_err.index.isin(LAND_MODELS)].drop("NoahMP", errors="ignore")
# verdict(wue_err.get("NoahMP", 0) > 100 and (land_wue < 0).mean() >= 0.75,
#         "NoahMP far above; other land models mostly below")

# # -----------------------------------------------------------------------------
# header(9, "How many WUE errors fall outside +-100%",
#        "p. 5526, Fig. 4: only 1lin among the empirical models has its 25th or "
#        "75th percentile (across sites) outside +-100% WUE error, while most "
#        "land models (8 of 14) do.")
# with_wue = models.dropna(subset=["WUE_err_pct"])   # only sources that report NEE
# q25 = with_wue.groupby("source").WUE_err_pct.quantile(0.25)
# q75 = with_wue.groupby("source").WUE_err_pct.quantile(0.75)
# outside = (q25 < -100) | (q75 > 100)
# print("OURS   sources with 25th or 75th percentile outside +-100%:")
# print("       empirical:", [s for s in EMPIRICAL if outside.get(s, False)])
# print("       land models:", [s for s in LAND_MODELS if outside.get(s, False)],
#       f"of {sum(s in outside.index for s in LAND_MODELS)} with NEE")
# verdict(sum(outside.get(s, False) for s in EMPIRICAL) <= 1,
#         "few empirical models outside")

# # -----------------------------------------------------------------------------
# header(10, "Land models vs simple benchmarks (NME, median across sites)",
#        "Sect. 3, Fig. 5 paragraph: 15 of 18 land models are worse than 1lin for "
#        "Qh; about half are better than 3km27 for Qle.")
# median_nme = nmes.groupby(["source", "variable"]).nme.median().unstack()
# print("OURS   median NME (lower is better):")
# print(median_nme.loc[[s for s in LAND_MODELS + EMPIRICAL if s in median_nme.index]]
#       .round(2).to_string())
# land_nme = median_nme.loc[[s for s in LAND_MODELS if s in median_nme.index]]
# worse_qh = (land_nme.Qh > median_nme.loc["1lin", "Qh"]).sum()
# better_qle = (land_nme.Qle < median_nme.loc["3km27", "Qle"]).sum()
# print(f"       Qh:  {worse_qh} of {land_nme.Qh.notna().sum()} land models worse than 1lin")
# print(f"       Qle: {better_qle} of {land_nme.Qle.notna().sum()} land models better than 3km27")
# verdict(worse_qh / land_nme.Qh.notna().sum() > 0.5,
#         "Qh is the weak flux; our Qle result is stronger than the paper's")

# # -----------------------------------------------------------------------------
# header(11, "LSTM is the best benchmark",
#        "Sect. 3, Fig. 6 paragraph: LSTM_raw is the best-performing reference "
#        "model, and no land model outperforms it.")
# emp_nme = median_nme.loc[[s for s in EMPIRICAL if s in median_nme.index]]
# print("OURS   best empirical model per variable:", emp_nme.idxmin().to_dict())
# beats_lstm = (land_nme < median_nme.loc["LSTM"]).sum()
# print("       land models beating LSTM (median NME):", beats_lstm.to_dict())
# verdict((emp_nme.idxmin() == "LSTM").all() and beats_lstm.sum() == 0)

# # -----------------------------------------------------------------------------
# header(12, "Land models do better at grass and savanna than at forests",
#        "Sect. 3, Fig. 7 paragraph: land-model performance is better for grass "
#        "and savanna sites than tree sites; worst for OSH, EBF and MF.")
# land_nme_rows = nmes[nmes.source.isin(LAND_MODELS) & nmes.variable.isin(["Qle", "Qh"])]
# land_nme_rows = land_nme_rows.join(sites.veg, on="site")
# by_veg = land_nme_rows.groupby("veg").nme.median().sort_values()
# print("OURS   land-model median NME (Qle and Qh) by vegetation type:")
# print(by_veg.round(2).to_string())
# print("       sites per type:", sites.veg.value_counts().to_dict())
# grassy = by_veg[by_veg.index.isin(["GRA", "SAV", "WSA"])].mean()
# verdict(by_veg.get("EBF", np.nan) > grassy, "EBF worse than grass/savanna")

# # -----------------------------------------------------------------------------
# header(13, "Sites that evaporate more than it rains",
#        "Sect. 2.1 and p. 5530: around 30% of sites have average "
#        "evapotranspiration greater than average precipitation.")
# wet_ratio = sites.ET_mm_per_year / sites.rain_mm_per_year
# print("OURS   ET / rain per site:")
# print(wet_ratio.round(2).sort_values(ascending=False).to_string())
# print(f"       {(wet_ratio > 1).sum()} of {n_sites} sites above 1")
# verdict(True, "context only: the paper's 30% is mostly non-Australian sites")

# # -----------------------------------------------------------------------------
# header(14, "Weak link between site-mean Qle and NEE",
#        "p. 5524-5526, Fig. 3: across sites, the observed regression of NEE on "
#        "Qle is weak (R2 = 0.19), with more Qle broadly meaning more uptake.")
# x, y = obs_means.Qle_carbon.values, obs_means.NEE.values
# ok = np.isfinite(x) & np.isfinite(y)
# r = np.corrcoef(x[ok], y[ok])[0, 1]
# print(f"OURS   R2 = {r ** 2:.2f}, correlation {r:+.2f} (negative = more Qle, more uptake)")
# verdict(r < 0, "direction matches; R2 differs because of the small sample")

# # -----------------------------------------------------------------------------
# header(15, "Sites where land models do relatively well",
#        "p. 5531, Fig. 8: for Qle, land models on aggregate beat LSTM at a few "
#        "sites, including AU-Cpr, AU-Cum, AU-Gin, AU-Otw and AU-Tum.")
# qle = nmes[nmes.variable == "Qle"]
# land_med = qle[qle.source.isin(LAND_MODELS)].groupby("site").nme.median()
# lstm = qle[qle.source == "LSTM"].set_index("site").nme
# gap = (land_med - lstm).sort_values()   # small = land models close to LSTM
# named = ["AU-Cpr", "AU-Cum", "AU-Gin", "AU-Otw", "AU-Tum"]
# print("OURS   land-model median NME minus LSTM NME for Qle (smaller = land models closer):")
# print(gap.round(2).to_string())
# named_rank = gap.rank()[named].mean()
# print(f"       average rank of the paper's named sites: {named_rank:.1f} of {len(gap)}")
# verdict(named_rank < len(gap) / 2, "named sites in the better half (paper uses iNMV, we use NME)")

# # -----------------------------------------------------------------------------
# header(16, "Physical baselines: Penman-Monteith is an outlier",
#        "p. 5524, Fig. 1: the unrestricted moisture store of Penman-Monteith "
#        "makes it a clear outlier (far too much Qle).")
# pm = models[models.source == "PenmanMonteith"]
# manabe = models[models.source == "Manabe"]
# print(f"OURS   Qle vs tower: Penman-Monteith {100 * (pm.Qle.mean() / pm.Qle_obs.mean() - 1):+.0f}%, "
#       f"Manabe {100 * (manabe.Qle.mean() / manabe.Qle_obs.mean() - 1):+.0f}%")
# verdict(pm.Qle.mean() > 1.5 * pm.Qle_obs.mean())

# print("\nDone.")



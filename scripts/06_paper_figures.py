"""Australian versions of the PLUMBER2 paper's figures, to compare side by side.

Paper: Abramowitz et al. (2024), Biogeosciences 21, 5517-5538.

Run from the repo root:
    python scripts/06_paper_figures.py

Writes to outputs/paper_figures/:
    fig1_qle_vs_qh.png        Fig. 1  site-mean Qle vs Qh per model
    fig2_ef_errors.png        Fig. 2  error in evaporative fraction, per model
    fig3_qle_vs_nee.png       Fig. 3  site-mean Qle vs NEE per model
    fig4_wue_errors.png       Fig. 4  error in water use efficiency, per model
    fig5_nme_vs_benchmarks.png Fig. 5-6 (simplified) NME of each land model vs 1lin, 3km27, LSTM
    fig7_nme_by_vegetation.png Fig. 7  land-model NME by vegetation type
    fig8_lstm_gain.png        Fig. 8  how much better LSTM is than the land models, per site,
                                      on axes of ET/rain vs dryness (PET/rain)

Same conventions as scripts/general.py: AU-Ync skipped (copy of AU-ASM),
averages over time steps the tower measured, NME instead of the paper's
dNMV/iNMV.
"""

import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

HARMONISED = Path("data/harmonised")
OUT = Path("outputs/paper_figures")
SKIP_SITES = ["AU-Ync"]
LATENT_HEAT = 2.45e6
SECONDS_PER_YEAR = 365.25 * 86400

LAND_MODELS = ["CABLE", "JULES", "CHTESSEL_1", "CLM5", "EntTBM", "GFDL", "MATSIRO",
               "MuSICA", "NoahMP", "ORCHIDEE2", "ORCHIDEE3", "QUINCY", "STEMMUS-SCOPE"]
EMPIRICAL = ["1lin", "3km27", "6km729", "6km729lag", "RF", "LSTM"]
PHYSICAL = ["Manabe", "PenmanMonteith"]
ALL_SOURCES = LAND_MODELS + EMPIRICAL + PHYSICAL

# One distinct colour per source; CABLE and JULES keep their usual colours.
palette = plt.get_cmap("tab20")
free = [i for i in range(20) if i not in (0, 1, 6, 7)]   # skip the blues and reds used by CABLE/JULES
COLOURS = {s: palette(free[i % len(free)]) for i, s in enumerate(LAND_MODELS[2:] + PHYSICAL)}
COLOURS.update({"CABLE": "#1f77b4", "JULES": "#d62728"})
COLOURS.update({s: str(0.25 + 0.1 * i) for i, s in enumerate(EMPIRICAL)})  # greys


def nme(model, obs):
    ok = np.isfinite(model) & np.isfinite(obs)
    if ok.sum() < 100:
        return np.nan
    return np.abs(model[ok] - obs[ok]).sum() / np.abs(obs[ok] - obs[ok].mean()).sum()


def wue(nee, qle):
    """-NEE / ET in umol of carbon per gram of water."""
    return -nee / (qle / LATENT_HEAT * 1000)


# -----------------------------------------------------------------------------
# Read every site once.
# -----------------------------------------------------------------------------
mean_rows, nme_rows, site_rows = [], [], []

for path in sorted(HARMONISED.glob("*.nc")):
    site = path.stem
    if site in SKIP_SITES:
        continue
    print(f"reading {site}")
    ds = xr.open_dataset(path)

    m_qle = ds.measured_Qle.values.astype(bool)
    m_qh = ds.measured_Qh.values.astype(bool)
    m_nee = ds.measured_NEE.values.astype(bool)
    energy_steps = m_qle & m_qh
    carbon_steps = m_qle & m_nee

    rain = np.nanmean(ds.Precip.values) * SECONDS_PER_YEAR            # mm/yr
    et = np.nanmean(ds.Qle.sel(source="obs").values) / LATENT_HEAT * SECONDS_PER_YEAR
    # PET: the Penman-Monteith baseline's ET (unlimited water), as in the paper's Fig. 8.
    pet = np.nan
    if "PenmanMonteith" in ds.source:
        pet = np.nanmean(ds.Qle.sel(source="PenmanMonteith").values) / LATENT_HEAT * SECONDS_PER_YEAR
    site_rows.append({"site": site, "veg": ds.attrs.get("IGBP_veg_short", "?"),
                      "ET_over_rain": et / rain, "PET_over_rain": pet / rain})

    obs = {v: ds[v].sel(source="obs").values for v in ["Qle", "Qh", "NEE"]}
    masks = {"Qle": m_qle, "Qh": m_qh, "NEE": m_nee}
    for source in ds.source.values:
        qle = ds.Qle.sel(source=source).values
        qh = ds.Qh.sel(source=source).values
        nee = ds.NEE.sel(source=source).values
        mean_rows.append({"site": site, "source": source,
                          "Qle": np.nanmean(qle[energy_steps]),
                          "Qh": np.nanmean(qh[energy_steps]),
                          "NEE": np.nanmean(nee[carbon_steps]),
                          "Qle_carbon": np.nanmean(qle[carbon_steps])})
        if source != "obs":
            for v in ["Qle", "Qh", "NEE"]:
                nme_rows.append({"site": site, "source": source, "variable": v,
                                 "nme": nme(ds[v].sel(source=source).values,
                                            np.where(masks[v], obs[v], np.nan))})

means = pd.DataFrame(mean_rows)
nmes = pd.DataFrame(nme_rows)
sites = pd.DataFrame(site_rows).set_index("site")

obs_means = means[means.source == "obs"].set_index("site")
models = means[means.source != "obs"].join(obs_means.drop(columns="source"),
                                           on="site", rsuffix="_obs")
models["EF_err"] = models.Qle / (models.Qle + models.Qh) - models.Qle_obs / (models.Qle_obs + models.Qh_obs)
models["WUE_err_pct"] = 100 * (wue(models.NEE, models.Qle_carbon)
                               - wue(models.NEE_obs, models.Qle_carbon_obs)) \
                        / wue(models.NEE_obs, models.Qle_carbon_obs)

OUT.mkdir(parents=True, exist_ok=True)
n = len(sites)


def save(fig, name):
    for ax in fig.axes:  # wrap long titles onto several lines
        ax.set_title("\n".join(textwrap.wrap(ax.get_title(), 110)), fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / name, dpi=130)
    plt.close(fig)
    print(f"wrote {OUT / name}")


def scatter_means(xvar, yvar, xlabel, ylabel, title, name, guides=True):
    """All-site mean of x vs y per source (big markers), site means (small grey dots)."""
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.scatter(obs_means[xvar], obs_means[yvar], s=12, color="0.8", label="sites (obs)", zorder=1)
    for s in ALL_SOURCES:
        g = models[models.source == s]
        if g.empty or g[[xvar, yvar]].isna().all().any():
            continue
        marker = "s" if s in EMPIRICAL else "^" if s in PHYSICAL else "o"
        ax.scatter(g[xvar].mean(), g[yvar].mean(), s=70, marker=marker,
                   color=COLOURS[s], edgecolor="k", lw=0.5, label=s, zorder=3)
    ox, oy = obs_means[xvar].mean(), obs_means[yvar].mean()
    ax.axvline(ox, color="k", lw=0.8, zorder=2)
    ax.axhline(oy, color="k", lw=0.8, zorder=2)
    ax.scatter(ox, oy, s=120, marker="s", color="k", label="obs", zorder=4)
    if guides:  # observed available energy (x + y = const) and Bowen ratio (x / y)
        xs = np.linspace(*ax.get_xlim(), 50)
        ax.plot(xs, ox + oy - xs, "k--", lw=0.6)
        ax.plot(xs, xs * oy / ox, "k:", lw=0.6)
        ax.set_xlim(xs[0], xs[-1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=7, ncol=2, loc="best")
    save(fig, name)


def boxplots(column, ylabel, title, name, hline=0, ylim=None, extra_lines=()):
    data, labels, colours = [], [], []
    for s in ALL_SOURCES:
        v = models.loc[models.source == s, column].dropna()
        if len(v):
            data.append(v.values)
            labels.append(s)
            colours.append(COLOURS[s])
    fig, ax = plt.subplots(figsize=(12, 5))
    box = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=True)
    for patch, c in zip(box["boxes"], colours):
        patch.set_facecolor(c)
    ax.axhline(hline, color="k", lw=0.8)
    for y in extra_lines:
        ax.axhline(y, color="m", lw=1.5, ls="--")
    ax.set_xticks(range(1, len(labels) + 1), labels, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel(ylabel)
    if ylim:
        ax.set_ylim(ylim)
    ax.set_title(title, fontsize=10)
    save(fig, name)


# -----------------------------------------------------------------------------
# Fig. 1: Qle vs Qh
# -----------------------------------------------------------------------------
scatter_means("Qh", "Qle", "Qh (W/m2)", "Qle (W/m2)",
              f"Like paper Fig. 1: mean Qle vs Qh over {n} AU sites (raw obs). "
              "Dashed = observed available energy, dotted = observed Bowen ratio.",
              "fig1_qle_vs_qh.png")

# Fig. 2: EF errors
boxplots("EF_err", "EF error (model - obs)",
         f"Like paper Fig. 2: error in site-mean evaporative fraction Qle/(Qle+Qh), {n} AU sites",
         "fig2_ef_errors.png")

# Fig. 3: Qle vs NEE
scatter_means("NEE", "Qle", "NEE (umol/m2/s, negative = uptake)", "Qle (W/m2)",
              f"Like paper Fig. 3: mean Qle vs NEE over {n} AU sites",
              "fig3_qle_vs_nee.png", guides=False)

# Fig. 4: WUE errors (%), with +-100% lines
boxplots("WUE_err_pct", "WUE error (% of obs)",
         f"Like paper Fig. 4: error in site-mean water use efficiency (-NEE/ET), {n} AU sites. "
         "Pink lines = +-100%",
         "fig4_wue_errors.png", ylim=(-250, 400), extra_lines=(-100, 100))

# -----------------------------------------------------------------------------
# Fig. 5 (simplified): each land model's median NME vs the three reference benchmarks
# -----------------------------------------------------------------------------
med = nmes.groupby(["source", "variable"]).nme.median().unstack()
fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=False)
for ax, var in zip(axes, ["Qle", "Qh", "NEE"]):
    lm = med.loc[[s for s in LAND_MODELS if s in med.index], var].dropna().sort_values()
    ax.barh(lm.index, lm.values, color=[COLOURS[s] for s in lm.index])
    for bench, ls in [("1lin", "--"), ("3km27", ":"), ("LSTM", "-")]:
        ax.axvline(med.loc[bench, var], color="k", ls=ls, lw=1, label=bench)
    ax.set_title(var)
    ax.set_xlabel("median NME across sites (lower is better)")
    ax.legend(fontsize=8)
fig.suptitle("Like paper Figs 5-6 (simplified, NME only): land models vs 1lin, 3km27 and LSTM", fontsize=10)
save(fig, "fig5_nme_vs_benchmarks.png")

# -----------------------------------------------------------------------------
# Fig. 7: land-model NME by vegetation type
# -----------------------------------------------------------------------------
lm_rows = nmes[nmes.source.isin(LAND_MODELS)].join(sites.veg, on="site")
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for ax, var in zip(axes, ["Qle", "Qh", "NEE"]):
    sub = lm_rows[lm_rows.variable == var].dropna(subset=["nme"])
    order = sub.groupby("veg").nme.median().sort_values().index
    ax.boxplot([sub[sub.veg == v].nme.values for v in order], showfliers=False)
    counts = sites.veg.value_counts()
    ax.set_xticks(range(1, len(order) + 1), [f"{v}\n(n={counts[v]})" for v in order], fontsize=8)
    lstm = nmes[(nmes.source == "LSTM") & (nmes.variable == var)].join(sites.veg, on="site")
    ax.plot(range(1, len(order) + 1), [lstm[lstm.veg == v].nme.median() for v in order],
            "k*", ms=10, label="LSTM median")
    ax.set_title(var)
    ax.set_ylabel("land-model NME across sites")
    ax.legend(fontsize=8)
fig.suptitle("Like paper Fig. 7: land-model NME by vegetation type (n = sites)", fontsize=10)
save(fig, "fig7_nme_by_vegetation.png")

# -----------------------------------------------------------------------------
# Fig. 8: how much better LSTM is than the median land model, per site
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
for ax, var in zip(axes, ["Qle", "Qh", "NEE"]):
    sub = nmes[nmes.variable == var]
    gain = sub[sub.source.isin(LAND_MODELS)].groupby("site").nme.median() \
        - sub[sub.source == "LSTM"].set_index("site").nme
    d = sites.join(gain.rename("gain")).dropna(subset=["gain"])
    sc = ax.scatter(d.PET_over_rain, d.ET_over_rain, c=d.gain, cmap="RdYlBu", vmin=-0.3, vmax=0.3,
                    s=60, edgecolor="k", lw=0.4)
    for s, row in d.iterrows():
        ax.annotate(s.replace("AU-", ""), (row.PET_over_rain, row.ET_over_rain), fontsize=6,
                    xytext=(3, 3), textcoords="offset points")
    ax.axhline(1, color="0.6", lw=0.6)
    ax.plot([0, 1], [0, 1], color="0.6", lw=0.6, ls="--")
    ax.set_xlabel("dryness index PET / rain (PET from Penman-Monteith)")
    ax.set_ylabel("ET / rain (obs)")
    ax.set_title(var)
    fig.colorbar(sc, ax=ax, label="median land-model NME - LSTM NME\n(red = land models better)")
fig.suptitle("Like paper Fig. 8: LSTM improvement over the median land model, per site", fontsize=10)
save(fig, "fig8_lstm_gain.png")

print("Done.")

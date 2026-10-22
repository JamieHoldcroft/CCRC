"""Generates the ET signature figures after having run step1_collect_results.py

Just reads the tables in hydrology/results/ and saves PNGs to hydrology/figures/
(results_fluxnet/ and figures_fluxnet/ when DATASET = "fluxnet" in et_signatures.py).
Works with CABLE only: any model missing from the tables is left out of every figure.

How the script is laid out:
    1. settings (folders, colours, site groups)
    2. load the tables
    3. numbers derived from the tables, used by several figures
    4. one function per figure (fig1 ... fig8)
    5. error maps
    6. water stress time series
    7. main(): calls every figure function in turn. Comment a line out there to
       skip a figure.
"""

import sys
from pathlib import Path

import matplotlib
import matplotlib.dates

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import pandas as pd
from scipy import stats


# ============================================================================
# 1. Settings
# ============================================================================
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import et_signatures as ets     # only for the DATASET switch (folder names)

OUT = HERE / f"results{ets.SUFFIX}"     # tables from step1_collect_results.py
FIGS = HERE / f"figures{ets.SUFFIX}"    # every figure is saved here
FIGS.mkdir(exist_ok=True)

# Map positions. True = take lat/lon for the maps from et_sites_lat_changed.csv,
# where 4 sites are nudged so their labels don't overlap. False = true positions.
# Either way the signatures come from et_sites.csv, which step1_collect_results.py
# rewrites every run. The nudged file only exists for the 22 PLUMBER2 sites.
USE_NUDGED_POSITIONS = ets.DATASET == "plumber2"

# Site groups. Each site's Koppen-Geiger climate (Rubel et al. 2017, 1986-2010,
# 100 arc-second map, looked up at the tower's lat/lon) decides its group:
#   Af/Am                -> Tropical rainforest (all three are rainforest towers)
#   Aw, and BSh north of the Tropic of Capricorn (23.4 S)  -> Tropical savanna
#   BWh, BSk, and BSh south of the Tropic                 -> Arid & semi-arid
#   Cfa, Cfb, Csa, Csb   -> Temperate & subtropical
# One exception, by vegetation: AU-Rob is tropical rainforest on the Atherton
# Tableland, which the map calls Cfa because of its 700 m altitude.
GROUP = {
    # Tropical rainforest
    "AU-Cow": "Tropical rainforest",       # Am  Cow Bay QLD, rainforest
    "AU-Ctr": "Tropical rainforest",       # Am  Cape Tribulation QLD, rainforest
    "AU-Rob": "Tropical rainforest",       # Cfa Robson Creek QLD, upland rainforest (exception above)
    # Tropical savanna
    "AU-How": "Tropical savanna",          # Aw  Howard Springs NT, eucalypt savanna
    "AU-Lit": "Tropical savanna",          # Aw  Litchfield NT, eucalypt savanna
    "AU-DaS": "Tropical savanna",          # Aw  Daly River Savanna NT
    "AU-DaP": "Tropical savanna",          # Aw  Daly River Pasture NT, cleared savanna
    "AU-Dry": "Tropical savanna",          # Aw  Dry River NT, savanna
    "AU-Ade": "Tropical savanna",          # Aw  Adelaide River NT, eucalypt savanna woodland
    "AU-RDF": "Tropical savanna",          # Aw  Red Dirt Melon Farm NT, cleared savanna
    "AU-Fog": "Tropical savanna",          # Aw  Fogg Dam NT, seasonally flooded wetland
    "AU-Stp": "Tropical savanna",          # BSh Sturt Plains NT (17 S), tropical grassland
    "AU-Fle": "Tropical savanna",          # BSh Fletcherview QLD (20 S), semi-arid tropical savanna
    # Arid & semi-arid
    "AU-ASM": "Arid & semi-arid",          # BWh Alice Springs Mulga NT
    "AU-TTE": "Arid & semi-arid",          # BWh Ti Tree East NT
    "AU-GWW": "Arid & semi-arid",          # BSh Great Western Woodlands WA (30 S)
    "AU-Emr": "Arid & semi-arid",          # BSh Emerald QLD (23.9 S), grassland
    "AU-Lon": "Arid & semi-arid",          # BSh Longreach QLD (23.5 S), Mitchell grass
    "AU-Cpr": "Arid & semi-arid",          # BSk Calperum SA, mallee
    "AU-Ync": "Arid & semi-arid",          # BSk Yanco NSW, grassland
    # Temperate & subtropical
    "AU-Cum": "Temperate & subtropical",   # Cfa Cumberland Plain NSW, dry sclerophyll forest
    "AU-Sam": "Temperate & subtropical",   # Cfa Samford QLD, pasture
    "AU-Rig": "Temperate & subtropical",   # Cfa Riggs Creek VIC, pasture
    "AU-Whr": "Temperate & subtropical",   # Cfa Whroo VIC, box-ironbark woodland
    "AU-MvA": "Temperate & subtropical",   # Cfa Myall Vale NSW, dryland pasture
    "AU-Wel": "Temperate & subtropical",   # Cfa Wellington NSW, pasture/woodland
    "AU-Ya1": "Temperate & subtropical",   # Cfa Yarramundi NSW, pasture
    "AU-Ya2": "Temperate & subtropical",   # Cfa Yarramundi NSW, irrigated pasture
    "AU-Tum": "Temperate & subtropical",   # Cfb Tumbarumba NSW, wet sclerophyll forest
    "AU-Wrr": "Temperate & subtropical",   # Cfb Warra TAS, wet eucalypt forest
    "AU-Otw": "Temperate & subtropical",   # Cfb Otway VIC, pasture
    "AU-Wac": "Temperate & subtropical",   # Cfb Wallaby Creek VIC, mountain ash forest
    "AU-Wom": "Temperate & subtropical",   # Cfb Wombat VIC, eucalypt forest
    "AU-APL": "Temperate & subtropical",   # Cfb Alpine Peatland VIC, 1670 m Sphagnum bog
    "AU-Sno": "Temperate & subtropical",   # Cfb Snow Gum NSW, 1616 m subalpine woodland
    "AU-Sil": "Temperate & subtropical",   # Cfb Silver Plains TAS, highland grassland
    "AU-Wgt": "Temperate & subtropical",   # Cfb Wedgetail TAS, highland grassland
    "AU-Gin": "Temperate & subtropical",   # Csa Gingin WA, banksia woodland
    "AU-Boy": "Temperate & subtropical",   # Csa Boyagin WA, wandoo woodland
    "AU-Rgf": "Temperate & subtropical",   # Csa Ridgefield WA, wheat
    "AU-Col": "Temperate & subtropical",   # Csb Collie WA, jarrah forest
}
GCOL = {"Tropical rainforest": "#1b7837", "Tropical savanna": "#d9a400",
        "Arid & semi-arid": "#b2182b", "Temperate & subtropical": "#2166ac"}
MCOL = {"obs": "k", "CABLE": "#1f77b4", "JULES": "#d62728"}

METRICS = [("ET_median", "1 Long-term median ET (mm/yr)"),
           ("CV_annual", "2 CV_annual"),
           ("CV_monthly", "5 CV_monthly"),
           ("WS", "6 WS (water stress)"),
           ("R", "8 R (event responsiveness)"),
           ("dry_wet_ratio", "Extra: dry/wet-season ET ratio")]

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10})
A4 = (8.27, 11.69)
A4L = (11.69, 8.27)


# ============================================================================
# 2. Load the tables (all written by step1_collect_results.py)
# ============================================================================
d = pd.read_csv(OUT / "et_sites.csv")                                 # one row per site x source
if USE_NUDGED_POSITIONS:
    nudged = pd.read_csv(OUT / "et_sites_lat_changed.csv").drop_duplicates("site").set_index("site")
    d["lat"] = d.site.map(nudged.lat)
    d["lon"] = d.site.map(nudged.lon)
clim = pd.read_csv(OUT / "et_monthly_clim.csv")                       # mean seasonal cycle
eb = pd.read_csv(OUT / "et_eb.csv").set_index("site")                 # tower energy gap
annual = pd.read_csv(OUT / "et_annual.csv")                           # WS each year
monthly = pd.read_csv(OUT / "et_monthly.csv", parse_dates=["month"])  # WS each month


# ============================================================================
# 3. Numbers derived from the tables
# ============================================================================
d["aridity"] = d.pet_mm_yr / d.rain_mm_yr      # PET / rain: depends on the PET estimate


missing = sorted(set(d.site) - set(GROUP))
if missing:
    raise SystemExit(f"Add these sites to GROUP at the top of this script: {', '.join(missing)}")
d["group"] = d.site.map(GROUP)

# One table per source, indexed by site
o = d[d.source == "obs"].set_index("site")
c = d[d.source == "CABLE"].set_index("site")
j = d[d.source == "JULES"].set_index("site")    # empty until JULES is in the tables
# The models in the tables, as (name, table) pairs. Every figure loops over these,
# so a missing JULES is simply left out.
MODELS = [(mod, m) for mod, m in [("CABLE", c), ("JULES", j)] if len(m)]
SOURCES = ["obs"] + [mod for mod, m in MODELS]
first = MODELS[0][1]                            # used where one model's value places a label
order = o.sort_values("aridity").index          # wettest to driest


def grid(n_panels, width, height_per_row, **kw):
    """4 columns and enough rows for n_panels (6 rows for the 22 PLUMBER2 sites)."""
    rows = int(np.ceil(n_panels / 4))
    return plt.subplots(rows, 4, figsize=(width, height_per_row * rows), squeeze=False, **kw)


def save_fig(fig, name):
    fig.savefig(FIGS / f"{name}.png", dpi=200)
    plt.close(fig)
    print(f"wrote {FIGS / name}.png")


def rel_err(m):
    """Relative annual ET error, (model - tower) / tower."""
    return (m.ET_median - o.ET_median) / o.ET_median


# Peak-month offsets (Fig 7) and group medians (Fig 5)
tsp = {mod: ((m.TSP - o.TSP + 6) % 12) - 6 for mod, m in MODELS}

grp = {}
for mod, m in MODELS:
    t = pd.DataFrame({"group": o.group, "Annual ET rel. error": rel_err(m),
                      "WS error": m.WS - o.WS, "R error": m.R - o.R,
                      "Dry/wet ratio error": m.dry_wet_ratio - o.dry_wet_ratio,
                      "CV_monthly error": m.CV_monthly - o.CV_monthly})
    grp[mod] = t.groupby("group").median()


# ============================================================================
# 4. One function per figure
# ============================================================================


def fig1_sites_climate():
    """Figure 1: site map and rain vs PET."""
    fig, axes = plt.subplots(1, 2, figsize=A4L)
    ax = axes[0]
    for s in o.index:
        aridity = o.loc[s, "aridity"]
        size = 40 if pd.isna(aridity) else 40 + 25 * aridity     # no aridity (no complete year of rain): smallest dot
        ax.scatter(o.loc[s, "lon"], o.loc[s, "lat"], s=size,
                   color=GCOL[GROUP[s]], edgecolor="k", lw=0.5, zorder=3)
        ax.annotate(s.replace("AU-", ""), (o.loc[s, "lon"], o.loc[s, "lat"]), fontsize=7,
                    xytext=(4, 3), textcoords="offset points")
    ax.set_xlim(112, 155)
    ax.set_ylim(-45, -10)
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")
    ax.set_title("Site locations (marker size = aridity)")
    ax.grid(alpha=0.3)
    for g, col in GCOL.items():
        ax.scatter([], [], color=col, edgecolor="k", label=g)
    ax.legend(fontsize=7, loc="lower left")
    ax = axes[1]
    for s in o.index:
        ax.scatter(o.loc[s, "rain_mm_yr"], o.loc[s, "pet_mm_yr"], color=GCOL[GROUP[s]],
                   edgecolor="k", lw=0.5, s=50)
        ax.annotate(s.replace("AU-", ""), (o.loc[s, "rain_mm_yr"], o.loc[s, "pet_mm_yr"]),
                    fontsize=7, xytext=(4, 3), textcoords="offset points")
    xx = np.geomspace(200, 4000, 50)
    ax.plot(xx, xx, "k--", lw=0.7, label="PET = rain")
    ax.set_xscale("log")
    ax.set_xlabel("median annual rain (mm/yr, log scale)")
    ax.set_ylabel("median annual PET (mm/yr)")
    ax.set_title("Climate: demand vs supply of water")
    ax.legend(fontsize=7)
    fig.suptitle("Figure 1. Where the sites are and what climates they span", fontsize=12, x=0.06, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save_fig(fig, "fig1_sites_climate")


def fig2_model_vs_tower():
    """Figure 2: model vs tower, one panel per signature."""
    fig, axes = plt.subplots(2, 3, figsize=A4L)
    for ax, (k, label) in zip(axes.ravel(), METRICS):
        sel = o.index if k != "CV_annual" else o.index[o.years >= 3]
        lo = min([o.loc[sel, k].min()] + [m.loc[sel, k].min() for mod, m in MODELS])
        hi = max([o.loc[sel, k].max()] + [m.loc[sel, k].max() for mod, m in MODELS])
        pad = 0.05 * (hi - lo)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "k--", lw=0.7)
        for mod, m in MODELS:
            ok = o.loc[sel, k].notna() & m.loc[sel, k].notna()      # sites with both values
            r = stats.pearsonr(o.loc[sel, k][ok], m.loc[sel, k][ok])[0] if ok.sum() >= 2 else np.nan
            ax.scatter(o.loc[sel, k], m.loc[sel, k], color=MCOL[mod], s=22, alpha=0.8,
                       label=f"{mod} (r = {r:.2f})")
        for s in sel:     # label sites where a model is far off
            if any(abs(m.loc[s, k] - o.loc[s, k]) > 0.25 * (hi - lo) for mod, m in MODELS):
                ax.annotate(s.replace("AU-", ""), (o.loc[s, k], max(m.loc[s, k] for mod, m in MODELS)),
                            fontsize=6.5, xytext=(2, 2), textcoords="offset points")
        ax.set_title(label)
        ax.set_xlabel("tower")
        ax.set_ylabel("model")
        ax.legend(fontsize=7)
    fig.suptitle("Figure 2. Model vs tower for each signature (one point per site; dashed = 1:1)",
                 fontsize=12, x=0.06, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save_fig(fig, "fig2_model_vs_tower")


def fig4_errors_vs_aridity():
    """Figure 4: errors against aridity."""
    fig, axes = plt.subplots(2, 2, figsize=A4L)
    panels = [("Annual ET relative error", lambda m: rel_err(m), "(model - tower) / tower"),
              ("WS error", lambda m: m.WS - o.WS, "model - tower"),
              ("R error (event responsiveness)", lambda m: m.R - o.R, "model - tower"),
              ("Dry/wet-season ratio error", lambda m: m.dry_wet_ratio - o.dry_wet_ratio, "model - tower")]
    for ax, (title, fn, ylab) in zip(axes.ravel(), panels):
        for mod, m in MODELS:
            e = fn(m)
            if (o.aridity.notna() & e.notna()).sum() >= 3:
                rho = stats.spearmanr(o.aridity, e, nan_policy="omit")
                label = f"{mod} (rho = {rho.statistic:.2f}, p = {rho.pvalue:.2f})"
            else:
                label = f"{mod} (too few sites for rho)"
            ax.scatter(o.aridity, e, color=MCOL[mod], s=22, alpha=0.8, label=label)
        for s in o.index:
            ax.annotate(s.replace("AU-", ""), (o.loc[s, "aridity"], fn(c)[s]), fontsize=6,
                        color="0.4", xytext=(2, 2), textcoords="offset points")
        ax.axhline(0, color="k", lw=0.7)
        ax.set_title(title)
        ax.set_xlabel("aridity  PET / rain")
        ax.set_ylabel(ylab)
        ax.legend(fontsize=7)
    fig.suptitle("Figure 4. Do the errors depend on how dry the site is?", fontsize=12, x=0.06, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save_fig(fig, "fig4_errors_vs_aridity")


def fig5_group_errors():
    """Figure 5: median error by site group."""
    fig, axes = plt.subplots(1, 5, figsize=A4L, sharey=True)
    groups = ["Tropical rainforest", "Tropical savanna", "Temperate & subtropical", "Arid & semi-arid"]
    for ax, col in zip(axes, grp["CABLE"].columns):
        y = np.arange(len(groups))
        for mod, shift in zip(grp, [0.2, -0.2]):     # CABLE above, JULES below
            ax.barh(y + shift, [grp[mod][col].get(g, np.nan) for g in groups], 0.38,
                    color=MCOL[mod], label=mod)
        ax.axvline(0, color="k", lw=0.7)
        ax.set_yticks(y, [g.replace(" & ", " &\n") for g in groups])
        ax.set_title(col, fontsize=9)
    axes[0].legend(fontsize=7)
    fig.suptitle("Figure 5. Median error by site group (model - tower)", fontsize=12, x=0.06, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save_fig(fig, "fig5_group_errors")


def fig6_seasonal_cycles():
    """Figure 6: mean seasonal cycle of ET at every site."""
    fig, axes = grid(len(order) + 1, A4[0], A4[1] / 6, sharex=True)   # +1 panel for the legend
    for ax, s in zip(axes.ravel(), order):
        cs = clim[clim.site == s]
        rain = cs[cs.source == "obs"].sort_values("month")
        ax2 = ax.twinx()
        ax2.bar(rain.month, rain.rain, color="#c6dbef", width=0.8, zorder=0)
        ax2.set_yticks([])
        ax.set_zorder(ax2.get_zorder() + 1)
        ax.patch.set_visible(False)
        for src in SOURCES:
            t = cs[cs.source == src].sort_values("month")
            ax.plot(t.month, t.ET, color=MCOL[src], lw=1.6 if src == "obs" else 1.1)
        ax.set_title(f"{s} ({o.loc[s, 'veg']}, PET/P {o.loc[s, 'aridity']:.1f})", fontsize=7.5)
        ax.tick_params(labelsize=6.5)
        ax.set_xticks([1, 4, 7, 10], ["J", "A", "J", "O"])
    for ax in axes.ravel()[len(order):]:
        ax.axis("off")
    ax = axes.ravel()[-1]
    for src in SOURCES:
        ax.plot([], [], color=MCOL[src], label="tower" if src == "obs" else src)
    ax.bar([], [], color="#c6dbef", label="rain (scaled)")
    ax.legend(loc="center", fontsize=8)
    fig.suptitle("Figure 6. Mean seasonal cycle of ET (mm/month), wettest site first", fontsize=11, x=0.06, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save_fig(fig, "fig6_seasonal_cycles")


def fig7_site_detail():
    """Figure 7: peak-month offsets and per-site error heatmap."""
    fig, axes = plt.subplots(1, 2, figsize=A4L, gridspec_kw={"width_ratios": [1, 1.6]})
    ax = axes[0]
    y = np.arange(len(order))
    shifts = [0.15, -0.15] if len(MODELS) == 2 else [0]      # side by side, or centred if one model
    for (mod, m), shift, marker in zip(MODELS, shifts, ["o", "^"]):
        ax.scatter(tsp[mod][order], y + shift, color=MCOL[mod], label=mod, s=30, marker=marker)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_yticks(y, [f"{s} (tower peak: {o.loc[s, 'TSP']:.0f})" for s in order], fontsize=7)
    ax.set_xticks(range(-4, 5))
    ax.invert_yaxis()                     # wettest site at the top, like the heatmap
    ax.set_xlabel("model peak month - tower peak month (months, circular)")
    ax.set_title("Timing of seasonal peak (TSP)")
    ax.legend(fontsize=7)
    ax.grid(axis="x", alpha=0.3)
    ax = axes[1]
    cols = [("Annual ET rel.", lambda m: rel_err(m), 0.5), ("WS", lambda m: m.WS - o.WS, 0.4),
            ("CV_monthly", lambda m: m.CV_monthly - o.CV_monthly, 0.4),
            ("R", lambda m: m.R - o.R, 0.4), ("Dry/wet", lambda m: m.dry_wet_ratio - o.dry_wet_ratio, 0.6)]
    mat, labels = [], []
    for name, fn, sc in cols:
        for mod, m in MODELS:
            mat.append((fn(m)[order] / sc).values)
            labels.append(f"{name}\n{mod[0]}")      # C or J
    mat = np.array(mat).T
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(labels)), labels, fontsize=6.5)
    ax.set_yticks(range(len(order)), order, fontsize=7)
    for i in range(mat.shape[0]):
        for k in range(mat.shape[1]):
            ax.text(k, i, f"{mat[i, k] * cols[k // len(MODELS)][2]:.2f}", ha="center", va="center", fontsize=5.5)
    ax.set_title("Errors per site (model - tower); "
                 + ", ".join(f"{mod[0]} = {mod}" for mod, m in MODELS))
    fig.colorbar(im, ax=ax, shrink=0.6, label="error / scale (red = model too high)")
    fig.suptitle("Figure 7. Site-by-site detail, wettest site at the top", fontsize=12, x=0.06, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save_fig(fig, "fig7_site_detail")


def fig8_energy_closure():
    """Figure 8: model errors vs tower energy gap; raw vs corrected tower ET."""
    fig, axes = plt.subplots(1, 2, figsize=A4L)
    ax = axes[0]
    for mod, m in MODELS:
        rho = stats.spearmanr(eb.loc[o.index, "eb_gap_pct"], rel_err(m), nan_policy="omit")
        ax.scatter(eb.loc[o.index, "eb_gap_pct"], rel_err(m), color=MCOL[mod], s=25,
                   label=f"{mod} (rho = {rho.statistic:.2f}, p = {rho.pvalue:.3f})")
    for s in o.index:
        ax.annotate(s.replace("AU-", ""), (eb.loc[s, "eb_gap_pct"], rel_err(MODELS[-1][1])[s]), fontsize=6.5,
                    color="0.4", xytext=(2, 2), textcoords="offset points")
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xlabel("tower energy gap (% extra energy in corrected fluxes)")
    ax.set_ylabel("annual ET relative error vs raw tower")
    ax.set_title("Model 'errors' vs how much energy the tower misses")
    ax.legend(fontsize=7)
    ax = axes[1]
    x = np.arange(len(order))
    ax.plot(x, o.loc[order, "ET_median"], "ks-", ms=4, label="tower (raw)")
    ax.plot(x, eb.loc[order, "ET_obs_cor"], "s--", color="0.5", ms=4, label="tower (energy-balance corrected)")
    for (mod, m), style in zip(MODELS, ["o-", "^-"]):
        ax.plot(x, m.loc[order, "ET_median"], style, color=MCOL[mod], ms=4, label=mod)
    ax.set_xticks(x, [s.replace("AU-", "") for s in order], rotation=70, fontsize=7)
    ax.set_ylabel("median annual ET (mm/yr)")
    ax.set_title("Annual ET: raw vs corrected tower, wettest to driest")
    ax.legend(fontsize=7)
    fig.suptitle("Figure 8. Part of the 'model error' is tower error", fontsize=12, x=0.06, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save_fig(fig, "fig8_energy_closure")


# ============================================================================
# 5. Error maps: map_<metric>.png, four maps of Australia per signature
# ============================================================================
# Maps are drawn with cartopy (install once with:  conda install -c conda-forge cartopy).
# LONLAT says "these numbers are longitude and latitude in degrees".
# MAP_PROJ is how the round Earth is flattened onto the page: a Lambert conformal
# projection centred on Australia, the usual choice for Australian maps.
LONLAT = ccrs.PlateCarree()
MAP_PROJ = ccrs.LambertConformal(central_longitude=134, central_latitude=-27,
                                 standard_parallels=(-18, -36),
                                 cutoff=10)   # cartopy's default cuts the map off at 30 S
                                              # (it assumes the north); 10 = stop at 10 N

# The signatures to map: (column in et_sites.csv, title, units of the error)
MAP_METRICS = [("ET_median", "1 Long-term median ET", "mm/yr"),
               ("CV_annual", "2 CV_annual", "-"),
               ("TSP", "4 TSP (month of peak ET)", "months"),
               ("CV_monthly", "5 CV_monthly", "-"),
               ("WS", "6 WS (water stress)", "-"),
               ("R", "8 R (event responsiveness)", "-")]


def month_difference(a, b):
    """a - b for months (1-12), going the short way round the calendar.
    December (12) minus January (1) is -1, not +11. The answer is -6 ... +5."""
    return ((a - b + 6) % 12) - 6


def draw_australia(ax):
    """Land, sea and coastline from cartopy; site dots go on top."""
    ax.set_extent([112, 155, -45, -9], crs=LONLAT)     # west, east, south, north
    ax.add_feature(cfeature.LAND, facecolor="#eeeeee", zorder=0)
    ax.add_feature(cfeature.OCEAN, facecolor="#f7fbff", zorder=0)
    ax.coastlines(resolution="50m", color="0.4", linewidth=0.7)
    ax.add_feature(cfeature.STATES, edgecolor="0.75", linewidth=0.4)   # state borders
    grid = ax.gridlines(draw_labels=True, x_inline=False, y_inline=False,
                        color="0.85", linewidth=0.5)
    grid.top_labels = False
    grid.right_labels = False
    grid.xlabel_style = {"size": 7}
    grid.ylabel_style = {"size": 7}


def error_map(ax, lon, lat, err, vmax, title, units):
    """One dot per site, coloured by the error. Red = too high, blue = too low."""
    draw_australia(ax)
    dots = ax.scatter(lon, lat, c=err, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                      s=70, edgecolor="k", linewidth=0.6, zorder=2, transform=LONLAT)
    cbar = plt.colorbar(dots, ax=ax, shrink=0.8)
    cbar.set_label(units)
    ax.set_title(title)


def sign_category(err_c, err_j):
    """Name the agreement between the two models' errors at each site."""
    labels = []
    for ec, ej in zip(err_c, err_j):
        if np.isnan(ec) or np.isnan(ej):
            labels.append("no value")
        elif ec == 0 and ej == 0:
            labels.append("both match obs")
        elif ec < 0 and ej < 0:
            labels.append("both too low")
        elif ec > 0 and ej > 0:
            labels.append("both too high")
        else:
            labels.append("disagree on sign")
    return pd.Series(labels, index=err_c.index)


SIGN_COLOURS = {"both too low": "#2166ac", "both too high": "#b2182b",
                "disagree on sign": "0.6", "both match obs": "white"}


def sign_map(ax, lon, lat, category, title):
    draw_australia(ax)
    for name, colour in SIGN_COLOURS.items():
        here = category == name
        if here.sum() == 0:
            continue
        ax.scatter(lon[here], lat[here], color=colour, s=70, edgecolor="k", linewidth=0.6,
                   zorder=2, transform=LONLAT, label=f"{name} ({here.sum()})")
    ax.legend(loc="lower left", fontsize=7, framealpha=0.9)
    ax.set_title(title)


def single_model_map(col, title, units, sites, lon, lat):
    """map_<metric>.png when the tables hold one model only: its error at each site."""
    mod, m = MODELS[0]
    if col == "TSP":
        err = month_difference(m.loc[sites, col], o.loc[sites, col])
    else:
        err = m.loc[sites, col] - o.loc[sites, col]
    if err.notna().sum() == 0:
        print(f"skipped map_{col}.png: no site has a value for it")
        return
    fig, ax = plt.subplots(figsize=(7, 5.5), subplot_kw={"projection": MAP_PROJ})
    error_map(ax, lon, lat, err, np.nanmax(np.abs(err.values)), f"{mod} error ({mod} - obs)", units)
    fig.suptitle(f"{title}: model errors at {len(sites)} sites", fontsize=13, weight="bold")
    fig.tight_layout()
    fig.savefig(FIGS / f"map_{col}.png", dpi=200)
    plt.close(fig)
    print(f"wrote {FIGS / f'map_{col}.png'}")


def map_figures():
    """One PNG per signature in MAP_METRICS."""
    for col, title, units in MAP_METRICS:
        # Which sites to show. CV_annual needs at least 3 complete years,
        # as everywhere else in this report.
        sites = o.index if col != "CV_annual" else o.index[o.years >= 3]
        lon, lat = o.loc[sites, "lon"], o.loc[sites, "lat"]

        # Without JULES there is only one map: CABLE's error.
        if len(MODELS) == 1:
            single_model_map(col, title, units, sites, lon, lat)
            continue

        # Step 1: the errors (model minus obs), and CABLE minus JULES.
        # For TSP the values are months, so subtract the short way round the year.
        if col == "TSP":
            err_c = month_difference(c.loc[sites, col], o.loc[sites, col])
            err_j = month_difference(j.loc[sites, col], o.loc[sites, col])
            err_cj = month_difference(c.loc[sites, col], j.loc[sites, col])
        else:
            err_c = c.loc[sites, col] - o.loc[sites, col]
            err_j = j.loc[sites, col] - o.loc[sites, col]
            err_cj = err_c - err_j          # the same as CABLE - JULES; obs cancels out

        if (err_c.notna() & err_j.notna()).sum() == 0:
            print(f"skipped map_{col}.png: no site has a value for it")
            continue

        # Step 2: one colour scale for the two model maps, so the same colour means
        # the same error in both. The scale is symmetric so white is exactly zero.
        vmax = np.nanmax(np.abs(np.concatenate([err_c.values, err_j.values])))
        vmax_cj = np.nanmax(np.abs(err_cj.values))

        # Step 3: draw the four maps
        fig, axes = plt.subplots(2, 2, figsize=(11, 9.5),
                                 subplot_kw={"projection": MAP_PROJ})   # every panel is a map
        error_map(axes[0, 0], lon, lat, err_c, vmax, "CABLE error (CABLE - obs)", units)
        error_map(axes[0, 1], lon, lat, err_j, vmax, "JULES error (JULES - obs)", units)
        error_map(axes[1, 0], lon, lat, err_cj, vmax_cj,
                  "CABLE - JULES (red: CABLE higher)", units)
        sign_map(axes[1, 1], lon, lat, sign_category(err_c, err_j),
                 "Do CABLE and JULES err the same way?")
        fig.suptitle(f"{title}: model errors at {len(sites)} sites", fontsize=13, weight="bold")
        fig.tight_layout()
        fig.savefig(FIGS / f"map_{col}.png", dpi=200)
        plt.close(fig)
        print(f"wrote {FIGS / f'map_{col}.png'}")


# ============================================================================
# 6. Water stress time series
# ============================================================================

# Water stress (WS) year by year. WS = 1 - ET/PET for each complete year,
# the same formula as signature 6 but for one year at a time (from
# step1_collect_results.py, et_annual.csv). Sites go wettest to driest.
def ws_annual_all():
    fig, axes = grid(len(order), 11.69, 14 / 6, sharey=True)
    for ax, site in zip(axes.flat, order):
        for source in SOURCES:
            one = annual[(annual.site == site) & (annual.source == source)]
            ax.plot(one.year, one.WS, marker="o", markersize=3, color=MCOL[source], label=source)
        ax.axhline(0, color="0.7", linewidth=0.6)
        ax.set_title(f"{site} ({o.loc[site, 'group']})", fontsize=8)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True, nbins=4))
        ax.ticklabel_format(axis="x", useOffset=False)     # show 2016, not "0 + 2.016e3"
        ax.tick_params(labelsize=7)
    for ax in axes.flat[len(order):]:
        ax.axis("off")                     # hide the spare panels (22 sites, 24 panels)
    for ax in axes[:, 0]:
        ax.set_ylabel("WS = 1 - ET/PET")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("Water stress each year, wettest site to driest\n"
                 "(0 = ET keeps up with PET, 1 = almost no ET)", fontsize=12, weight="bold")
    fig.tight_layout()
    fig.savefig(FIGS / "ws_annual.png", dpi=200)
    plt.close(fig)
    print(f"wrote {FIGS / 'ws_annual.png'}")


# Water stress month by month: WS = 1 - ET/PET for each complete month
# (et_monthly.csv from step1_collect_results.py). 12x more points than the yearly
# plot, so you can see the seasons and individual droughts. In winter PET is
# small, so WS swings far below zero when ET is larger than PET; the y axis is
# cut at -1 and those months run off the bottom of the panel.
def ws_monthly_all():
    fig, axes = grid(len(order), 11.69, 14 / 6, sharey=True)
    for ax, site in zip(axes.flat, order):
        for source in SOURCES:
            one = monthly[(monthly.site == site) & (monthly.source == source)]
            ws = one.set_index("month").WS
            if ws.empty:
                continue                     # no complete month for this source
            # Fill in missing (incomplete) months with NaN so the line breaks there
            # instead of joining across the gap.
            all_months = pd.date_range(ws.index.min(), ws.index.max(), freq="MS")
            ws = ws.reindex(all_months)
            ax.plot(ws.index, ws.values, color=MCOL[source], linewidth=0.8, label=source)
        ax.axhline(0, color="0.7", linewidth=0.6)
        ax.set_ylim(-1, 1)
        ax.set_title(f"{site} ({o.loc[site, 'group']})", fontsize=8)
        ax.xaxis.set_major_locator(matplotlib.dates.YearLocator(base=2))
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y"))
        ax.tick_params(labelsize=7)
    for ax in axes.flat[len(order):]:
        ax.axis("off")
    for ax in axes[:, 0]:
        ax.set_ylabel("WS = 1 - ET/PET")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("Water stress each month, wettest site to driest\n"
                 "(0 = ET keeps up with PET, 1 = almost no ET; y axis cut at -1)",
                 fontsize=12, weight="bold")
    fig.tight_layout()
    fig.savefig(FIGS / "ws_monthly.png", dpi=200)
    plt.close(fig)
    print(f"wrote {FIGS / 'ws_monthly.png'}")


# The same two water stress plots, but one PNG per site, so each site can be
# looked at on its own at a readable size. Saved into figures/ws_timeseries/:
#   <SITE>_WS_annual.png    WS each year   (from et_annual.csv)
#   <SITE>_WS_monthly.png   WS each month  (from et_monthly.csv)
def ws_per_site():
    WS_DIR = FIGS / "ws_timeseries"
    WS_DIR.mkdir(exist_ok=True)          # make the subfolder if it isn't there yet

    for site in order:
        # Yearly
        fig, ax = plt.subplots(figsize=(8, 4))
        for source in SOURCES:
            one = annual[(annual.site == site) & (annual.source == source)]
            ax.plot(one.year, one.WS, marker="o", markersize=4, color=MCOL[source], label=source)
        ax.axhline(0, color="0.7", linewidth=0.6)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        ax.ticklabel_format(axis="x", useOffset=False)
        ax.set_ylabel("WS = 1 - ET/PET")
        ax.set_title(f"{site} ({o.loc[site, 'group']}): water stress each year")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(WS_DIR / f"{site}_WS_annual.png", dpi=200)
        plt.close(fig)

        # Monthly (same gap handling and y limits as the all-sites plot above)
        fig, ax = plt.subplots(figsize=(12, 4))
        for source in SOURCES:
            one = monthly[(monthly.site == site) & (monthly.source == source)]
            ws = one.set_index("month").WS
            if ws.empty:
                continue                     # no complete month for this source
            all_months = pd.date_range(ws.index.min(), ws.index.max(), freq="MS")
            ws = ws.reindex(all_months)
            ax.plot(ws.index, ws.values, color=MCOL[source], linewidth=0.9, label=source)
        ax.axhline(0, color="0.7", linewidth=0.6)
        ax.set_ylim(-1, 1)
        ax.xaxis.set_major_locator(matplotlib.dates.YearLocator())
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y"))
        ax.set_ylabel("WS = 1 - ET/PET")
        ax.set_title(f"{site} ({o.loc[site, 'group']}): water stress each month (y axis cut at -1)")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(WS_DIR / f"{site}_WS_monthly.png", dpi=200)
        plt.close(fig)

    print(f"wrote {2 * len(order)} per-site WS plots to {WS_DIR}")



# ============================================================================
# 7. Main: make every figure. Comment out a line to skip that figure.
# ============================================================================
def main():
    fig1_sites_climate()
    fig2_model_vs_tower()
    fig4_errors_vs_aridity()
    fig5_group_errors()
    fig6_seasonal_cycles()
    fig7_site_detail()
    fig8_energy_closure()
    map_figures()
    ws_annual_all()
    ws_monthly_all()
    ws_per_site()


if __name__ == "__main__":
    main()

import textwrap
from pathlib import Path

import matplotlib
import matplotlib.dates

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"          # tables from collect_results.py
FIGS = HERE / "figures"         # every figure page is also saved here
FIGS.mkdir(exist_ok=True)

###############################################    CHANGE THIS TO 'et_sites.csv' FOR TRUE LAT/LONG    ############################################
# d = pd.read_csv(OUT / "et_sites.csv")
d = pd.read_csv(OUT / "et_sites_lat_changed.csv")

clim = pd.read_csv(OUT / "et_monthly_clim.csv")
eb = pd.read_csv(OUT / "et_eb.csv").set_index("site")

d["aridity"] = d.pet_mm_yr / d.rain_mm_yr
GROUP = {
    "AU-Cow": "Tropical rainforest", "AU-Ctr": "Tropical rainforest", "AU-Rob": "Tropical rainforest",
    "AU-How": "Tropical savanna", "AU-Lit": "Tropical savanna", "AU-DaS": "Tropical savanna",
    "AU-Dry": "Tropical savanna", "AU-DaP": "Tropical savanna", "AU-Stp": "Tropical savanna",
    "AU-ASM": "Arid & semi-arid", "AU-TTE": "Arid & semi-arid", "AU-GWW": "Arid & semi-arid",
    "AU-Cpr": "Arid & semi-arid", "AU-Emr": "Arid & semi-arid",
    "AU-Tum": "Temperate & subtropical", "AU-Wrr": "Temperate & subtropical",
    "AU-Whr": "Temperate & subtropical", "AU-Rig": "Temperate & subtropical",
    "AU-Otw": "Temperate & subtropical", "AU-Cum": "Temperate & subtropical",
    "AU-Sam": "Temperate & subtropical", "AU-Gin": "Temperate & subtropical",
}
GCOL = {"Tropical rainforest": "#1b7837", "Tropical savanna": "#d9a400",
        "Arid & semi-arid": "#b2182b", "Temperate & subtropical": "#2166ac"}
MCOL = {"obs": "k", "CABLE": "#1f77b4", "JULES": "#d62728"}
d["group"] = d.site.map(GROUP)

o = d[d.source == "obs"].set_index("site")
c = d[d.source == "CABLE"].set_index("site")
j = d[d.source == "JULES"].set_index("site")
order = o.sort_values("aridity").index          # wettest to driest
SHORT = list(o.index[o.years < 3])              # sites with < 3 complete years

METRICS = [("ET_median", "1 Long-term median ET (mm/yr)"),
           ("CV_annual", "2 CV_annual"),
           ("CV_monthly", "5 CV_monthly"),
           ("WS", "6 WS (water stress)"),
           ("R", "8 R (event responsiveness)"),
           ("dry_wet_ratio", "Extra: dry/wet-season ET ratio")]

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10})
A4 = (8.27, 11.69)
A4L = (11.69, 8.27)


# --------------------------------------------------------------------------
# Text-page helper
# --------------------------------------------------------------------------
class Writer:
    def __init__(self, pdf):
        self.pdf = pdf
        self.fig = None
        self.y = 0
        self.page = 0

    def new_page(self):
        if self.fig is not None:
            self._finish()
        self.fig = plt.figure(figsize=A4)
        self.y = 0.95

    def _finish(self):
        self.page += 1
        self.fig.text(0.5, 0.025, f"{self.page}", ha="center", fontsize=8, color="0.4")
        self.pdf.savefig(self.fig)
        plt.close(self.fig)
        self.fig = None

    def _need(self, h):
        if self.fig is None or self.y - h < 0.05:
            self.new_page()

    def title(self, text, size=20):
        self._need(0.06)
        self.fig.text(0.08, self.y, text, fontsize=size, weight="bold", va="top")
        self.y -= 0.05 if size >= 16 else 0.035

    def heading(self, text):
        self._need(0.06)
        self.y -= 0.01
        self.fig.text(0.08, self.y, text, fontsize=13, weight="bold", va="top", color="#1f3b73")
        self.y -= 0.032

    def para(self, text, size=9.5, width=105, indent=0.08, bullet=None, gap=0.008):
        lines = textwrap.wrap(text, width=width - (4 if bullet else 0))
        h = len(lines) * 0.0157 * size / 9.5 + gap
        self._need(h)
        x = indent + (0.025 if bullet else 0)
        if bullet:
            self.fig.text(indent + 0.005, self.y, bullet, fontsize=size, va="top")
        for line in lines:
            self.fig.text(x, self.y, line, fontsize=size, va="top")
            self.y -= 0.0157 * size / 9.5
        self.y -= gap

    def bullets(self, items, **kw):
        for it in items:
            self.para(it, bullet="•", **kw)

    def close(self):
        if self.fig is not None:
            self._finish()


def save_fig(pdf, w, fig):
    if w.fig is not None:          # finish any open text page first, so pages stay in order
        w.close()
    w.page += 1
    fig.text(0.5, 0.012, f"{w.page}", ha="center", fontsize=8, color="0.4")
    pdf.savefig(fig)
    fig.savefig(FIGS / f"page{w.page:02d}.png", dpi=200)
    plt.close(fig)


def caption(fig, text, y=0.085, width=150):
    fig.text(0.06, y, "\n".join(textwrap.wrap(text, width)), fontsize=8.5, va="top", style="italic")


# --------------------------------------------------------------------------
# Numbers used in the text
# --------------------------------------------------------------------------
def rel_err(m):
    return (m.ET_median - o.ET_median) / o.ET_median


stats_rows = []
for mod, m in [("CABLE", c), ("JULES", j)]:
    for k, label in METRICS:
        sel = o.index if k != "CV_annual" else o.index[o.years >= 3]
        x, y = o.loc[sel, k], m.loc[sel, k]
        e = y - x
        stats_rows.append({"Model": mod, "Metric": label, "n": len(sel),
                           "r (across sites)": stats.pearsonr(x, y)[0],
                           "median error": e.median(), "mean error": e.mean(),
                           "RMSE": np.sqrt((e ** 2).mean())})
stats_tab = pd.DataFrame(stats_rows)

h2h_rows = []
for k, label in METRICS:
    sel = o.index if k != "CV_annual" else o.index[o.years >= 3]
    ec = (c.loc[sel, k] - o.loc[sel, k]).abs()
    ej = (j.loc[sel, k] - o.loc[sel, k]).abs()
    p = stats.wilcoxon(ej, ec).pvalue
    h2h_rows.append({"Metric": label, "n": len(sel), "JULES closer": int((ej < ec).sum()),
                     "CABLE closer": int((ec < ej).sum()),
                     "median |err| CABLE": ec.median(), "median |err| JULES": ej.median(),
                     "Wilcoxon p": p})
h2h = pd.DataFrame(h2h_rows)

tsp_c = ((c.TSP - o.TSP + 6) % 12) - 6
tsp_j = ((j.TSP - o.TSP + 6) % 12) - 6

grp = {}
for mod, m in [("CABLE", c), ("JULES", j)]:
    t = pd.DataFrame({"group": o.group, "Annual ET rel. error": rel_err(m),
                      "WS error": m.WS - o.WS, "R error": m.R - o.R,
                      "Dry/wet ratio error": m.dry_wet_ratio - o.dry_wet_ratio,
                      "CV_monthly error": m.CV_monthly - o.CV_monthly})
    grp[mod] = t.groupby("group").median()

rho_cj = stats.spearmanr(rel_err(c), rel_err(j))
rho_r_c = stats.spearmanr(o.aridity, c.R - o.R)
rho_r_j = stats.spearmanr(o.aridity, j.R - o.R)
rho_et_c = stats.spearmanr(o.aridity, rel_err(c))
rho_et_j = stats.spearmanr(o.aridity, rel_err(j))
rho_lai_c = stats.spearmanr(o.lai_mean, rel_err(c))
rho_lai_j = stats.spearmanr(o.lai_mean, rel_err(j))
rho_eb_c = stats.spearmanr(eb.loc[o.index, "eb_gap_pct"], rel_err(c))
rho_eb_j = stats.spearmanr(eb.loc[o.index, "eb_gap_pct"], rel_err(j))
relc_cor = (c.ET_median - eb.ET_obs_cor) / eb.ET_obs_cor
relj_cor = (j.ET_median - eb.ET_obs_cor) / eb.ET_obs_cor
rho_ws_c = stats.spearmanr(c.ET_median - o.ET_median, c.WS - o.WS)
rho_ws_j = stats.spearmanr(j.ET_median - o.ET_median, j.WS - o.WS)

print("numbers ready")

# --------------------------------------------------------------------------
# Build the PDF
# --------------------------------------------------------------------------
pdf = PdfPages(HERE / "et_signatures_report.pdf")
w = Writer(pdf)

# ---- Title page and key findings ----
w.new_page()
w.title("Evapotranspiration signatures:", 20)
w.title("CABLE and JULES vs flux towers at 22 Australian sites", 15)
w.para("Jamie (CCRC placement), October 2026. Signatures from Table 2 (p. 6) of Gardiya "
       "Weligamage et al. (2025, EGUsphere, doi:10.5194/egusphere-2025-3373), computed with "
       "et_signatures.py on the harmonised PLUMBER2 site files. AU-Ync is excluded "
       "(its tower files are a copy of AU-ASM).", size=9, gap=0.02)
w.heading("Key findings")
w.bullets([
    f"Both models get the big picture right. Across sites, annual ET ranges from about 260 to "
    f"1,390 mm/yr and both models track it (cross-site correlation r = "
    f"{stats.pearsonr(o.ET_median, c.ET_median)[0]:.2f} for CABLE, "
    f"{stats.pearsonr(o.ET_median, j.ET_median)[0]:.2f} for JULES). CABLE's typical annual ET is "
    f"{100 * rel_err(c).median():.0f}% of the tower's value off, JULES's "
    f"{100 * rel_err(j).median():+.0f}%.",
    f"Their errors are shared. CABLE's and JULES's annual-ET errors rise and fall together "
    f"across sites (Spearman rho = {rho_cj.statistic:.2f}). At 13 of 22 sites both models are "
    f"off in the same direction. This points to shared model structure, or to the forcing and "
    f"tower data, more than to one model's quirks.",
    "The wettest sites are missed most. At the two lowland tropical rainforests (AU-Cow, "
    "AU-Ctr) both models evaporate 26-35% less than the tower, and 12-17% less at upland "
    "AU-Rob. At AU-Ctr the tower measures 1,390 mm/yr against 948 (CABLE) and 984 (JULES).",
    "Tropical savannas collapse in the dry season. At AU-How, AU-DaS and AU-Dry the real "
    "savanna's ET in its three driest months is 41-53% of its wettest-month ET; CABLE's is "
    "16-26%. This dry-season collapse also exaggerates CABLE's seasonal swing (CV_monthly "
    "+0.21 in this group). JULES does better but has the same tendency.",
    f"Arid sites: totals right, dynamics wrong. Both models match arid-site annual ET (it is "
    f"capped by rainfall), but both respond far more strongly to individual rain events than "
    f"the tower does (e.g. R = 0.47/0.49 vs 0.10 at AU-Cpr). The excess responsiveness grows "
    f"with dryness (rho = {rho_r_c.statistic:.2f}, p = {rho_r_c.pvalue:.2f} for CABLE).",
    "A clear missing process at AU-Gin. Gingin (Banksia woodland on deep sand) is the only "
    "site that evaporates more in its dry season than its wet season (ratio 1.27, peak in "
    "November), presumably from shallow groundwater. Neither model can do this: ratios 0.52 "
    "and 0.61, with peaks 2-3 months early.",
    f"JULES is better in some respects, not all. JULES is closer on water stress (15 of 22 "
    f"sites, Wilcoxon p = {h2h.loc[3, 'Wilcoxon p']:.3f}) and on dry-season persistence "
    f"(17 of 22, p = {h2h.loc[5, 'Wilcoxon p']:.3f}). For annual totals it wins 14 of 22 "
    f"(p = {h2h.loc[0, 'Wilcoxon p']:.2f}, not significant), and there is no difference for "
    f"CV_annual, CV_monthly or R. JULES gets the peak month exactly right at "
    f"{int((tsp_j == 0).sum())} sites against CABLE's {int((tsp_c == 0).sum())}; CABLE "
    f"tends to peak early.",
    f"Which tower data you compare against matters. JULES's apparent over-estimates sit at "
    f"sites where the tower misses most energy (rho = {rho_eb_j.statistic:.2f}, p = "
    f"{rho_eb_j.pvalue:.3f}). Against energy-balance-corrected tower ET, both models "
    f"under-estimate (median {100 * relc_cor.median():.0f}% CABLE, "
    f"{100 * relj_cor.median():.0f}% JULES).",
], size=9.3)

# ---- Data and methods ----
w.new_page()
w.heading("1. Data and method")
w.para("Data. One harmonised file per site (data/harmonised/<SITE>.nc) holds the tower "
       "observations and the CABLE and JULES PLUMBER2 runs on the same half-hourly axis. "
       "Both models were run 'out of the box' with the tower's weather as forcing.")
w.para("ET. Latent heat flux is converted to evapotranspiration with a temperature-dependent "
       "latent heat of vaporisation, lambda = (2.501 - 0.00237 T) x 10^6 J/kg, summed to days, "
       "months and years. As in the paper, the tower ET is gap-filled and negative values are "
       "excluded (here: negative daily totals). A month or year counts only if 90% of its "
       "days are present.")
w.para("PET. The paper uses Morton wet-environment PET from SILO. We use a substitute: "
       "Priestley-Taylor PET (alpha = 1.26) from the tower's net radiation and air "
       "temperature, the same PET for the tower and both models. This affects metric 6 (WS) "
       "and the aridity index used to order sites.")
w.para("Metrics. From Table 2 we use 1 (long-term median annual ET), 2 (CV of annual ET), "
       "4 (TSP, month of peak ET), 5 (CV of monthly ET), 6 (WS = (PET - AET)/PET) and 8 (R, "
       "event responsiveness: correlation between rain on days with more than 10 mm and the "
       "maximum daily ET in the following 0-10 days). Metrics 3 and 7 were left out on "
       "purpose. One extra diagnostic was added for this report: the dry/wet-season ratio = "
       "mean monthly ET in the 3 calendar months with least rain divided by the 3 with most. "
       "It measures whether a site keeps evaporating through its dry season.")
w.para("Comparison. For each metric we compare model and tower across the 22 sites "
       "(correlation, median error, RMSE), count where each model is closer, and test the "
       "difference in absolute error with a paired Wilcoxon signed-rank test. With only 22 "
       "sites, p-values are indicative, not definitive. TSP is compared as a circular month "
       "offset (December is one month from January). CV_annual uses only the 17 sites with "
       f"at least 3 complete years; the short sites are {', '.join(SHORT)}.")
w.para("Site groups. To read patterns, sites are grouped into tropical rainforest (AU-Cow, "
       "AU-Ctr, AU-Rob), tropical savanna (AU-How, AU-Lit, AU-DaS, AU-Dry, AU-DaP, AU-Stp), "
       "arid & semi-arid (AU-ASM, AU-TTE, AU-GWW, AU-Cpr, AU-Emr) and temperate & subtropical "
       "(AU-Tum, AU-Wrr, AU-Whr, AU-Rig, AU-Otw, AU-Cum, AU-Sam, AU-Gin). The aridity index is "
       "PET / rain: below 1 is wet (energy-limited), above 2 is dry (water-limited).")

# ---- Site overview table ----
fig, ax = plt.subplots(figsize=A4L)
ax.axis("off")
cols = ["Site", "Name", "Veg", "Group", "Lat", "Rain\nmm/yr", "PET\nmm/yr", "Aridity\nPET/P",
        "LAI", "Years", "Rain\nevents", "Energy\ngap %", "Tower ET\nmm/yr"]
cell = []
for s in order:
    r = o.loc[s]
    short_group = {"Tropical rainforest": "Rainforest", "Tropical savanna": "Savanna",
                   "Arid & semi-arid": "Arid", "Temperate & subtropical": "Temperate"}[GROUP[s]]
    cell.append([s, str(r.site_name)[:16], r.veg, short_group, f"{r.lat:.1f}",
                 f"{r.rain_mm_yr:.0f}", f"{r.pet_mm_yr:.0f}", f"{r.aridity:.2f}",
                 f"{r.lai_mean:.1f}", f"{int(r.years)}", f"{int(r.n_events)}",
                 f"{eb.loc[s, 'eb_gap_pct']:.0f}", f"{r.ET_median:.0f}"])
tab = ax.table(cellText=cell, colLabels=cols, loc="upper center", cellLoc="center",
               colWidths=[0.07, 0.13, 0.05, 0.09] + [0.074] * 9)
tab.auto_set_font_size(False)
tab.set_fontsize(7.5)
tab.scale(1, 1.35)
for (ri, ci), cl in tab.get_celld().items():
    if ri == 0:
        cl.set_facecolor("#dde4f0")
        cl.set_text_props(weight="bold")
    elif cell[ri - 1][0] in SHORT:
        cl.set_facecolor("#f4f4f4")
ax.set_title("Table 1. The 22 sites, ordered from wettest to driest (aridity = PET / rain). "
             "Grey rows: fewer than 3 complete years.", fontsize=10, loc="left")
caption(fig, "Energy gap % = how much more energy the tower's energy-balance-corrected fluxes "
        "contain than the raw ones (Qle_cor + Qh_cor vs Qle + Qh, measured steps): a larger "
        "gap means the tower misses more energy. LAI = mean satellite leaf area index from the "
        "tower file. Tower ET = metric 1 (gap-filled tower).", y=0.075)
save_fig(pdf, w, fig)

# ---- Fig 1: map and climate ----
fig, axes = plt.subplots(1, 2, figsize=A4L)
ax = axes[0]
for s in o.index:
    ax.scatter(o.loc[s, "lon"], o.loc[s, "lat"], s=40 + 25 * o.loc[s, "aridity"],
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
caption(fig, "The sites cover monsoon tropics (north), wet tropical rainforest (north-east), the "
        "arid interior and temperate south. Sites left of the dashed line receive less rain "
        "than the atmosphere could evaporate (water-limited); the rainforests on the right are "
        "energy-limited. There is no coastline: positions are plotted as raw latitude/longitude.")
fig.tight_layout(rect=(0, 0.12, 1, 0.95))
save_fig(pdf, w, fig)

# ---- Fig 2: model vs obs scatter for each metric ----
fig, axes = plt.subplots(2, 3, figsize=A4L)
for ax, (k, label) in zip(axes.ravel(), METRICS):
    sel = o.index if k != "CV_annual" else o.index[o.years >= 3]
    lo = min(o.loc[sel, k].min(), c.loc[sel, k].min(), j.loc[sel, k].min())
    hi = max(o.loc[sel, k].max(), c.loc[sel, k].max(), j.loc[sel, k].max())
    pad = 0.05 * (hi - lo)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "k--", lw=0.7)
    for mod, m in [("CABLE", c), ("JULES", j)]:
        ax.scatter(o.loc[sel, k], m.loc[sel, k], color=MCOL[mod], s=22, alpha=0.8,
                   label=f"{mod} (r = {stats.pearsonr(o.loc[sel, k], m.loc[sel, k])[0]:.2f})")
    for s in sel:
        if abs(c.loc[s, k] - o.loc[s, k]) > 0.25 * (hi - lo) or abs(j.loc[s, k] - o.loc[s, k]) > 0.25 * (hi - lo):
            ax.annotate(s.replace("AU-", ""), (o.loc[s, k], max(c.loc[s, k], j.loc[s, k])),
                        fontsize=6.5, xytext=(2, 2), textcoords="offset points")
    ax.set_title(label)
    ax.set_xlabel("tower")
    ax.set_ylabel("model")
    ax.legend(fontsize=7)
fig.suptitle("Figure 2. Model vs tower for each signature (one point per site; dashed = 1:1)",
             fontsize=12, x=0.06, ha="left")
caption(fig, "Points on the dashed line are perfect. Labelled sites are the largest misses. Annual "
        "ET and WS are well ordered across sites; R and CV_annual much less so. CV_annual uses "
        "only the 17 sites with at least 3 complete years.")
fig.tight_layout(rect=(0, 0.11, 1, 0.95))
save_fig(pdf, w, fig)

# ---- Table 2: statistics and head-to-head ----
fig, axes = plt.subplots(2, 1, figsize=A4L, gridspec_kw={"height_ratios": [1.4, 1]})
ax = axes[0]
ax.axis("off")
t = stats_tab.copy()
cell = [[r.Model, r.Metric, r.n, f"{r['r (across sites)']:.2f}", f"{r['median error']:+.3g}",
         f"{r['mean error']:+.3g}", f"{r.RMSE:.3g}"] for _, r in t.iterrows()]
tb = ax.table(cellText=cell, colLabels=["Model", "Metric", "n", "r across sites",
                                        "median error", "mean error", "RMSE"],
              loc="upper center", cellLoc="center",
              colWidths=[0.08, 0.3, 0.05, 0.13, 0.14, 0.14, 0.12])
tb.auto_set_font_size(False)
tb.set_fontsize(8)
tb.scale(1, 1.25)
for (ri, ci), cl in tb.get_celld().items():
    if ri == 0:
        cl.set_facecolor("#dde4f0")
        cl.set_text_props(weight="bold")
    elif cell[ri - 1][0] == "JULES":
        cl.set_facecolor("#fbeaea")
ax.set_title("Table 2a. How well each model reproduces each signature across the 22 sites "
             "(error = model - tower)", fontsize=10, loc="left")
ax = axes[1]
ax.axis("off")
cell = [[r.Metric, r.n, r["JULES closer"], r["CABLE closer"], f"{r['median |err| CABLE']:.3g}",
         f"{r['median |err| JULES']:.3g}", f"{r['Wilcoxon p']:.3f}"] for _, r in h2h.iterrows()]
cell.append(["4 TSP (exact peak month)", 22, int((tsp_j == 0).sum()), int((tsp_c == 0).sum()),
             f"{tsp_c.abs().median():.1f} months", f"{tsp_j.abs().median():.1f} months", "-"])
tb = ax.table(cellText=cell, colLabels=["Metric", "n", "sites JULES closer", "sites CABLE closer",
                                        "median |err| CABLE", "median |err| JULES", "Wilcoxon p"],
              loc="upper center", cellLoc="center",
              colWidths=[0.3, 0.05, 0.12, 0.12, 0.14, 0.14, 0.1])
tb.auto_set_font_size(False)
tb.set_fontsize(8)
tb.scale(1, 1.3)
for (ri, ci), cl in tb.get_celld().items():
    if ri == 0:
        cl.set_facecolor("#dde4f0")
        cl.set_text_props(weight="bold")
    elif ci == 6 and cell[ri - 1][6] != "-" and float(cell[ri - 1][6]) < 0.05:
        cl.set_facecolor("#c7e9c0")
ax.set_title("Table 2b. Head to head: at how many sites is each model closer to the tower? "
             "\nGreen = significant at p < 0.05. For TSP the counts are sites with the exact peak month.",
             fontsize=10, loc="left")
fig.tight_layout()
save_fig(pdf, w, fig)

# ---- Fig 3: Budyko ----
fig, ax = plt.subplots(figsize=(9, 7))
phi = np.linspace(0.2, 4.6, 200)
budyko = np.sqrt(phi * np.tanh(1 / phi) * (1 - np.exp(-phi)))
ax.plot(phi, budyko, color="0.5", lw=1.2, label="Budyko curve")
ax.plot([0, 1], [0, 1], color="0.7", lw=0.8, ls=":")
ax.axhline(1, color="0.7", lw=0.8, ls=":")
for s in order:
    x = o.loc[s, "aridity"]
    ys = [o.loc[s, "ET_median"] / o.loc[s, "rain_mm_yr"], c.loc[s, "ET_median"] / o.loc[s, "rain_mm_yr"],
          j.loc[s, "ET_median"] / o.loc[s, "rain_mm_yr"]]
    ax.plot([x, x], [min(ys), max(ys)], color="0.8", lw=1, zorder=1)
    for y, mod, mk in zip(ys, ["obs", "CABLE", "JULES"], ["s", "o", "^"]):
        ax.scatter(x, y, color=MCOL[mod], marker=mk, s=35, zorder=3)
    ax.annotate(s.replace("AU-", ""), (x, ys[0]), fontsize=7, xytext=(4, -9), textcoords="offset points")
for mod, mk in [("obs", "s"), ("CABLE", "o"), ("JULES", "^")]:
    ax.scatter([], [], color=MCOL[mod], marker=mk, label=mod if mod != "obs" else "tower")
ax.set_xlabel("aridity index  PET / rain")
ax.set_ylabel("evaporative index  ET / rain")
ax.set_xlim(0, 4.6)
ax.set_ylim(0, 1.4)
ax.legend(fontsize=8, loc="lower right")
ax.set_title("Figure 3. Budyko view: how much of the rain each site and model evaporates", loc="left")
caption(fig, "Each site has one x position (its climate) and three points (tower, CABLE, JULES), "
        "joined by a grey line. The Budyko curve is the classic expectation for catchments. Wet "
        "sites (left) are under-evaporated by both models, most at the rainforests and AU-How/"
        "AU-Lit. Dry sites (right) sit near ET/rain = 1 for everyone, because rain is the cap. "
        "Points above 1 (e.g. AU-Whr tower, JULES at several arid sites) evaporate more than it "
        "rains, which needs stored or extra water.", y=0.115, width=120)
fig.tight_layout(rect=(0, 0.16, 1, 1))
save_fig(pdf, w, fig)

# ---- Fig 4: errors vs aridity ----
fig, axes = plt.subplots(2, 2, figsize=A4L)
panels = [("Annual ET relative error", lambda m: rel_err(m), "(model - tower) / tower"),
          ("WS error", lambda m: m.WS - o.WS, "model - tower"),
          ("R error (event responsiveness)", lambda m: m.R - o.R, "model - tower"),
          ("Dry/wet-season ratio error", lambda m: m.dry_wet_ratio - o.dry_wet_ratio, "model - tower")]
for ax, (title, fn, ylab) in zip(axes.ravel(), panels):
    for mod, m in [("CABLE", c), ("JULES", j)]:
        e = fn(m)
        rho = stats.spearmanr(o.aridity, e)
        ax.scatter(o.aridity, e, color=MCOL[mod], s=22, alpha=0.8,
                   label=f"{mod} (rho = {rho.statistic:.2f}, p = {rho.pvalue:.2f})")
    for s in o.index:
        ax.annotate(s.replace("AU-", ""), (o.loc[s, "aridity"], fn(c)[s]), fontsize=6,
                    color="0.4", xytext=(2, 2), textcoords="offset points")
    ax.axhline(0, color="k", lw=0.7)
    ax.set_title(title)
    ax.set_xlabel("aridity  PET / rain")
    ax.set_ylabel(ylab)
    ax.legend(fontsize=7)
fig.suptitle("Figure 4. Do the errors depend on how dry the site is?", fontsize=12, x=0.06, ha="left")
caption(fig, "rho = Spearman rank correlation of error with aridity across the 22 sites. Labels "
        "are placed at CABLE's point. Annual ET goes from under-estimated at wet sites to "
        "slightly over-estimated at dry sites (significant for CABLE). Excess rain "
        "responsiveness (R) grows clearly with dryness for both models.")
fig.tight_layout(rect=(0, 0.11, 1, 0.95))
save_fig(pdf, w, fig)

# ---- Fig 5: group summary ----
fig, axes = plt.subplots(1, 5, figsize=A4L, sharey=True)
groups = ["Tropical rainforest", "Tropical savanna", "Temperate & subtropical", "Arid & semi-arid"]
for ax, col in zip(axes, grp["CABLE"].columns):
    y = np.arange(len(groups))
    ax.barh(y + 0.2, [grp["CABLE"].loc[g, col] for g in groups], 0.38, color=MCOL["CABLE"], label="CABLE")
    ax.barh(y - 0.2, [grp["JULES"].loc[g, col] for g in groups], 0.38, color=MCOL["JULES"], label="JULES")
    ax.axvline(0, color="k", lw=0.7)
    ax.set_yticks(y, [g.replace(" & ", " &\n") for g in groups])
    ax.set_title(col, fontsize=9)
axes[0].legend(fontsize=7)
fig.suptitle("Figure 5. Median error by site group (model - tower)", fontsize=12, x=0.06, ha="left")
caption(fig, "Annual ET relative error is a fraction (-0.3 = 30% too little). The rainforests stand "
        "out for annual ET and WS; the savannas for the dry/wet ratio and CV_monthly (CABLE "
        "especially); the arid sites for R. Groups have 3-8 sites, so these are indicative.")
fig.tight_layout(rect=(0, 0.11, 1, 0.93))
save_fig(pdf, w, fig)

# ---- Fig 6: seasonal cycles small multiples ----
fig, axes = plt.subplots(6, 4, figsize=A4, sharex=True)
for ax, s in zip(axes.ravel(), order):
    cs = clim[clim.site == s]
    rain = cs[cs.source == "obs"].sort_values("month")
    ax2 = ax.twinx()
    ax2.bar(rain.month, rain.rain, color="#c6dbef", width=0.8, zorder=0)
    ax2.set_yticks([])
    ax.set_zorder(ax2.get_zorder() + 1)
    ax.patch.set_visible(False)
    for src in ["obs", "CABLE", "JULES"]:
        t = cs[cs.source == src].sort_values("month")
        ax.plot(t.month, t.ET, color=MCOL[src], lw=1.6 if src == "obs" else 1.1)
    ax.set_title(f"{s} ({o.loc[s, 'veg']}, PET/P {o.loc[s, 'aridity']:.1f})", fontsize=7.5)
    ax.tick_params(labelsize=6.5)
    ax.set_xticks([1, 4, 7, 10], ["J", "A", "J", "O"])
for ax in axes.ravel()[len(order):]:
    ax.axis("off")
ax = axes.ravel()[-1]
for src in ["obs", "CABLE", "JULES"]:
    ax.plot([], [], color=MCOL[src], label="tower" if src == "obs" else src)
ax.bar([], [], color="#c6dbef", label="rain (scaled)")
ax.legend(loc="center", fontsize=8)
fig.suptitle("Figure 6. Mean seasonal cycle of ET (mm/month), wettest site first", fontsize=11, x=0.06, ha="left")
fig.text(0.06, 0.035, "\n".join(textwrap.wrap(
    "Blue bars: mean monthly rain (own scale). Note the dry-season gap between the tower and "
    "CABLE at the northern savannas (How, DaS, Dry, Lit), the dry-season peak at Gin that "
    "neither model has, and the rainforest under-estimate in every month.", 120)),
    fontsize=7.5, style="italic", va="top")
fig.tight_layout(rect=(0, 0.045, 1, 0.97))
save_fig(pdf, w, fig)

# ---- Fig 7: TSP offsets and error heatmap ----
fig, axes = plt.subplots(1, 2, figsize=A4L, gridspec_kw={"width_ratios": [1, 1.6]})
ax = axes[0]
y = np.arange(len(order))
ax.scatter(tsp_c[order], y + 0.15, color=MCOL["CABLE"], label="CABLE", s=30)
ax.scatter(tsp_j[order], y - 0.15, color=MCOL["JULES"], label="JULES", s=30, marker="^")
ax.axvline(0, color="k", lw=0.8)
ax.set_yticks(y, [f"{s} (tower peak: {int(o.loc[s, 'TSP'])})" for s in order], fontsize=7)
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
    for mod, m in [("C", c), ("J", j)]:
        mat.append((fn(m)[order] / sc).values)
        labels.append(f"{name}\n{mod}")
mat = np.array(mat).T
im = ax.imshow(mat, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
ax.set_xticks(range(len(labels)), labels, fontsize=6.5)
ax.set_yticks(range(len(order)), order, fontsize=7)
for i in range(mat.shape[0]):
    for k in range(mat.shape[1]):
        ax.text(k, i, f"{mat[i, k] * cols[k // 2][2]:.2f}", ha="center", va="center", fontsize=5.5)
ax.set_title("Errors per site (model - tower); C = CABLE, J = JULES")
fig.colorbar(im, ax=ax, shrink=0.6, label="error / scale (red = model too high)")
fig.suptitle("Figure 7. Site-by-site detail, wettest site at the top", fontsize=12, x=0.06, ha="left")
caption(fig, "Left: CABLE's seasonal peak is early or on time almost everywhere (it is never more "
        "than one month late); JULES is exactly right at most sites but far off at a few short "
        "or unusual records (Emr, Lit, Gin). Right: numbers are the errors; colours are scaled "
        "per column so each metric is readable. At the top (rainforests), annual ET is blue (too "
        "little) and WS is red (too much water stress).")
fig.tight_layout(rect=(0, 0.11, 1, 0.95))
save_fig(pdf, w, fig)

# ---- Fig 8: energy closure ----
fig, axes = plt.subplots(1, 2, figsize=A4L)
ax = axes[0]
for mod, m in [("CABLE", c), ("JULES", j)]:
    rho = stats.spearmanr(eb.loc[o.index, "eb_gap_pct"], rel_err(m))
    ax.scatter(eb.loc[o.index, "eb_gap_pct"], rel_err(m), color=MCOL[mod], s=25,
               label=f"{mod} (rho = {rho.statistic:.2f}, p = {rho.pvalue:.3f})")
for s in o.index:
    ax.annotate(s.replace("AU-", ""), (eb.loc[s, "eb_gap_pct"], rel_err(j)[s]), fontsize=6.5,
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
ax.plot(x, c.loc[order, "ET_median"], "o-", color=MCOL["CABLE"], ms=4, label="CABLE")
ax.plot(x, j.loc[order, "ET_median"], "^-", color=MCOL["JULES"], ms=4, label="JULES")
ax.set_xticks(x, [s.replace("AU-", "") for s in order], rotation=70, fontsize=7)
ax.set_ylabel("median annual ET (mm/yr)")
ax.set_title("Annual ET: raw vs corrected tower, wettest to driest")
ax.legend(fontsize=7)
fig.suptitle("Figure 8. Part of the 'model error' is tower error", fontsize=12, x=0.06, ha="left")
caption(fig, "Towers typically measure less energy than the radiation budget says is available. "
        "Where that gap is large (Emr, Cum, Wrr, Rig, Rob, Sam) the raw tower ET is probably too "
        "low, and JULES's apparent over-estimates there shrink or reverse against the corrected "
        "ET. Against corrected ET both models under-estimate at most wet and mid sites.")
fig.tight_layout(rect=(0, 0.11, 1, 0.94))
save_fig(pdf, w, fig)

# ---- Detailed findings ----
w.new_page()
w.heading("2. Findings in detail")
w.para("2.1 Annual totals (metric 1)", size=10.5)
w.para(f"Both models order the sites correctly: annual ET correlates with the tower at r = "
       f"{stats.pearsonr(o.ET_median, c.ET_median)[0]:.2f} (CABLE) and "
       f"{stats.pearsonr(o.ET_median, j.ET_median)[0]:.2f} (JULES). This is largely because "
       f"climate sets most of the between-site differences, and both models are forced with the "
       f"real weather. The errors are not random, though. Relative error rises with aridity "
       f"(CABLE rho = {rho_et_c.statistic:.2f}, p = {rho_et_c.pvalue:.2f}; JULES rho = "
       f"{rho_et_j.statistic:.2f}, p = {rho_et_j.pvalue:.2f}) and falls with leaf area "
       f"(CABLE rho = {rho_lai_c.statistic:.2f}, JULES rho = {rho_lai_j.statistic:.2f}): the "
       f"leafiest, wettest sites are under-evaporated, the dry sparse sites roughly right or "
       f"slightly over. The worst cases are the lowland tropical rainforests AU-Cow and AU-Ctr (-26% to -35%), AU-Lit "
       f"(CABLE -40%), AU-How (CABLE -34%), AU-Tum (both -26% to -29%) and, in the other "
       f"direction, AU-Sam (both +55% to +60%) and AU-Emr (CABLE +49%, JULES +18%).")
w.para("2.2 Shared errors", size=10.5)
w.para(f"CABLE's and JULES's annual-ET errors are strongly correlated across sites (rho = "
       f"{rho_cj.statistic:.2f}). At 8 sites both under-estimate and at 5 both over-estimate. "
       f"Two very different models failing the same way at the same sites suggests a common "
       f"cause. Candidates: (a) a process both lack, such as deep roots or groundwater; (b) "
       f"PLUMBER2's out-of-the-box vegetation and soil parameters; (c) problems in the tower "
       f"data itself (Section 2.6). Shared errors are exactly the 'many models fail' signal "
       f"Gab described as most valuable.")
w.para("2.3 Seasonality (metrics 4 and 5, and the dry/wet ratio)", size=10.5)
w.para(f"Peak timing is mostly right: JULES matches the tower's peak month at "
       f"{int((tsp_j == 0).sum())} sites and CABLE at {int((tsp_c == 0).sum())}; within one "
       f"month, {int((tsp_j.abs() <= 1).sum())} and {int((tsp_c.abs() <= 1).sum())}. CABLE's "
       f"misses are almost all early (mean offset {tsp_c.mean():+.1f} months): it evaporates "
       f"hardest early in the wet season and fades before the tower does. The seasonal swing "
       f"(CV_monthly) is too big in both models (median +{(c.CV_monthly - o.CV_monthly).median():.2f} "
       f"CABLE, +{(j.CV_monthly - o.CV_monthly).median():.2f} JULES), and for CABLE that excess is "
       f"tied to the dry season: sites where CABLE's dry/wet ratio is too low are those where its "
       f"CV_monthly is too high (rho = -0.58). In other words the exaggerated swing is a "
       f"dry-season problem, not a wet-season one.")
w.para("2.4 Dry-season persistence (extra diagnostic)", size=10.5)
w.para("The tropical savannas show the clearest model failure in this analysis. In the dry "
       "season the tower at AU-How still evaporates about half of its wet-season rate "
       "(ratio 0.53); CABLE drops to 0.26 and JULES to 0.40. AU-DaS (0.43 vs 0.20/0.23) and "
       "AU-Dry (0.41 vs 0.16/0.19) are similar. Savanna trees reach deep soil water through the "
       "dry season; a model with too-shallow roots or too-strict water-stress response runs dry. "
       f"JULES is significantly closer on this ratio (17 of 22 sites, p = "
       f"{h2h.loc[5, 'Wilcoxon p']:.3f}). The same signature appears at AU-Gin, more extremely: "
       f"the tower's dry/wet ratio is 1.27 (more ET in the dry summer), against 0.52 and 0.61 "
       f"for the models. Gingin sits on deep sand above a shallow aquifer, and the vegetation is "
       f"known to use groundwater; neither model can.")
w.para("2.5 Water stress and event responsiveness (metrics 6 and 8)", size=10.5)
w.para(f"Because the tower and both models share the same PET, WS is almost the mirror image "
       f"of the annual-ET error (rho = {rho_ws_c.statistic:.2f} for CABLE, "
       f"{rho_ws_j.statistic:.2f} for JULES), so it adds little new information here. It does "
       f"expose the PET substitute: at AU-Cow and AU-Ctr the tower evaporates more than our "
       f"Priestley-Taylor PET (negative WS), which is physically possible for wet canopies and "
       f"warm advected air but suggests PET is low there. R is the most interesting dynamic "
       f"signature. At wet sites R is near zero for everyone. At dry sites the tower's R stays "
       f"modest (0.10-0.34), but both models' R climbs to 0.4-0.55: modelled ET jumps after "
       f"every big storm and then falls, while the real ecosystem buffers rain in deeper soil "
       f"and releases it slowly. The excess grows with aridity (CABLE rho = "
       f"{rho_r_c.statistic:.2f}, p = {rho_r_c.pvalue:.3f}; JULES rho = {rho_r_j.statistic:.2f}, "
       f"p = {rho_r_j.pvalue:.3f}).")
w.para("2.6 The tower data are part of the answer", size=10.5)
w.para(f"The raw tower fluxes miss between -7% and +65% of the available energy (Table 1). "
       f"JULES's relative error correlates with that gap (rho = {rho_eb_j.statistic:.2f}, p = "
       f"{rho_eb_j.pvalue:.3f}; CABLE rho = {rho_eb_c.statistic:.2f}): the sites where JULES "
       f"seems to over-evaporate (AU-Sam, AU-Emr, AU-Wrr, AU-Cum) are among those with the "
       f"largest gaps, where the raw tower ET is probably too low. Against "
       f"energy-balance-corrected ET, both models under-evaporate (median "
       f"{100 * relc_cor.median():.0f}% CABLE, {100 * relj_cor.median():.0f}% JULES). The paper "
       f"itself uses the raw (uncorrected) OzFlux ET, so our main results follow it, but any "
       f"claim about the sign of a model bias should be checked against both versions.")

# ---- Site by site notes ----
notes = {
    "AU-Cow": "Wet rainforest. Both models 26-35% too little ET in every month; tower ET above PT PET.",
    "AU-Ctr": "Wet rainforest. Largest absolute miss (-440 mm/yr CABLE). CABLE peaks a month early.",
    "AU-Wrr": "Tall wet forest, 2 years only. JULES +36% ET, CABLE about right; WS very low in JULES.",
    "AU-Rob": "Upland rainforest, 37% energy gap. Both ~12-17% low; peaks 1-2 months early.",
    "AU-Tum": "Tall temperate forest. Both 26-29% too little ET, worse against corrected ET; summer peak too sharp.",
    "AU-Sam": "Subtropical pasture. Both +55-60% ET; tower ET only half the rain and 37% energy gap: suspect tower.",
    "AU-DaP": "Monsoon pasture. Totals within 3-11%; dry-season ET a little too high in both models.",
    "AU-Lit": "Savanna, 2 years. CABLE -40% ET; JULES peak 4 months late (short record).",
    "AU-Otw": "Temperate pasture, 2 years. Close on everything; a good 'control' site.",
    "AU-How": "Savanna, 15 years. CABLE -34%, JULES -18%; dry-season collapse (0.26/0.40 vs 0.53). Key case.",
    "AU-DaS": "Savanna. JULES total right; both collapse in dry season (0.20/0.23 vs 0.43).",
    "AU-Cum": "Dry sclerophyll forest, 43% energy gap. CABLE -14%, JULES +4%; R far too low in CABLE.",
    "AU-Dry": "Savanna. Dry-season collapse (0.16/0.19 vs 0.41); seasonal swing far too big.",
    "AU-Rig": "Temperate grassland. Totals right; interannual variability 2-3x too high in both.",
    "AU-Stp": "Tussock grassland. Totals within ~12%; dry-season ET too low; swing too big.",
    "AU-Emr": "Semi-arid grassland, 2 years, 65% energy gap. CABLE +49%; JULES peak 4 months off.",
    "AU-Gin": "Banksia woodland on sand. Tower peaks in dry summer (groundwater); neither model can.",
    "AU-Whr": "Box woodland, 2 years. Tower ET above rain (1.25); models R too high.",
    "AU-TTE": "Arid shrubland. Very close on totals, seasonality and WS; R slightly high.",
    "AU-GWW": "Semi-arid woodland. Both +24% ET and R ~0.54 vs 0.21: too responsive to rain.",
    "AU-ASM": "Mulga. Totals close (JULES +32%); interannual variability too low; JULES R high.",
    "AU-Cpr": "Mallee. Totals exact; R 0.47/0.49 vs 0.10: ET reacts to every storm in the models.",
}
fig, ax = plt.subplots(figsize=A4)
ax.axis("off")
cell = [[s, GROUP[s].replace(" & ", "/"), "\n".join(textwrap.wrap(notes[s], 78))] for s in order]
tb = ax.table(cellText=cell, colLabels=["Site", "Group", "What stands out"], loc="upper center",
              cellLoc="left", colWidths=[0.1, 0.17, 0.73])
tb.auto_set_font_size(False)
tb.set_fontsize(7.2)
for (ri, ci), cl in tb.get_celld().items():
    cl.set_height(0.038 if ri else 0.025)
    if ri == 0:
        cl.set_facecolor("#dde4f0")
        cl.set_text_props(weight="bold")
ax.set_title("Table 3. Site by site, wettest to driest", fontsize=11, loc="left")
fig.tight_layout()
save_fig(pdf, w, fig)

# ---- Caveats and next steps ----
w.new_page()
w.heading("3. Caveats")
w.bullets([
    "PET is a substitute (Priestley-Taylor from tower net radiation) for the paper's Morton "
    "PET. It is likely too low at the wet tropical sites (tower ET exceeds it). WS and the "
    "aridity ordering depend on it; the other metrics do not.",
    "Tower ET is gap-filled (as in the paper), and 11-51% of half-hours at each site are "
    "gap-filled. The raw tower data also miss energy (Section 2.6).",
    f"Five sites have only 2 complete years ({', '.join(SHORT)}). Their CV_annual is excluded "
    "and their other metrics are noisy.",
    "R uses a 10 mm/day rain threshold; the paper does not state its value. Events less than "
    "10 days apart share their windows.",
    "22 sites is a small sample. Wilcoxon and Spearman p-values are indicative; with 6 metrics "
    "tested, one 'significant' result could be chance.",
    "The model runs are PLUMBER2 'out of the box' runs. Which CABLE and JULES configurations "
    "these are (CABLE vs CABLE-POP, JULES GL9 vs GL9_LAI) is still to be confirmed.",
    "AU-Ync is excluded: its tower files contain AU-ASM data. Note that et_signatures.py "
    "currently loops over every file in data/harmonised, including AU-Ync.",
])
w.heading("4. What this suggests doing next")
w.bullets([
    "Get Morton PET from SILO for each site and rerun, to remove the PET caveat and match the paper.",
    "Rerun every metric against energy-balance-corrected tower ET and report both; it changes "
    "the sign of JULES's bias at several sites.",
    "Use the dry-season result as a case study: the AU-How, AU-DaS and AU-Dry dry seasons, "
    "with soil moisture and transpiration from the CABLE and JULES outputs, to test whether "
    "root depth or the water-stress function is the cause.",
    "AU-Gin as a groundwater test case, and the arid sites (AU-Cpr, AU-GWW) for the "
    "too-responsive-to-rain behaviour: how fast does modelled ET fall after a storm compared "
    "with the tower (dry-down time scale)?",
    "Check AU-Sam and AU-Emr tower data quality before treating them as model failures.",
    "Repeat the signatures for the other 11 PLUMBER2 models to see which errors are shared "
    "by many models (a general process gap) and which are specific to CABLE or JULES.",
])
w.close()
pdf.close()
print(f"wrote {HERE / 'et_signatures_report.pdf'} ({w.page} pages)")


# ==========================================================================
# Extra figures (PNG only, not in the PDF): error maps and WS time series
# ==========================================================================
# These are saved straight into hydrology/figures/:
#   map_<metric>.png    four maps of Australia per signature, one dot per site
#   ws_timeseries.png   water stress year by year, obs vs CABLE vs JULES

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


for col, title, units in MAP_METRICS:
    # Which sites to show. CV_annual needs at least 3 complete years,
    # as everywhere else in this report.
    sites = o.index if col != "CV_annual" else o.index[o.years >= 3]
    lon, lat = o.loc[sites, "lon"], o.loc[sites, "lat"]

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


# --------------------------------------------------------------------------
# Water stress (WS) year by year. WS = 1 - ET/PET for each complete year,
# the same formula as signature 6 but for one year at a time (from
# collect_results.py, et_annual.csv). Sites go wettest to driest.
# --------------------------------------------------------------------------
annual = pd.read_csv(OUT / "et_annual.csv")
fig, axes = plt.subplots(6, 4, figsize=(11.69, 14), sharey=True)
for ax, site in zip(axes.flat, order):
    for source in ["obs", "CABLE", "JULES"]:
        one = annual[(annual.site == site) & (annual.source == source)]
        ax.plot(one.year, one.WS, marker="o", markersize=3, color=MCOL[source], label=source)
    ax.axhline(0, color="0.7", linewidth=0.6)
    ax.set_title(f"{site} ({GROUP[site]})", fontsize=8)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True, nbins=4))
    ax.ticklabel_format(axis="x", useOffset=False)     # show 2016, not "0 + 2.016e3"
    ax.tick_params(labelsize=7)
for ax in axes.flat[len(order):]:
    ax.axis("off")                     # 22 sites, 24 panels: hide the spare two
for ax in axes[:, 0]:
    ax.set_ylabel("WS = 1 - ET/PET")
axes[0, 0].legend(fontsize=7)
fig.suptitle("Water stress each year, wettest site to driest\n"
             "(0 = ET keeps up with PET, 1 = almost no ET)", fontsize=12, weight="bold")
fig.tight_layout()
fig.savefig(FIGS / "ws_timeseries.png", dpi=200)
plt.close(fig)
print(f"wrote {FIGS / 'ws_timeseries.png'}")


# --------------------------------------------------------------------------
# Water stress month by month: WS = 1 - ET/PET for each complete month
# (et_monthly.csv from collect_results.py). 12x more points than the yearly
# plot, so you can see the seasons and individual droughts. In winter PET is
# small, so WS swings far below zero when ET is larger than PET; the y axis is
# cut at -1 and those months run off the bottom of the panel.
# --------------------------------------------------------------------------
monthly = pd.read_csv(OUT / "et_monthly.csv", parse_dates=["month"])
fig, axes = plt.subplots(6, 4, figsize=(11.69, 14), sharey=True)
for ax, site in zip(axes.flat, order):
    for source in ["obs", "CABLE", "JULES"]:
        one = monthly[(monthly.site == site) & (monthly.source == source)]
        ws = one.set_index("month").WS
        # Fill in missing (incomplete) months with NaN so the line breaks there
        # instead of joining across the gap.
        all_months = pd.date_range(ws.index.min(), ws.index.max(), freq="MS")
        ws = ws.reindex(all_months)
        ax.plot(ws.index, ws.values, color=MCOL[source], linewidth=0.8, label=source)
    ax.axhline(0, color="0.7", linewidth=0.6)
    ax.set_ylim(-1, 1)
    ax.set_title(f"{site} ({GROUP[site]})", fontsize=8)
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


# --------------------------------------------------------------------------
# The same two water stress plots, but one PNG per site, so each site can be
# looked at on its own at a readable size. Saved into figures/WS_timeseries/:
#   <SITE>_WS_annual.png    WS each year   (from et_annual.csv)
#   <SITE>_WS_monthly.png   WS each month  (from et_monthly.csv)
# --------------------------------------------------------------------------
WS_DIR = FIGS / "WS_timeseries"
WS_DIR.mkdir(exist_ok=True)          # make the subfolder if it isn't there yet

for site in order:
    # Yearly
    fig, ax = plt.subplots(figsize=(8, 4))
    for source in ["obs", "CABLE", "JULES"]:
        one = annual[(annual.site == site) & (annual.source == source)]
        ax.plot(one.year, one.WS, marker="o", markersize=4, color=MCOL[source], label=source)
    ax.axhline(0, color="0.7", linewidth=0.6)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.ticklabel_format(axis="x", useOffset=False)
    ax.set_ylabel("WS = 1 - ET/PET")
    ax.set_title(f"{site} ({GROUP[site]}): water stress each year")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(WS_DIR / f"{site}_WS_annual.png", dpi=200)
    plt.close(fig)

    # Monthly (same gap handling and y limits as the all-sites plot above)
    fig, ax = plt.subplots(figsize=(12, 4))
    for source in ["obs", "CABLE", "JULES"]:
        one = monthly[(monthly.site == site) & (monthly.source == source)]
        ws = one.set_index("month").WS
        all_months = pd.date_range(ws.index.min(), ws.index.max(), freq="MS")
        ws = ws.reindex(all_months)
        ax.plot(ws.index, ws.values, color=MCOL[source], linewidth=0.9, label=source)
    ax.axhline(0, color="0.7", linewidth=0.6)
    ax.set_ylim(-1, 1)
    ax.xaxis.set_major_locator(matplotlib.dates.YearLocator())
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y"))
    ax.set_ylabel("WS = 1 - ET/PET")
    ax.set_title(f"{site} ({GROUP[site]}): water stress each month (y axis cut at -1)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(WS_DIR / f"{site}_WS_monthly.png", dpi=200)
    plt.close(fig)

print(f"wrote {2 * len(order)} per-site WS plots to {WS_DIR}")

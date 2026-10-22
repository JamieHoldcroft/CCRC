import et_signatures as ets     
import xarray as xr    
import pandas as pd 
import argparse

# parser = argparse.ArgumentParser()
# parser.add_argument("var", nargs="*")

# args =parser.parse_args()
# VARIABLES = args.var or ["Qle", "Qh", "Precip", "Rnet"]

# rows = []
# for path in sorted(ets.HARMONISED.glob("*.nc")):
#     ds = xr.open_dataset(path)
#     for var in VARIABLES:
#         gap = ~ets.is_measured(ds, var)
#         t = gap.index
#         rows.append({
#             "site": path.stem, "var": var,
#             "pct_gapfilled": 100 * gap.mean(),
#             "pct_day": 100 * gap[(t.hour >= 6) & (t.hour < 18)].mean(),
#             "pct_night": 100 * gap[(t.hour < 6) | (t.hour >= 18)].mean(),
#             "day/night": (100 * gap[(t.hour >= 6) & (t.hour < 18)].mean()) / (100 * gap[(t.hour < 6) | (t.hour >= 18)].mean())
#         })
# pd.DataFrame(rows).to_csv("hydrology/gap_filled.csv", index=False)

# d = pd.read_csv("hydrology/gap_filled.csv")
# print(d.to_string())

############################################################################
# Days, months and years lost once gap-filled data are removed.
# Uses exactly the same filters as the analysis (et_signatures.py), so the
# results depend on MIN_MEASURED_PER_DAY and MIN_COMPLETE there.
# A day is kept only if the tower, CABLE and JULES ET and PET all have a value
# (ets.analysis_days); a month or year is kept if it is complete (ets.totals).

from pathlib import Path
import matplotlib
matplotlib.use("Agg")                          # save figures without opening windows
import matplotlib.pyplot as plt

PLOT_DIR = Path(__file__).resolve().parent / f"figures{ets.SUFFIX}" / "gap_filled"
PLOT_DIR.mkdir(parents=True, exist_ok=True)

lost_rows = []
for path in sorted(ets.HARMONISED.glob("*.nc")):
    if path.stem in ets.SKIP_SITES:
        continue
    ds = xr.open_dataset(path)
    et_by_source, pet = ets.analysis_days(ds)
    et = et_by_source["obs"]                       # every source has the same days

    days_total = len(et)                           # every day in the record
    months_total = len(et.resample("MS").size())   # every month in the record
    years_total = len(et.resample("YS").size())    # every year in the record (first and last may be partial)

    days_kept = et.notna().sum()
    months_kept = len(ets.totals(et, "MS"))
    years_kept = len(ets.totals(et, "YS"))

    lost_rows.append({
        "site": path.stem,
        "days_lost_%": 100 * (1 - days_kept / days_total),
        "months_lost_%": 100 * (1 - months_kept / months_total),
        "years_lost_%": 100 * (1 - years_kept / years_total),
        "days_kept": days_kept, "months_kept": months_kept,
        "years_kept": years_kept, "years_total": years_total,
    })

    # # Plot: daily tower ET over the whole record. Grey = every day, gap-filled
    # # data included. Black dots = the days the analysis keeps (measured data only).
    # step = ds.attrs["timestep_s"]
    # lam = ets.latent_heat(ds.Tair_C.to_series())
    # et_all = (ds.Qle.sel(source="obs").to_series() / lam * step).resample("1D").sum(min_count=1)

    # fig, ax = plt.subplots(figsize=(14, 4))
    # ax.plot(et_all.index, et_all.values, color="0.75", linewidth=0.6, label="all days (gap-filled included)")
    # ax.plot(et.index, et.values, color="k", linewidth=0.8, label="days kept (measured only)")
    # # ax.plot(et.index, et.values, "o", color="k", markersize=1.5, label="days kept (measured only)")
    # ax.set_ylabel("tower ET (mm/day)")
    # ax.set_title(f"{path.stem}: {days_kept} of {days_total} days kept "
    #              f"({100 * days_kept / days_total:.0f}%), {years_kept} complete years")
    # ax.legend(loc="upper right")
    # ax.grid(alpha=0.3)
    # fig.tight_layout()
    # fig.savefig(PLOT_DIR / f"{path.stem}_ET_kept.png", dpi=150)
    # plt.close(fig)

lost = pd.DataFrame(lost_rows)
print()
print(f"Lost after removing gap-filled data (MIN_MEASURED_PER_DAY = {ets.MIN_MEASURED_PER_DAY}, "
      f"MIN_COMPLETE = {ets.MIN_COMPLETE})")
print(lost.round(1).to_string(index=False))








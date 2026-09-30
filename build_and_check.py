# run from root and saves output .png to /outputs.
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ccrc.events import _daily, _seasonal_percentile, add_event_flags, summarise_events
from ccrc.loaders import combine, load_jules, load_site


met = 'data/AU-Cum_2013-2018_OzFlux_Met.nc'
flux = 'data/AU-Cum_2013-2018_OzFlux_Flux.nc'

obs = load_site(met, flux)

jules = load_jules('data/JULES/local_AU-Cum_fluxnet2015_GAL9.AU-Cum.nc')

# list allows for any number of models
ds = add_event_flags(combine(obs, [jules]))
print(ds.time.values[0])
print("\nEvents:\n", summarise_events(ds))
for name in ["drought", "heatwave", "high_vpd"]:
    print(f"{name}: {ds[name].attrs['definition']}")

for var, units in [("Qle", "W/m2"), ("Qh", "W/m2"), ("GPP", "umol/m2/s")]:
    print(f"\nMean {var} ({units}), measured obs only:")
    measured = ds[f"measured_{var}"]
    for name in ["normal", "drought", "heatwave", "high_vpd"]:
        q = ds[var].where(ds[name] & measured).mean("time").astype(float)
        print(f"  {name:9s}", dict(zip(ds.source.values.tolist(), [round(x, 1) for x in q.values.tolist()])))

# Time series with events shaded, to check the flags by eye.
rain90 = _daily(ds.Precip_mm, "sum").rolling(90, min_periods=90).sum()
tmean = _daily(ds.Tair_C, "mean")
vpd = _daily(ds.VPD_kPa.where(ds.SWdown > 50), "mean")
panels = [
    ("drought", "tan", rain90, "90-day rain (mm)",
     _seasonal_percentile(rain90, 10, 15), "seasonal p10"),
    ("heatwave", "salmon", tmean, "Daily mean Tair (degC)",
     pd.Series(np.nanpercentile(tmean, 95), index=tmean.index), "p95 of record"),
    ("high_vpd", "plum", vpd, "Daytime VPD (kPa)",
     pd.Series(np.nanpercentile(vpd, 90), index=vpd.index), "p90 of record"),
]
fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
for ax, (flag, colour, series, label, thresh, thresh_label) in zip(axes, panels):
    ax.plot(series.index, series.values, lw=0.7, color="k")
    ax.plot(thresh.index, thresh.values, lw=1, color="tab:blue", ls="--",
            label=thresh_label)
    f = ds[flag].resample(time="1D").first()
    ax.fill_between(f.time.values, 0, 1, where=f.values, color=colour, alpha=0.6,
                    transform=ax.get_xaxis_transform(), label=flag)
    ax.set_ylabel(label)
    ax.legend(loc="upper right", fontsize=8)
fig.suptitle("AU-Cum event flags (dashed = threshold)")
fig.tight_layout()
fig.savefig("outputs/au_cum_events.png")
print("\nSaved outputs/au_cum_events.png")

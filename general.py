import os
import numpy as np
import pandas as pd
import xarray as xr

flux_data = 'data/JULES/local_AU-Cum_fluxnet2015_GAL9.AU-Cum.nc'
met_data = 'data/AU-Cum_2013-2018_OzFlux_Met.nc'

ds = xr.open_dataset(flux_data)
print(ds['time'].values[0])
# for name, value in ds.attrs.items():
#     print(f'{name}: {value}')

# ds = ds.squeeze()

# ds.Tair.sel(time='2017-03').plot()

# t = ds.time.to_index()
# print('first / last stamp :', t[0], '/', t[-1])
# print('number of steps    :', len(t))
# print('monotonic, dupes   :', t.is_monotonic_increasing, int(t.duplicated().sum()))
# print('\nspacing between stamps:')
# print(pd.Series(t).diff().value_counts().head())
# print('\nsteps per calendar year (by stamp):')
# print(pd.Series(t.year).value_counts().sort_index().to_string())
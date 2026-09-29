import os
import numpy as np
import pandas as pd
import xarray as xr

flux_data = 'Data/AU-Cum_2013-2018_OzFlux_Flux.nc'
met_data = 'Data/AU-Cum_2013-2018_OzFlux_Met.nc'

ds = xr.open_dataset(met_data)
for name, value in ds.attrs.items():
    print(f'{name}: {value}')
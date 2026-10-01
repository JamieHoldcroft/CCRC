import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# ds = xr.open_dataset("data/harmonised/AU-Cum.nc")

ds = xr.open_dataset("data/models/LSTM/AU-Cum.nc")
ds1 = xr.open_dataset("data/models/1lin/AU-Cum.nc")

print(ds.data_vars)
print(ds1.data_vars)


models = [
    'obs',
    'CABLE',
    'CHTESSEL_1',
    'CLM5',
    'EntTBM',
    'GFDL',
    'JULES',
    'MATSIRO',
    'MuSICA',
    'NoahMP',
    'ORCHIDEE2',
    'ORCHIDEE3',
    'QUINCY',
    'STEMMUS-SCOPE'
]

# --------------------------------------------------
# Select the time period
# --------------------------------------------------

# start = 0
# end = 50000

# # Observations
# obs = ds.sel(source="obs").Qle.isel(time=slice(start, end)).values

# results = []

# for model in models[1:]:

#     # Model Qle
#     model_qle = ds.sel(source=model).Qle.isel(
#         time=slice(start, end)
#     ).values

#     # Remove timesteps where either observation or model is NaN
#     valid = np.isfinite(obs) & np.isfinite(model_qle)

#     obs_valid = obs[valid]
#     model_valid = model_qle[valid]

#     # -----------------------------
#     # Basic statistics
#     # -----------------------------

#     obs_mean = np.mean(obs_valid)
#     model_mean = np.mean(model_valid)

#     # Bias
#     bias = np.mean(model_valid - obs_valid)

#     # Relative bias (%)
#     relative_bias = 100 * bias / obs_mean

#     # Mean Absolute Error
#     mae = np.mean(np.abs(model_valid - obs_valid))

#     # Root Mean Square Error
#     rmse = np.sqrt(
#         np.mean((model_valid - obs_valid) ** 2)
#     )

#     # Correlation coefficient
#     correlation = np.corrcoef(
#         obs_valid,
#         model_valid
#     )[0, 1]

#     # Standard deviation
#     obs_std = np.std(obs_valid)
#     model_std = np.std(model_valid)

#     # Min / max
#     obs_min = np.min(obs_valid)
#     obs_max = np.max(obs_valid)

#     model_min = np.min(model_valid)
#     model_max = np.max(model_valid)

#     # Number of valid observations
#     n = len(obs_valid)

#     # Store results
#     results.append({
#         "Model": model,
#         "N": n,

#         "Obs Mean (W/m²)": obs_mean,
#         "Model Mean (W/m²)": model_mean,

#         "Bias (W/m²)": bias,
#         "Relative Bias (%)": relative_bias,

#         "MAE (W/m²)": mae,
#         "RMSE (W/m²)": rmse,

#         "Correlation": correlation,

#         "Obs Std (W/m²)": obs_std,
#         "Model Std (W/m²)": model_std,

#         "Obs Min (W/m²)": obs_min,
#         "Obs Max (W/m²)": obs_max,

#         "Model Min (W/m²)": model_min,
#         "Model Max (W/m²)": model_max
#     })


# # --------------------------------------------------
# # Convert to DataFrame
# # --------------------------------------------------

# results_df = pd.DataFrame(results)

# # Display
# print(results_df.to_string(index=False))

# ds = xr.open_dataset("data/flux/AU-Cum_2013-2018_OzFlux_Flux.nc")
# print(ds.data_vars)

# for model in models:
#     print(f"Name: {model}: ")
#     if not (model == "ORCHIDEE2" or model == "ORCHIDEE3"):
#         ds = xr.open_dataset(f"data/models/{model}/AU-Cum.nc")
#         print(ds.data_vars)
#     else:
#         ds = xr.open_dataset(f"data/models/{model}/AU-Cum.nc")
#         print(ds.data_vars)

# flux_data = 'data/JULES/local_AU-Cum_fluxnet2015_GAL9.AU-Cum.nc'
# met_data = 'data/AU-Cum_2013-2018_OzFlux_Met.nc'

# ds = xr.open_dataset(flux_data)
# print(ds['time'].values[0])


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
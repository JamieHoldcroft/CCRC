# Notes:

Getting started code notes:
1. Reading variables from Xarray Dataset: `ds['Tair']`, `ds.data_vars['Tair']` and `ds.Tair` all give the same DataArray. `ds['Tair']` is the safest, because ds.Tair breaks if a variable name clashes with a method such as `mean`.
2. Dot `.` access looks through data variables, then coordinates, then attrs. Ambiguous though so best to use `ds.attrs['timestep_s']` for attributes. Also `.name` exists on a DataArray (`ds.Tair.name` gives `'Tair'`) but not on a Dataset.
3. Dimensions say which axes exist. Coordinates say what each position on an axis means, so you select by label (.sel(source='JULES')) instead of by position number.
4. `ds.Qle.mean('time')` collapses each row into the mean value giving e.g `array([61.92756, 63.79816], dtype=float32)`
5. Same ways to reach a value in dataset
    1. `ds.Qle.sel(source='JULES', time='2013-01-01 00:30')` # by label (preferred)
    2. `ds.Qle.isel(source=1, time=1)` # by position, xarray
    3. `ds.Qle.values[1, 1]` # raw numpy 
6. Models missing flux sites: QUINCY: Ync, ORCHIDEE2: Rob and Wrr, ORCHIDEE3: Rob and Wrr


# General notes
1. Jules starts outputs at 00:30:00 and flux tower observations start 00:00:00 (off by 30 mins)
2. 

## Variables in the combined AU-Cum dataset (`ds`)

Every `X_qc` is the quality flag for `X`: 0 = measured by the tower; 20, 30 or 50 = gap-filled or estimated. The dims column says whether the variable has one value per time step `(time)`, one value per source and time step `(source, time)` for obs vs JULES, or a single number `()`.

| Variable | Category | Dims | Units | What it is |
|---|---|---|---|---|
| CO2air | Weather | time | µmol/mol (ppm) | CO₂ concentration in the air |
| LWdown | Weather / radiation | time | W/m² | Heat (longwave) radiation coming down from sky and clouds |
| SWdown | Weather / radiation | time | W/m² | Sunlight arriving; >50 is used to mean "daytime" |
| Precip | Weather | time | kg/m²/s | Rainfall rate |
| RH | Weather | time | % | Relative humidity |
| Qair | Weather | time | kg/kg | Specific humidity: mass of water vapour per mass of air |
| Tair | Weather | time | K | Air temperature |
| VPD | Weather | time | hPa | Vapour pressure deficit: how dry the air is |
| Wind | Weather | time | m/s | Wind speed |
| Psurf | Weather | time | Pa | Air pressure at the surface |
| latitude, longitude | Site info | () | degrees | Tower location (-33.61, 150.7) |
| reference_height | Site info | () | m | Height of the tower's sensors |
| canopy_height | Site info | () | m | Height of the forest canopy |
| elevation | Site info | () | m | Height above sea level |
| LAI | Vegetation | time | m²/m² | Leaf area index from MODIS satellite: leaf area per ground area |
| LAI_alternative | Vegetation | time | m²/m² | LAI from a second satellite product (Copernicus) |
| Resp | Carbon (obs only) | time | µmol/m²/s | Ecosystem respiration: CO₂ breathed out by plants and soil |
| LWup | Radiation (obs only) | time | W/m² | Heat radiation emitted up by the surface |
| SWup | Radiation (obs only) | time | W/m² | Sunlight reflected back up |
| GPP_LL | Carbon (obs only) | time | µmol/m²/s | GPP from the Lasslop light-response method (= −NEE + respiration) |
| GPP_LL_all | Carbon (obs only) | time | µmol/m²/s | GPP fully modelled by the Lasslop (2010) method |
| GPP_SOLO | Carbon (obs only) | time | µmol/m²/s | GPP estimated with the SOLO neural-network method |
| Ustar | Turbulence (obs only) | time | m/s | Friction velocity: how turbulent the air is; low values make fluxes unreliable |
| Qle_cor | Energy (obs only) | time | W/m² | Qle adjusted so the energy balance closes |
| Qh_cor | Energy (obs only) | time | W/m² | Qh adjusted so the energy balance closes |
| Tair_C | Added by the build | time | °C | Tair converted from K |
| VPD_kPa | Added by the build | time | kPa | VPD converted from hPa |
| Precip_mm | Added by the build | time | mm per step | Rain fallen in each half-hour |
| measured_Qle, _Qh, _Rnet, _Qg, _NEE, _GPP | Added by the build | time | True/False | True where that obs flux was really measured (qc = 0) |
| **Qle** | **Compared flux** | source, time | W/m² | Latent heat: energy used to evaporate water |
| **Qh** | **Compared flux** | source, time | W/m² | Sensible heat: energy that heats the air |
| **Rnet** | **Compared flux** | source, time | W/m² | Net radiation: total energy the surface keeps (≈ Qle + Qh + Qg) |
| **Qg** | **Compared flux** | source, time | W/m² | Ground heat flux: energy into the soil |
| **NEE** | **Compared flux** | source, time | µmol/m²/s | Net CO₂ flow; positive = released to the air. NaN for JULES |
| **GPP** | **Compared flux** | source, time | µmol/m²/s | CO₂ taken up by photosynthesis. Tower GPP is calculated from NEE, not measured directly |
| drought | Event flag | time | True/False | 90-day rain below the seasonal 10th percentile, runs of 14+ days |
| heatwave | Event flag | time | True/False | Bureau Excess Heat Factor > 0 (3-day mean T above the site's p95) |
| heatwave_ehf | Event flag | time | °C² | Heatwave severity (Excess Heat Factor); 0 outside heatwaves |
| high_vpd | Event flag | time | True/False | Daytime-mean VPD above the 90th percentile |
| normal | Event flag | time | True/False | None of the three events |

Quality flags present: CO2air_qc, LWdown_qc, SWdown_qc, Precip_qc, RH_qc, Qair_qc, Tair_qc, VPD_qc, Wind_qc, Psurf_qc, Resp_qc, NEE_qc, Qle_qc, Qg_qc, Qh_qc, LWup_qc, Rnet_qc, SWup_qc, GPP_LL_qc, GPP_LL_all_qc, GPP_qc, GPP_SOLO_qc, Ustar_qc.

# Open questions

Domain or data uncertainties that a result depends on. Each says what was done
for now, so the choice can be revisited. Evidence lives in `outputs/build/`
(`alignment.csv`, `<SITE>_corrections.csv`) and `outputs/inventory/`.

## Model outputs

1. **NoahMP NEE units and sign.** The file says `g m-2 s-1`. First read as g
   CO2 (as carbon its midday uptake looked 6x the obs). PLUMBER2 (p. 5526,
   Figs 3-4) says Noah-MP over-predicts carbon uptake with a very high WUE bias;
   read as grams of carbon it does exactly that (more uptake than obs at 73% of
   AU sites, WUE error +191%), read as CO2 it does not (-21%).
   *Decision (2026-10-02):* grams of carbon. Sign still flipped (uptake-positive).
2. **CHTESSEL_ERA5_3 / _4 forcing.** Their net shortwave matches the tower
   sunlight only at r = 0.87-0.95, with lags that change sign between sites
   (CHTESSEL_1 matches exactly). They look forced by ERA5, not the tower.
   *Decision (2026-10-01):* excluded from the build (`role: excluded`).
3. **BEPS timing cannot be verified.** BEPS writes no forcing or radiation, so
   only the weak Qle-vs-sunlight check is possible (lag 1-3 steps vs 0-1 for
   the obs). *Decision (2026-10-01):* excluded from the build (`role: excluded`).
4. **EntTBM is 60 min late.** Its SWnet copies the obs sunlight exactly but one
   hour late at all four starter sites. *For now:* shifted by -60 min.
5. **ORCHIDEE GPP/NEE per plant type.** First summed with `vegetfrac` as
   weights, which made arid-site GPP 5-20x below the obs. Each type's GPP turns
   out to be roughly proportional to its vegetfrac already, so the values look
   pre-scaled by cover. *Now:* plain sum over plant types. INFERRED; confirm
   with ORCHIDEE documentation.
6. **ORCHIDEE3 Qg sign** is opposite to ORCHIDEE2's. *For now:* flipped.
7. **MuSICA NEE** barely tracks the obs (r = -0.07..0.35, negative at AU-DaS,
   AU-Dry, AU-TTE) and ranges over -39..37 umol/m2/s at AU-Cum, while its GPP
   looks normal. *Decision (2026-10-01):* dropped (`drop_vars: [NEE]`).
8. **Sign fixes** (correlation with obs was negative before): CHTESSEL GPP,
   GFDL NEE, NoahMP NEE, STEMMUS-SCOPE NEE, ORCHIDEE3 Qg. Plausible, but
   worth confirming against each model's documentation.
9. **Which runs are `CABLE` and `JULES`?** CABLE or CABLE-POP; JULES_GL9 or
   JULES_GL9_LAI? And how does each get its LAI (needed to pick a fair pair)?
10. **GFDL AU-Ync** starts 11 h after the obs, and its SWdown matches the obs at
    no lag within +-15 h (best r = 0.41): not a simple shift. Excluded at that site.
13. **QUINCY AU-Wrr Qh** is broken (mean 190 W/m2 vs 5 in the obs, peaks at
    22:00, range -627..1518). QUINCY excluded at that site.
14. **STEMMUS-SCOPE NEE at AU-TTE** correlates at -0.21 with the obs (positive
    at every other site), so not a sign error. A very dry site with small NEE;
    probably a genuine model failure. Kept.
15. **MATSIRO AU-Ync** has no value for 2.4% of steps. Kept (gap is NaN).
11. **Missing sites:** ORCHIDEE2/3 have no AU-Rob or AU-Wrr; QUINCY has no AU-Ync.

## Observations

12. **VPD.** The Met file's VPD is labelled hPa but holds kPa at most sites. At
    AU-ASM, AU-Cow, AU-TTE, AU-Ync and AU-Rob it does not match the other
    humidity variables at all, and RH disagrees with Qair there too.
    *For now:* `VPD_kPa` is recomputed from Tair, Qair and Psurf (what the
    models are forced with). The file value is kept as `VPD_file`.

16. **AU-Ync tower files are Alice Springs data.** `AU-Ync_2011-2017_OzFlux_Met.nc`
    and `_Flux.nc` hold exactly the same SWdown, Tair, Precip, VPD, LWdown and Qle
    as AU-ASM (array-identical); only lat/lon differ, and `site_name` in both is
    "AliceSpringsMulga". Most models' AU-Ync runs were forced with this ASM data
    (CABLE, NoahMP, STEMMUS-SCOPE, MuSICA SWdown r = 1.0 with ASM), while GFDL's
    was not (r = 0.02), which is why GFDL "failed" there. AU-Ync is therefore a
    duplicate of AU-ASM with Yanco coordinates. *Recommendation:* exclude AU-Ync
    from all analysis and tell Anna/Gab (PLUMBER2 v1.0 dataset issue). The
    GFDL AU-Ync exclusion becomes moot.

import os

SITES = [
    "AU-ASM",
    "AU-Cow",
    "AU-Cpr",
    "AU-Ctr",
    "AU-Cum",
    "AU-DaP",
    "AU-DaS",
    "AU-Dry",
    "AU-Emr",
    "AU-Gin",
    "AU-GWW",
    "AU-How",
    "AU-Lit",
    "AU-Otw",
    "AU-Rig",
    "AU-Rob",
    "AU-Sam",
    "AU-Stp",
    "AU-TTE",
    "AU-Tum",
    "AU-Whr",
    "AU-Wrr",
    "AU-Ync"
]

MODELS = {
    "BEPS": "BEPS_{}.nc",
    "CABLE": "{}_out.nc",
    "CLM5": "{}-CLM5.nc",
    "EntTBM": "{}_NASAEnt_Phase1.nc",
    "GFDL": "{}.nc",
    "JULES": "local_{0}_fluxnet2015_GAL9.{0}.nc",
    "LPJ-GUESS": "LPJ-GUESS_{}_plumber.nc",
    "NoahMP": "NoahMPv401_UAlb_r1a_{}.nc",
    "QUINCY": "QUINCY_{}.nc",
    # Empirical benchmarks and simple physical baselines (added 2026-10-01)
    "1lin": "Emp1lin_{}.nc",
    "3km27": "3km27_{}_robust_output.nc",
    "6km729": "6km729_{}_robust_output.nc",
    "6km729lag": "6km729shortlag_{}_robust_output.nc",
    "LSTM": "{}_lstm_output.nc",
    "RF": "{}_RF_output.nc",
    "Manabe": "local_{0}_fluxnet2015_Manabe.{0}.nc",
    "PenmanMonteith": "local_{0}_fluxnet2015_faopm.{0}.nc"
}

for model_name, file_path in MODELS.items():
    for site in SITES:
        old_name = f"data/models/{model_name}/{file_path.format(site)}"
        new_name = f"data/models/{model_name}/{site}.nc"
        try:
            os.rename(old_name, new_name)
            print("File renamed successfully!")
        except FileNotFoundError:
            print("The source file does not exist.")
        except FileExistsError:
            print("A file with the new name already exists.")
            
CHTESSEL_1 = [
    "chtessel_exp1_insitu_AU-ASM_2011-2017.nc",
    "chtessel_exp1_insitu_AU-Cow_2010-2015.nc",
    "chtessel_exp1_insitu_AU-Cpr_2011-2017.nc",
    "chtessel_exp1_insitu_AU-Ctr_2010-2017.nc",
    "chtessel_exp1_insitu_AU-Cum_2013-2018.nc",
    "chtessel_exp1_insitu_AU-DaP_2009-2012.nc",
    "chtessel_exp1_insitu_AU-DaS_2010-2017.nc",
    "chtessel_exp1_insitu_AU-Dry_2011-2015.nc",
    "chtessel_exp1_insitu_AU-Emr_2012-2013.nc",
    "chtessel_exp1_insitu_AU-Gin_2012-2017.nc",
    "chtessel_exp1_insitu_AU-GWW_2013-2017.nc",
    "chtessel_exp1_insitu_AU-How_2003-2017.nc",
    "chtessel_exp1_insitu_AU-Lit_2016-2017.nc",
    "chtessel_exp1_insitu_AU-Otw_2009-2010.nc",
    "chtessel_exp1_insitu_AU-Rig_2011-2016.nc",
    "chtessel_exp1_insitu_AU-Rob_2014-2017.nc",
    "chtessel_exp1_insitu_AU-Sam_2011-2017.nc",
    "chtessel_exp1_insitu_AU-Stp_2010-2017.nc",
    "chtessel_exp1_insitu_AU-TTE_2013-2017.nc",
    "chtessel_exp1_insitu_AU-Tum_2002-2017.nc",
    "chtessel_exp1_insitu_AU-Whr_2015-2016.nc",
    "chtessel_exp1_insitu_AU-Wrr_2016-2017.nc",
    "chtessel_exp1_insitu_AU-Ync_2011-2017.nc"
]

CHTESSEL_ERA5_3 = [
    "chtessel_exp3_ol_era5_AU-ASM_2011-2017.nc",
    "chtessel_exp3_ol_era5_AU-Cow_2010-2015.nc",
    "chtessel_exp3_ol_era5_AU-Cpr_2011-2017.nc",
    "chtessel_exp3_ol_era5_AU-Ctr_2010-2017.nc",
    "chtessel_exp3_ol_era5_AU-Cum_2013-2018.nc",
    "chtessel_exp3_ol_era5_AU-DaP_2009-2012.nc",
    "chtessel_exp3_ol_era5_AU-DaS_2010-2017.nc",
    "chtessel_exp3_ol_era5_AU-Dry_2011-2015.nc",
    "chtessel_exp3_ol_era5_AU-Emr_2012-2013.nc",
    "chtessel_exp3_ol_era5_AU-Gin_2012-2017.nc",
    "chtessel_exp3_ol_era5_AU-GWW_2013-2017.nc",
    "chtessel_exp3_ol_era5_AU-How_2003-2017.nc",
    "chtessel_exp3_ol_era5_AU-Lit_2016-2017.nc",
    "chtessel_exp3_ol_era5_AU-Otw_2009-2010.nc",
    "chtessel_exp3_ol_era5_AU-Rig_2011-2016.nc",
    "chtessel_exp3_ol_era5_AU-Rob_2014-2017.nc",
    "chtessel_exp3_ol_era5_AU-Sam_2011-2017.nc",
    "chtessel_exp3_ol_era5_AU-Stp_2010-2017.nc",
    "chtessel_exp3_ol_era5_AU-TTE_2013-2017.nc",
    "chtessel_exp3_ol_era5_AU-Tum_2002-2017.nc",
    "chtessel_exp3_ol_era5_AU-Whr_2015-2016.nc",
    "chtessel_exp3_ol_era5_AU-Wrr_2016-2017.nc",
    "chtessel_exp3_ol_era5_AU-Ync_2011-2017.nc"
]

CHTESSEL_ERA5_4 = [
    "chtessel_exp4_ol_era5_AU-ASM_2011-2017.nc",
    "chtessel_exp4_ol_era5_AU-Cow_2010-2015.nc",
    "chtessel_exp4_ol_era5_AU-Cpr_2011-2017.nc",
    "chtessel_exp4_ol_era5_AU-Ctr_2010-2017.nc",
    "chtessel_exp4_ol_era5_AU-Cum_2013-2018.nc",
    "chtessel_exp4_ol_era5_AU-DaP_2009-2012.nc",
    "chtessel_exp4_ol_era5_AU-DaS_2010-2017.nc",
    "chtessel_exp4_ol_era5_AU-Dry_2011-2015.nc",
    "chtessel_exp4_ol_era5_AU-Emr_2012-2013.nc",
    "chtessel_exp4_ol_era5_AU-Gin_2012-2017.nc",
    "chtessel_exp4_ol_era5_AU-GWW_2013-2017.nc",
    "chtessel_exp4_ol_era5_AU-How_2003-2017.nc",
    "chtessel_exp4_ol_era5_AU-Lit_2016-2017.nc",
    "chtessel_exp4_ol_era5_AU-Otw_2009-2010.nc",
    "chtessel_exp4_ol_era5_AU-Rig_2011-2016.nc",
    "chtessel_exp4_ol_era5_AU-Rob_2014-2017.nc",
    "chtessel_exp4_ol_era5_AU-Sam_2011-2017.nc",
    "chtessel_exp4_ol_era5_AU-Stp_2010-2017.nc",
    "chtessel_exp4_ol_era5_AU-TTE_2013-2017.nc",
    "chtessel_exp4_ol_era5_AU-Tum_2002-2017.nc",
    "chtessel_exp4_ol_era5_AU-Whr_2015-2016.nc",
    "chtessel_exp4_ol_era5_AU-Wrr_2016-2017.nc",
    "chtessel_exp4_ol_era5_AU-Ync_2011-2017.nc"
]

MATSIRO = [
    "AU-ASM_2011-2017_OzFlux_out.nc",
    "AU-Cow_2010-2015_OzFlux_out.nc",
    "AU-Cpr_2011-2017_OzFlux_out.nc",
    "AU-Ctr_2010-2017_OzFlux_out.nc",
    "AU-Cum_2013-2018_OzFlux_out.nc",
    "AU-DaP_2009-2012_OzFlux_out.nc",
    "AU-DaS_2010-2017_OzFlux_out.nc",
    "AU-Dry_2011-2015_OzFlux_out.nc",
    "AU-Emr_2012-2013_OzFlux_out.nc",
    "AU-Gin_2012-2017_OzFlux_out.nc",
    "AU-GWW_2013-2017_OzFlux_out.nc",
    "AU-How_2003-2017_OzFlux_out.nc",
    "AU-Lit_2016-2017_OzFlux_out.nc",
    "AU-Otw_2009-2010_OzFlux_out.nc",
    "AU-Rig_2011-2016_OzFlux_out.nc",
    "AU-Rob_2014-2017_OzFlux_out.nc",
    "AU-Sam_2011-2017_OzFlux_out.nc",
    "AU-Stp_2010-2017_OzFlux_out.nc",
    "AU-TTE_2013-2017_OzFlux_out.nc",
    "AU-Tum_2002-2017_OzFlux_out.nc",
    "AU-Whr_2015-2016_OzFlux_out.nc",
    "AU-Wrr_2016-2017_OzFlux_out.nc",
    "AU-Ync_2011-2017_OzFlux_out.nc"
]

MuSICA = [
    "AU-ASM_2011-2017_OzFlux-musica.nc",
    "AU-Cow_2010-2015_OzFlux-musica.nc",
    "AU-Cpr_2011-2017_OzFlux-musica.nc",
    "AU-Ctr_2010-2017_OzFlux-musica.nc",
    "AU-Cum_2013-2018_OzFlux-musica.nc",
    "AU-DaP_2009-2012_OzFlux-musica.nc",
    "AU-DaS_2010-2017_OzFlux-musica.nc",
    "AU-Dry_2011-2015_OzFlux-musica.nc",
    "AU-Emr_2012-2013_OzFlux-musica.nc",
    "AU-Gin_2012-2017_OzFlux-musica.nc",
    "AU-GWW_2013-2017_OzFlux-musica.nc",
    "AU-How_2003-2017_OzFlux-musica.nc",
    "AU-Lit_2016-2017_OzFlux-musica.nc",
    "AU-Otw_2009-2010_OzFlux-musica.nc",
    "AU-Rig_2011-2016_OzFlux-musica.nc",
    "AU-Rob_2014-2017_OzFlux-musica.nc",
    "AU-Sam_2011-2017_OzFlux-musica.nc",
    "AU-Stp_2010-2017_OzFlux-musica.nc",
    "AU-TTE_2013-2017_OzFlux-musica.nc",
    "AU-Tum_2002-2017_OzFlux-musica.nc",
    "AU-Whr_2015-2016_OzFlux-musica.nc",
    "AU-Wrr_2016-2017_OzFlux-musica.nc",
    "AU-Ync_2011-2017_OzFlux-musica.nc"
]

ORCHIDEE2 = [
    "AU-ASMNewSchemeNewPFT.202111HIST_20110101_20171231_1M_sechiba_history.nc",
    "AU-CowNewSchemeNewPFT.202111HIST_20100101_20151231_1M_sechiba_history.nc",
    "AU-CprNewSchemeNewPFT.202111HIST_20110101_20171231_1M_sechiba_history.nc",
    "AU-CtrNewSchemeNewPFT.202111HIST_20100101_20171231_1M_sechiba_history.nc",
    "AU-CumNewSchemeNewPFT.202111HIST_20130101_20181231_1M_sechiba_history.nc",
    "AU-DaPNewSchemeNewPFT.202111HIST_20090101_20121231_1M_sechiba_history.nc",
    "AU-DaSNewSchemeNewPFT.202111HIST_20100101_20171231_1M_sechiba_history.nc",
    "AU-DryNewSchemeNewPFT.202111HIST_20110101_20151231_1M_sechiba_history.nc",
    "AU-EmrNewSchemeNewPFT.202111HIST_20120101_20131231_1M_sechiba_history.nc",
    "AU-GinNewSchemeNewPFT.202111HIST_20120101_20171231_1M_sechiba_history.nc",
    "AU-GWWNewSchemeNewPFT.202111HIST_20130101_20171231_1M_sechiba_history.nc",
    "AU-HowNewSchemeNewPFT.202111HIST_20030101_20171231_1M_sechiba_history.nc",
    "AU-LitNewSchemeNewPFT.202111HIST_20160101_20171231_1M_sechiba_history.nc",
    "AU-OtwNewSchemeNewPFT.202111.hourlyHIST_20090101_20101231_1M_sechiba_history.nc",
    "AU-RigNewSchemeNewPFT.202111HIST_20110101_20161231_1M_sechiba_history.nc",
    "AU-SamNewSchemeNewPFT.202111HIST_20110101_20171231_1M_sechiba_history.nc",
    "AU-StpNewSchemeNewPFT.202111HIST_20100101_20171231_1M_sechiba_history.nc",
    "AU-TTENewSchemeNewPFT.202111HIST_20130101_20171231_1M_sechiba_history.nc",
    "AU-TumNewSchemeNewPFT.202111.hourlyHIST_20020101_20171231_1M_sechiba_history.nc",
    "AU-WhrNewSchemeNewPFT.202111HIST_20150101_20161231_1M_sechiba_history.nc",
    "AU-YncNewSchemeNewPFT.202111HIST_20110101_20171231_1M_sechiba_history.nc"
]

ORCHIDEE3 = [
    "AU-ASMplumber.202307HIST_20110101_20171231_1M_sechiba_history.nc",
    "AU-Cowplumber.202307HIST_20100101_20151231_1M_sechiba_history.nc",
    "AU-Cprplumber.202307HIST_20110101_20171231_1M_sechiba_history.nc",
    "AU-Ctrplumber.202307HIST_20100101_20171231_1M_sechiba_history.nc",
    "AU-Cumplumber.202307HIST_20130101_20181231_1M_sechiba_history.nc",
    "AU-DaPplumber.202307HIST_20090101_20121231_1M_sechiba_history.nc",
    "AU-DaSplumber.202307HIST_20100101_20171231_1M_sechiba_history.nc",
    "AU-Dryplumber.202307HIST_20110101_20151231_1M_sechiba_history.nc",
    "AU-Emrplumber.202307HIST_20120101_20131231_1M_sechiba_history.nc",
    "AU-Ginplumber.202307HIST_20120101_20171231_1M_sechiba_history.nc",
    "AU-GWWplumber.202307HIST_20130101_20171231_1M_sechiba_history.nc",
    "AU-Howplumber.202307HIST_20030101_20171231_1M_sechiba_history.nc",
    "AU-Litplumber.202307HIST_20160101_20171231_1M_sechiba_history.nc",
    "AU-Otwplumber1h.202307HIST_20090101_20101231_1M_sechiba_history.nc",
    "AU-Rigplumber.202307HIST_20110101_20161231_1M_sechiba_history.nc",
    "AU-Samplumber.202307HIST_20110101_20171231_1M_sechiba_history.nc",
    "AU-Stpplumber.202307HIST_20100101_20171231_1M_sechiba_history.nc",
    "AU-TTEplumber.202307HIST_20130101_20171231_1M_sechiba_history.nc",
    "AU-Tumplumber1h.202307HIST_20020101_20171231_1M_sechiba_history.nc",
    "AU-Whrplumber.202307HIST_20150101_20161231_1M_sechiba_history.nc",
    "AU-Yncplumber.202307HIST_20110101_20171231_1M_sechiba_history.nc"
]

SDGVM = [
    "SDGVM_PLUMMER2_AU-ASM_2011_2017.nc",
    "SDGVM_PLUMMER2_AU-Cow_2010_2015.nc",
    "SDGVM_PLUMMER2_AU-Cpr_2011_2017.nc",
    "SDGVM_PLUMMER2_AU-Ctr_2010_2017.nc",
    "SDGVM_PLUMMER2_AU-Cum_2013_2018.nc",
    "SDGVM_PLUMMER2_AU-DaP_2009_2012.nc",
    "SDGVM_PLUMMER2_AU-DaS_2010_2017.nc",
    "SDGVM_PLUMMER2_AU-Dry_2011_2015.nc",
    "SDGVM_PLUMMER2_AU-Emr_2012_2013.nc",
    "SDGVM_PLUMMER2_AU-Gin_2012_2017.nc",
    "SDGVM_PLUMMER2_AU-GWW_2013_2017.nc",
    "SDGVM_PLUMMER2_AU-How_2003_2017.nc",
    "SDGVM_PLUMMER2_AU-Lit_2016_2017.nc",
    "SDGVM_PLUMMER2_AU-Otw_2009_2010.nc",
    "SDGVM_PLUMMER2_AU-Rig_2011_2016.nc",
    "SDGVM_PLUMMER2_AU-Rob_2014_2017.nc",
    "SDGVM_PLUMMER2_AU-Sam_2011_2017.nc",
    "SDGVM_PLUMMER2_AU-Stp_2010_2017.nc",
    "SDGVM_PLUMMER2_AU-TTE_2013_2017.nc",
    "SDGVM_PLUMMER2_AU-Tum_2002_2017.nc",
    "SDGVM_PLUMMER2_AU-Whr_2015_2016.nc",
    "SDGVM_PLUMMER2_AU-Wrr_2016_2017.nc",
    "SDGVM_PLUMMER2_AU-Ync_2011_2017.nc"
]

STEMMUS_SCOPE = [
    "AU-ASM_2022-07-28-1043_STEMMUS_SCOPE.nc",
    "AU-Cow_2022-07-28-1043_STEMMUS_SCOPE.nc",
    "AU-Cpr_2022-07-28-1043_STEMMUS_SCOPE.nc",
    "AU-Ctr_2022-07-28-1044_STEMMUS_SCOPE.nc",
    "AU-Cum_2022-07-15-1101_STEMMUS_SCOPE.nc",
    "AU-DaP_2022-07-15-1100_STEMMUS_SCOPE.nc",
    "AU-DaS_2022-07-28-1044_STEMMUS_SCOPE.nc",
    "AU-Dry_2022-07-28-1043_STEMMUS_SCOPE.nc",
    "AU-Emr_2022-07-28-1043_STEMMUS_SCOPE.nc",
    "AU-Gin_2022-07-28-1043_STEMMUS_SCOPE.nc",
    "AU-GWW_2022-07-28-1043_STEMMUS_SCOPE.nc",
    "AU-How_2022-07-28-1045_STEMMUS_SCOPE.nc",
    "AU-Lit_2022-07-15-1059_STEMMUS_SCOPE.nc",
    "AU-Otw_2022-07-15-1100_STEMMUS_SCOPE.nc",
    "AU-Rig_2022-07-28-1043_STEMMUS_SCOPE.nc",
    "AU-Rob_2023-03-21-1452_STEMMUS_SCOPE.nc",
    "AU-Sam_2022-07-28-1044_STEMMUS_SCOPE.nc",
    "AU-Stp_2022-07-28-1044_STEMMUS_SCOPE.nc",
    "AU-TTE_2022-07-28-1043_STEMMUS_SCOPE.nc",
    "AU-Tum_2022-07-28-1044_STEMMUS_SCOPE.nc",
    "AU-Whr_2022-07-28-1043_STEMMUS_SCOPE.nc",
    "AU-Wrr_2023-03-24-0715_STEMMUS_SCOPE.nc",
    "AU-Ync_2022-07-15-1101_STEMMUS_SCOPE.nc"
]

ALL_ARRAYS = {
    "CHTESSEL_1": CHTESSEL_1,
    "CHTESSEL_ERA5_3": CHTESSEL_ERA5_3,
    "CHTESSEL_ERA5_4": CHTESSEL_ERA5_4,
    "MATSIRO": MATSIRO,
    "MuSICA": MuSICA,
    "ORCHIDEE2": ORCHIDEE2,
    "ORCHIDEE3": ORCHIDEE3,
    "SDGVM": SDGVM,
    "STEMMUS-SCOPE": STEMMUS_SCOPE,
}

for model, files in ALL_ARRAYS.items():
    for file_path in files:
        for site in SITES:
            if site in file_path:
                old_name = f"data/models/{model}/{file_path}"
                new_name = f"data/models/{model}/{site}.nc"
                try:
                    os.rename(old_name, new_name)
                    print("File renamed successfully!")
                except FileNotFoundError:
                    print("The source file does not exist.")
                except FileExistsError:
                    print("A file with the new name already exists.")
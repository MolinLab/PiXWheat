# 01_Preparation

These scripts harmonise German and Austrian winter-wheat phenology observations and extract daily weather and ERA5-Land variables. Run them from the repository root, or set `PHENOLOGY_PROJECT_DIR` to the absolute repository path. Scripts live under `Scripts/01_Preparation/`; generated and supplied data live under `Data/`.

## Directory layout

```text
Data/
├── DWD/          # DWD phenology outputs and HYRAS downloads
├── PEP725/       # PEP725 input and Austrian phenology QC outputs
├── GeoSphere/    # Austrian AGES/WheatVIZ input and SPARTACUS data
├── Era5/         # Annual ERA5-Land NetCDF and site-level extraction
├── BRIWECS/      # BRIWECS site-coordinate input
├── SoilGrids/     # SoilGrids export
└── Combine_DWD_AT.R # Create Machine learning Dataset
```
## Scripts and order

| Step | Script | Main input and output |
| --- | --- | --- |
| German phenology | `DWD/DWD_phenology_extraction.R` | Downloads DWD annual-reporter observations; writes `Data/DWD/DWD_phenology_16_06_final.csv` and QC tables. |
| German weather | `DWD/get_dwd_data.sh` | Downloads annual HYRAS NetCDF files to `Data/DWD/{tasmax,tasmin,pr,rsds}/`. |
| German site weather | `DWD/extract_weather_dwd.R` | Extracts HYRAS values at DWD sites and writes annual CSVs beneath the HYRAS variable directories. |
| Austrian phenology | `GeoSphere/phenological_qc_at.R` | Reads `Data/PEP725/pep725_with_geo_labels.csv` and `Data/GeoSphere/WVIZ_AGES_pheno.csv`; writes `Data/PEP725/pep_at_filtered_final_19_06.csv` and QC tables. |
| Austrian weather | `GeoSphere/get_geosphere_data.sh` | Downloads SPARTACUS NetCDF files to `Data/GeoSphere/{TN,TX,RR,SA}/`. |
| Austrian weather extraction | `GeoSphere/extract_weather_geosphere.R` | Extracts daily SPARTACUS values at Austrian observation locations. |
| ERA5-Land annual data | `Era5/extract_radiation_era5_with_soilwater_new.py` | Builds `Data/Era5/era5land_daily_yearly/era5land_daily_<year>.nc`. |
| ERA5-Land site extraction | `DWD/extract_and_buildphenology_ger_with_soil.R`, `GeoSphere/extract_and_buildphenology_geosphere_with_soil.R`, `BRIWECS/extract_BRIWECS.R` | Extracts the annual ERA5-Land NetCDF variables at German, Austrian and BRIWECS locations, respectively. |
| Soil hydraulic limits | `SoilGrids/prepare_soilgrids.R` | Reads `Data/SoilGrids/SoilGrids_Final_Fixed_NoNulls_17_06.csv`; writes `Data/Soil_hydraulic_ERA5_0_100.csv` and `Data/Soil_data_16_07.csv`. |
| German process predictors | `DWD/calculate_predictors_germany.R` | Combines phenology, extracted weather and soil hydraulic limits; writes `Data/DWD/DWD_phenology_GDD_era5_test1.csv`. |
| Austrian process predictors | `GeoSphere/calculate_predictors_austria.R` | Combines observation dates, extracted weather and the soil hydraulic table; writes `Data/GeoSphere/at_phenology_GDD_era5_test.csv`. |
| BRIWECS HYRAS weather | `BRIWECS/extract_briwecs_hyras_2014_2022.R` | Extracts HYRAS temperature and precipitation for six supplied sites in 2014–2022; writes `Data/BRIWECS/DWD_HYRAS_2014_2022_all_locations.csv`. |
| Dataset | `Combine_DWD_AT.R` | Combines data and creates final Dataset; writes `pheno_df_CV.csv`.|
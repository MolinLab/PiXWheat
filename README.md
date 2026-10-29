# PiXWheat
**P**rocess-**i**nformed predictors and **X**AI for **Wheat** development

# Explainable AI for winter wheat phenology

Code accompanying the manuscript **“Explainable Artificial Intelligence for Predicting Wheat Phenological Development Across Multi-Environment Agricultural Datasets”** by Joel J. Milek, Luca G. Bernardini, Lukas J. Koppensteiner, Andreas Holzinger and Eva M. Molin.

## Study overview

The study predicts the BBCH developmental stage of winter wheat on a supplied date using environmental information available up to that date. The German and Austrian development dataset comprises 251,828 observations from 45,631 site–season–sowing environments, 3,084 physical sites and 45 harvest years (1981–2025). Independent BRIWECS trials assess transfer to a different observation source and the accuracy of predicted milestone dates.

The analyses compare developmental clocks, environmental predictor representations and four learners (GAM, MARS, random forest and XGBoost). The selected model is evaluated under five transfer schemes and examined with out-of-fold TreeSHAP. A targeted GAM comparison tests whether explicitly modelling the thermal time by radiation interaction improves prediction. The manuscript and supplementary methods define the study cohort, predictors, validation schemes and reported estimands.

The prediction task uses **observed weather through each prediction date**. Expanding-origin evaluation restricts model training to earlier seasons; it does not assess a weather-forecast-driven operational forecast. Predicting BBCH on a supplied date and predicting the first date on which a milestone is reached are separate evaluation tasks.

## Repository layout

```text
.
├── README.md
├── Scripts/
│   ├── 01_Preparation/
│   │   ├── DWD/
│   │   ├── GeoSphere/
│   │   ├── Era5/
│   │   ├── SoilGrids/
│   │   ├── BRIWECS/
│   │   └── combine_DWD_AT.R
│   ├── 02_Predictor_Representation_Learner_Comparison/
│   ├── 03_Transfer_Selected_Model/
│   ├── 04_SHAP_Analysis/
│   └── 05_External_Validation/
├── Data/                         # supplied and generated data
├── benchmark_results/            # generated results from folders 02 and 03
├── shap_results/                 # generated results from folder 04
└── external_validation/          # generated results from folder 05
```

Each numbered script folder has its own README with its inputs, execution options and outputs. `Data/` and the three result directories represent the working layout; their contents depend on which data can be redistributed and which analyses have been run. The source-specific preparation scripts are in the `DWD`, `GeoSphere`, `Era5`, `SoilGrids` and `BRIWECS` subfolders of `01_Preparation`.

## Data and analysis workflow

Development observations come from the German Weather Service (DWD), PEP725, AGES and WheatVIZ. Weather inputs include DWD HYRAS, GeoSphere Austria SPARTACUS and ERA5-Land; soil information includes SoilGrids-derived inputs. BRIWECS trial records are used only for the independent assessment. Consult the manuscript and the folder READMEs for source versions, filters, preparation steps and attribution.

| Folder | Main scripts | Purpose |
| --- | --- | --- |
| `01_Preparation` | Source-specific scripts and `combine_DWD_AT.R` | Harmonise observations, extract weather and soil inputs, calculate predictors and create `Data/pheno_df_CV.csv`. |
| `02_Predictor_Representation_Learner_Comparison` | `Predictor_Representation_benchmark.R`, `Learner_selection_benchmark.R` | Compare developmental clocks, predictor representations and learners. |
| `03_Transfer_Selected_Model` | `multischeme_generalisation_benchmark.R` | Compare the selected learner with a calibrated thermal-time baseline across five transfer schemes. |
| `04_SHAP_Analysis` | `SHAP_XGBoost.R` | Reconstruct held-out XGBoost predictions and calculate TreeSHAP attributions and interactions. |
| `05_External_Validation` | `calculate_briwecs_environment_targets_paper_style.R`, `Validate_XGBoost_Briwecs.R` | Construct independent environment-level targets and evaluate stage and milestone-date predictions; run the targeted GAM interaction comparison. |

The main representation contains cumulative weather predictors from sowing to the observation date. A second representation combines process-informed variables with recent weather summaries. These predictors incorporate crop-modelling response functions, but they are not outputs of a fully calibrated crop simulation model. The transfer script supplied selects the cumulative representation; the external-validation script can evaluate cumulative or hybrid predictors when the corresponding transfer configuration and model cache are available.

## Running the analyses

Run the scripts from the repository root with R and the packages specified in each script. Set `PHENOLOGY_PROJECT_DIR` to an absolute repository path when running from elsewhere. Prepare or place the required source data under `Data/` as described in `Scripts/01_Preparation/README.md`, then produce `Data/pheno_df_CV.csv` with `combine_DWD_AT.R`. The analytical stages can then be invoked as follows:

```bash
Rscript Scripts/02_Predictor_Representation_Learner_Comparison/Predictor_Representation_benchmark.R
Rscript Scripts/02_Predictor_Representation_Learner_Comparison/Learner_selection_benchmark.R
Rscript Scripts/03_Transfer_Selected_Model/multischeme_generalisation_benchmark.R

PHENOLOGY_GENERALISATION_RUN_ID=multischeme_xgboost_cumulative_weather_paper \
  Rscript Scripts/04_SHAP_Analysis/SHAP_XGBoost.R

Rscript Scripts/05_External_Validation/calculate_briwecs_environment_targets_paper_style.R
PHENOLOGY_GENERALISATION_RUN_ID=multischeme_xgboost_cumulative_weather_paper \
PHENOLOGY_BRIWECS_FEATURE_SET=cumulative \
  Rscript Scripts/05_External_Validation/Validate_XGBoost_Briwecs.R
```

The transfer benchmark writes to `benchmark_results/<run ID>/`; SHAP and external validation must reference the **same saved benchmark run** and its development input. The SHAP script's original default run ID differs from the transfer script's current default, so the example sets it explicitly. The SHAP script calculates only spatial-block and expanding-origin attributions by default; set `PHENOLOGY_SHAP_SCHEMES=spatial_block,leave_one_year_out,expanding_origin,environmental_zone,joint_space_time` to request all five schemes.

External validation requires its BRIWECS inputs and a compatible development-only model cache; if that cache is unavailable, set `PHENOLOGY_BRIWECS_ALLOW_MODEL_REFIT=TRUE` to refit from development data. 
Change run IDs when inputs, predictors or model specifications change so checkpoints are not reused across different experiments.

## Validation and interpretation

Spatial blocks, environmental-zone exclusion, leave-one-year-out, expanding origin and joint space–time exclusion assess distinct transfer settings. In particular, leave-one-year-out can use both earlier and later years for training, whereas expanding origin uses only earlier years.

Overall MAE and RMSE weight observations equally. Environment–observed-stage macro metrics weight each environment–stage unit equally; the macro-RMSE used in the analysis is `sqrt(mean(environment_stage_MSE))`. Independent BRIWECS assessment additionally reports errors for the dates at which BBCH milestones are first predicted to occur. TreeSHAP values describe attribution within the fitted model, not causal effects of individual environmental variables.

## Reproducibility, data access and citation

The scripts save configurations, fold definitions, out-of-fold predictions, performance tables and session information in run-specific result folders. Use the matching cohort file, benchmark run ID and package versions when reconstructing a manuscript result. The targeted GAM comparison is implemented in `Validate_XGBoost_Briwecs.R`; its external evaluation uses the independent BRIWECS targets.

Citation: #########

**Manuscript:** Milek, J. J., Bernardini, L. G., Koppensteiner, L. J., Holzinger, A. & Molin, E. M. *Explainable Artificial Intelligence for Predicting Wheat Phenological Development Across Multi-Environment Agricultural Datasets*.

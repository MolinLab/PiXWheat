# PiXWheat
Explainable AI for Predicting Wheat Development Using Process-Informed Predictors

Code accompanying the manuscript **“Explainable Artificial Intelligence for Predicting Wheat Phenological Development Across Multi-Environment Agricultural Datasets”** by Joel Milek, Luca Bernardini, Lukas J. Koppensteiner, Andreas Holzinger and Eva M. Molin.

## Study overview

We predict the developmental stage of winter wheat on an **observation date**, expressed on the BBCH scale, from environmental information available up to that date. Phenological records from Germany and Austria are combined with weather and soil data. The analysis compares developmental-clock formulations, predictor representations and four statistical learners; evaluates transfer under five validation schemes; and examines the selected XGBoost model using out-of-fold TreeSHAP. Independent BRIWECS trials assess external stage and milestone-date accuracy.

The modelling dataset contains 251,828 observations from 45,631 site–season–sowing environments, 3,084 physical sites and 45 harvest years (1981–2025). The manuscript and its supplementary methods give the definitive inclusion rules, predictor definitions, validation design and results.

**Prediction task:** Given the date and environmental history up to that date, estimate BBCH stage. An accurate stage estimate on a supplied date does not necessarily imply an accurate estimate of the date when a BBCH milestone is first reached. The expanding-origin experiment restricts training to preceding seasons but uses *observed* weather through the prediction date; it is not a weather-forecast-driven operational forecast.

## Data and predictors

Phenological observations originate from the German Weather Service (DWD) annual crop reporter network, PEP725, Austrian AGES and WheatVIZ data. BRIWECS trials are used for independent assessment. Weather inputs include DWD HYRAS, GeoSphere Austria SPARTACUS and ERA5-Land; soil information includes SoilGrids and derived soil-water variables. 

The main comparison uses two environmental representations:

- **Cumulative weather:** variables accumulated from sowing to the observation date, including beta-response thermal time, radiation and precipitation.
- **Hybrid process and recent weather:** process-informed development variables together with recent weather summaries.

The process-informed variables adapt published temperature, photoperiod, vernalisation and water-response functions. They are predictors, **not outputs of a fully calibrated crop simulation model**. No cultivar identity is included as an explicit model input.

## Analysis scripts

The filenames below are those used in the study workflow. Place the scripts in `scripts/`, or adjust the example commands to match the final repository layout.

| Stage | Script | Role |
| --- | --- | --- |
| Phenology quality control | `phenological_qc_at.R` and source-specific preparation scripts | Harmonise and check source observations. |
| Daily predictors | `calc_cumGDD_ger_point3_corrected.R`, `calc_cumGDD_austria_point3_corrected.R` | Build weather and process-informed predictors for Germany and Austria. |
| Combined modelling data | `combine_DWD_AT.R` | Combine source-specific prepared tables. |
| Development-clock and feature comparison | `benchmark_features_point4_hypotheses.R` | Compare clock formulations and predictor representations under spatial and temporal blocks. |
| Learner selection | `benchmark_learner_selection.R` | Compare GAM, MARS, random forest and XGBoost on the selected representations, with training-fold tuning as applicable. |
| Transfer assessment | `multischeme_generalisation_benchmark.R` | Evaluate the frozen cumulative representation and selected learner against a calibrated thermal-time baseline under five transfer schemes. |
| Model explanation | `shap_xgboost_phenology.R` | Reconstruct fold models, verify held-out predictions and calculate TreeSHAP values and interactions. |
| Independent trials | `validate_xgboost_briwecs.R` | Refit using development data and evaluate BRIWECS stage and milestone-date predictions. |
| Figures | `plot_point4_and_learner_selection.R`, `plot_generalisation_benchmark.R` | Render comparison and transfer figures from saved results. |

Some preprocessing scripts have source-specific paths and settings. Read each script's **USER CONFIGURATION** block before running it. 

## Running the core analyses

The modelling scripts are R scripts intended to run from the command line. They require the prepared input tables and the packages listed in their `required_packages` declarations. Principal packages include `data.table`, `mlr3`, `mlr3learners`, `mlr3tuning`, `mlr3mbo`, `mgcv`, `earth`, `ranger`, `xgboost`, `sf` and `blockCV`. 

On Linux or macOS, set `PHENOLOGY_PROJECT_DIR` to the directory containing the prepared data and run the stages in order:

```bash
export PHENOLOGY_PROJECT_DIR=/path/to/phenology-project

Rscript scripts/benchmark_features_point4_hypotheses.R
Rscript scripts/benchmark_learner_selection.R
Rscript scripts/multischeme_generalisation_benchmark.R
Rscript scripts/shap_xgboost_phenology.R
Rscript scripts/validate_xgboost_briwecs.R
```

To calculate TreeSHAP for all five transfer schemes rather than the script's default subset:

```bash
export PHENOLOGY_SHAP_SCHEMES=spatial_block,leave_one_year_out,expanding_origin,environmental_zone,joint_space_time
Rscript scripts/shap_xgboost_phenology.R
```

The explanation stage requires the transfer benchmark's saved fold manifest, model specifications and out-of-fold predictions. It checks that reconstructed fold predictions match the stored predictions before accepting SHAP values. Set run IDs or result directories consistently if you change the defaults.

## Validation and metrics

The five transfer schemes answer different questions: spatial blocks hold out new locations; leave-one-year-out holds out a season while allowing earlier and later seasons in training; expanding-origin uses only preceding seasons; environmental-zone validation holds out a zone; and joint space–time validation excludes both the test spatial block and temporal block from training. 

Performance is reported in BBCH units. Overall MAE and RMSE weight observations equally. Environment–observed-stage macro-MAE gives each environment–stage unit equal weight. 

See the supplementary methods for the exact estimands, clustering and uncertainty calculations. BRIWECS analyses additionally report milestone-date errors in days.

## Reproducibility and citation

The transfer script writes fold specifications, out-of-fold predictions, performance tables, a run configuration and `session_info.txt` to a run-specific directory under `benchmark_results/`.

**Citation:** Milek, J., Bernardini, L., Koppensteiner, L. J., Holzinger, A. & Molin, E. M. *Explainable Artificial Intelligence for Predicting Wheat Phenological Development Across Multi-Environment Agricultural Datasets*.

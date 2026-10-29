# ==============================================================================
# Learner-selection benchmark for winter-wheat phenology
# ==============================================================================
# Select among GAM, MARS, random forest and XGBoost after the predictor
# representation comparison. Cumulative weather is the primary representation;
# the process-plus-recent-weather set is a sensitivity analysis.
# Spatial and contiguous temporal outer folds share a cohort within this
# script. MARS, random forest and XGBoost use site-grouped inner tuning folds.
# The GAM estimates smoothness by REML within each outer training set.
# Select using environment-stage macro-MAE; predictions are not clipped.
# Run from the repository root or set PHENOLOGY_PROJECT_DIR to its path.

# ------------------------------------------------------------------------------
# 1. USER CONFIGURATION
# ------------------------------------------------------------------------------
PROJECT_DIR <- Sys.getenv("PHENOLOGY_PROJECT_DIR", unset = getwd())
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = TRUE)
INPUT_FILE <- Sys.getenv(
  "PHENOLOGY_LEARNER_INPUT",
  unset = file.path(PROJECT_DIR, "Data", "pheno_df_CV.csv")
)

# Change RUN_ID or CHECKPOINT_VERSION whenever predictors, formulas, folds,
# learners or tuning spaces change. Do not point this run at old BAM checkpoints.
RUN_ID <- Sys.getenv(
  "PHENOLOGY_LEARNER_BENCHMARK_RUN_ID",
  unset = "learner_selection_cumulative_hybrid_paper"
)
OUTPUT_DIR <- file.path(PROJECT_DIR, "benchmark_results", RUN_ID)
CHECKPOINT_DIR <- file.path(OUTPUT_DIR, "checkpoints")

SEED <- 123L
TARGET <- "BBCH"

# The script first tries SITE_GROUP_COLUMN and then the listed fallbacks.
# This must identify a physical site, not a site-year environment.
SITE_GROUP_COLUMN <- Sys.getenv("PHENOLOGY_SITE_COLUMN", unset = "Standort")
SITE_COLUMN_CANDIDATES <- unique(c(
  SITE_GROUP_COLUMN, "Standort", "Location"
))

# Run both to check whether the learner ranking is representation-dependent.
# To minimize runtime, set this to "cumulative_weather" only.
FEATURE_SETS_TO_RUN <- c(
  "cumulative_weather",
  "hybrid_process_recent"
)
PRIMARY_FEATURE_SET <- "cumulative_weather"

# Evaluate spatial and temporal transfer during learner selection.
CV_SCHEMES <- c("spatial", "temporal_block")
TEMPORAL_BLOCKS <- 5L
SPATIAL_FOLDS <- 5L
SPATIAL_BLOCK_SIZE_M <- 119946.2
SPATIAL_SELECTION <- "random"

# Inner tuning is site-grouped for every outer job.
INNER_GROUP_FOLDS <- 5L

# Maximum tuning evaluations for the three tuned learners.
TUNING_EVALS_RF <- 40L
TUNING_EVALS_MARS <- 30L
TUNING_EVALS_XGB <- 60L

RF_NUM_TREES <- 1000L
GAM_FAMILY <- "gaussian"  # "scat" or "gaussian"

WORKERS_REQUESTED <- 3L
RESUME_FROM_CHECKPOINTS <- TRUE
SAVE_TUNING_ARCHIVES <- TRUE
CHECKPOINT_VERSION <- 1L

# Use "results_only" after all outer jobs have completed when only the final
# summaries need to be regenerated. This mode reads outer_results_compact.rds
# and never calls run_outer_job(), so it cannot refit or retune a learner.
RESULTS_MODE <- Sys.getenv(
  "PHENOLOGY_LEARNER_RESULTS_MODE",
  unset = "run_or_resume"
)
if (!RESULTS_MODE %in% c("run_or_resume", "results_only")) {
  stop(
    "PHENOLOGY_LEARNER_RESULTS_MODE must be either 'run_or_resume' or ",
    "'results_only'."
  )
}

BOOTSTRAP_REPLICATES <- 2000L

# Predictions are never clipped. These thresholds only flag numerical failure.
EXTREME_PREDICTION_LOWER <- -100
EXTREME_PREDICTION_UPPER <- 200

# ------------------------------------------------------------------------------
# 2. PACKAGES AND REPRODUCIBILITY
# ------------------------------------------------------------------------------

required_packages <- c(
  "data.table", "R6", "mlr3", "mlr3learners", "mlr3extralearners",
  "mlr3tuning", "mlr3mbo", "paradox", "future", "future.apply",
  "parallelly", "mgcv", "earth", "ranger", "xgboost",
  "mlr3spatiotempcv", "sf", "blockCV"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install these packages before running the benchmark: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(R6)
  library(mlr3)
  library(mlr3learners)
  library(mlr3extralearners)
  library(mlr3tuning)
  library(mlr3mbo)
  library(paradox)
  library(future)
  library(future.apply)
  library(parallelly)
  library(mgcv)
  library(mlr3spatiotempcv)
})

RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

workers <- max(
  1L,
  min(WORKERS_REQUESTED, parallelly::availableCores(omit = 1L))
)

# Outer jobs are parallel. Each fitted learner is restricted to one thread.
data.table::setDTthreads(1L)

# ------------------------------------------------------------------------------
# 3. PRE-SPECIFIED PREDICTOR SETS AND GAM FORMULAS
# ------------------------------------------------------------------------------

all_feature_sets <- list(
  cumulative_weather = c(
    "cumDTU_beta",
    "sowing_DOY",
    "photoperiod_civil",
    "cumRR",
    "cumRad",
    "cumCDD",
    "cumHDD",
    "paw_fc33_0_100_mm"
  ),
  
  hybrid_process_recent = c(
    "cumdevfun_phase",
    "sowing_DOY",
    "FTSW_0_100_proxy",
    "paw_fc33_0_100_mm",
    "sumRR_30",
    "sumRad_30",
    "sumCDD_30",
    "sumHDD_30",
    "sumPEvap_30",
    "climaticWaterBalance_30",
    "meanSoilWater_0_100_30"
  )
)


unknown_feature_sets <- setdiff(
  FEATURE_SETS_TO_RUN,
  names(all_feature_sets)
)
if (length(unknown_feature_sets)) {
  stop(
    "Unknown FEATURE_SETS_TO_RUN values: ",
    paste(unknown_feature_sets, collapse = ", ")
  )
}
if (!PRIMARY_FEATURE_SET %in% FEATURE_SETS_TO_RUN) {
  stop("PRIMARY_FEATURE_SET must be included in FEATURE_SETS_TO_RUN.")
}

feature_sets <- all_feature_sets[FEATURE_SETS_TO_RUN]

# k = 50 is used only for the main developmental clock. Smaller marginal bases
# and interactions reduce unsupported edge behavior. Every task feature appears
# in its corresponding GAM formula.
all_gam_formulas <- list(
  cumulative_weather =
    BBCH ~
    s(cumDTU_beta, k = 50, bs = "ts") +
    s(sowing_DOY, k = 12, bs = "ts") +
    s(photoperiod_civil, k = 12, bs = "ts") +
    s(cumRR, k = 12, bs = "ts") +
    s(cumRad, k = 12, bs = "ts") +
    s(cumCDD, k = 12, bs = "ts") +
    s(cumHDD, k = 12, bs = "ts") +
    s(paw_fc33_0_100_mm, k = 12, bs = "ts") +
    ti(
      cumDTU_beta, paw_fc33_0_100_mm,
      k = c(12, 5), bs = c("ts", "ts")
    )+
    ti(
      cumDTU_beta, sowing_DOY,
      k = c(12, 5), bs = c("ts", "ts")
    ),
  
  hybrid_process_recent =
    BBCH ~
    s(cumdevfun_phase, k = 50, bs = "ts") +
    s(sowing_DOY, k = 12, bs = "ts") +
    s(FTSW_0_100_proxy, k = 12, bs = "ts") +
    s(paw_fc33_0_100_mm, k = 12, bs = "ts") +
    s(sumRR_30, k = 12, bs = "ts") +
    s(sumRad_30, k = 12, bs = "ts") +
    s(sumCDD_30, k = 12, bs = "ts") +
    s(sumHDD_30, k = 12, bs = "ts") +
    s(sumPEvap_30, k = 12, bs = "ts") +
    s(climaticWaterBalance_30, k = 12, bs = "ts") +
    s(meanSoilWater_0_100_30, k = 12, bs = "ts") +
    ti(
      cumdevfun_phase,
      FTSW_0_100_proxy,
      k = c(12, 5),
      bs = c("ts", "ts")) +
    ti(
      cumdevfun_phase,
      sowing_DOY,
      k = c(12, 5),
      bs = c("ts", "ts")
    )
)

gam_formulas <- all_gam_formulas[FEATURE_SETS_TO_RUN]
stopifnot(identical(names(feature_sets), names(gam_formulas)))

# ------------------------------------------------------------------------------
# 4. UTILITIES
# ------------------------------------------------------------------------------

safe_id <- function(x) {
  gsub("[^A-Za-z0-9_-]+", "_", as.character(x))
}

resolve_column <- function(data, candidates, label) {
  hits <- candidates[candidates %in% names(data)]
  if (!length(hits)) {
    stop(
      "Could not resolve ", label, ". Expected one of: ",
      paste(candidates, collapse = ", ")
    )
  }
  hits[1L]
}

params_to_text <- function(values) {
  if (is.null(values) || !length(values)) return(NA_character_)

  value_text <- vapply(values, function(value) {
    if (inherits(value, "formula")) {
      paste(deparse(value, width.cutoff = 500L), collapse = " ")
    } else if (inherits(value, "family")) {
      value$family
    } else if (is.atomic(value) && length(value) <= 20L) {
      paste(value, collapse = ",")
    } else {
      paste(class(value), collapse = "/")
    }
  }, character(1))

  paste(paste0(names(value_text), "=", value_text), collapse = "; ")
}

prediction_to_table <- function(prediction) {
  row_ids <- prediction$row_ids
  truth <- prediction$truth
  response <- prediction$response

  if (is.null(row_ids) || is.null(truth) || is.null(response)) {
    stop(
      "Prediction object does not expose row_ids, truth and response. Class: ",
      paste(class(prediction), collapse = "/")
    )
  }

  if (is.data.frame(response)) response <- as.matrix(response)
  if (is.matrix(response) && ncol(response) != 1L) {
    stop(
      "Expected one response column, but prediction returned ",
      ncol(response), "."
    )
  }

  response <- as.numeric(response)
  truth <- as.numeric(truth)

  if (
    length(row_ids) != length(truth) ||
      length(truth) != length(response)
  ) {
    stop(
      "Prediction fields have inconsistent lengths: row_ids=",
      length(row_ids), ", truth=", length(truth),
      ", response=", length(response), "."
    )
  }
  if (
    anyNA(row_ids) ||
      any(!is.finite(truth)) ||
      any(!is.finite(response))
  ) {
    stop("Prediction contains missing or non-finite values.")
  }

  data.table(
    row_ids = as.integer(row_ids),
    truth = truth,
    response = response
  )
}

write_table <- function(x, filename) {
  fwrite(x, file.path(OUTPUT_DIR, filename), na = "NA")
}

make_contiguous_year_blocks <- function(year, k) {
  unique_years <- sort(unique(as.integer(year)))
  if (length(unique_years) < k) {
    stop(
      "Only ", length(unique_years),
      " season years are available for ", k, " temporal blocks."
    )
  }

  year_block <- ceiling(
    seq_along(unique_years) * k / length(unique_years)
  )
  year_block <- pmin(year_block, k)
  lookup <- setNames(year_block, unique_years)
  unname(lookup[as.character(as.integer(year))])
}

capture_condition <- function(store, condition) {
  unique(c(store, conditionMessage(condition)))
}

# ------------------------------------------------------------------------------
# 5. LOAD DATA AND DEFINE ONE COMMON COMPARISON COHORT
# ------------------------------------------------------------------------------

if (!file.exists(INPUT_FILE)) {
  stop("Input file does not exist: ", INPUT_FILE)
}

message("Reading: ", INPUT_FILE)
raw_data <- fread(
  INPUT_FILE,
  na.strings = c("", "NA", "NaN", "Inf", "-Inf"),
  showProgress = TRUE
)
raw_data[, .source_row_id := .I]

required_features <- unique(unlist(feature_sets, use.names = FALSE))

site_column <- resolve_column(
  raw_data,
  SITE_COLUMN_CANDIDATES,
  "physical-site grouping column"
)
environment_column <- resolve_column(
  raw_data,
  c("env_id", "traj_id"),
  "site-season environment ID"
)
year_column <- resolve_column(
  raw_data,
  c("SeasonYear", "season_year", "year"),
  "season year"
)
date_column <- resolve_column(
  raw_data,
  c("Date", "phen_date", "date"),
  "observation date"
)
longitude_column <- resolve_column(
  raw_data,
  c("lon", "longitude"),
  "longitude"
)
latitude_column <- resolve_column(
  raw_data,
  c("lat", "latitude"),
  "latitude"
)

missing_model_columns <- setdiff(
  c(TARGET, required_features),
  names(raw_data)
)
if (length(missing_model_columns)) {
  stop(
    "Required model columns are absent: ",
    paste(missing_model_columns, collapse = ", ")
  )
}

raw_data[, .site_id := trimws(as.character(get(site_column)))]
raw_data[, .env_id := trimws(as.character(get(environment_column)))]
raw_data[, .season_year := as.integer(get(year_column))]
raw_data[, .date := as.IDate(get(date_column))]
raw_data[, .lon := as.numeric(get(longitude_column))]
raw_data[, .lat := as.numeric(get(latitude_column))]

raw_data[!nzchar(.site_id), .site_id := NA_character_]
raw_data[!nzchar(.env_id), .env_id := NA_character_]

required_common_columns <- unique(c(
  TARGET,
  required_features,
  ".site_id",
  ".env_id",
  ".season_year",
  ".date",
  ".lon",
  ".lat"
))

non_numeric_features <- required_features[
  !vapply(raw_data[, ..required_features], is.numeric, logical(1))
]
if (length(non_numeric_features)) {
  stop(
    "All predictors must be numeric. Non-numeric predictors: ",
    paste(non_numeric_features, collapse = ", ")
  )
}

missingness_audit <- rbindlist(lapply(
  required_common_columns,
  function(column) {
    value <- raw_data[[column]]
    data.table(
      column = column,
      n_missing = sum(is.na(value)),
      proportion_missing = mean(is.na(value))
    )
  }
))

complete_mask <- complete.cases(raw_data[, ..required_common_columns])
finite_columns <- unique(c(
  TARGET,
  required_features,
  ".season_year",
  ".lon",
  ".lat"
))
for (column in finite_columns) {
  complete_mask <- complete_mask & is.finite(raw_data[[column]])
}

cohort <- copy(raw_data[complete_mask])
cohort[, .obs_id := .I]

cohort_flow <- data.table(
  step = c("input", "complete_common_cohort", "removed"),
  n_rows = c(
    nrow(raw_data),
    nrow(cohort),
    nrow(raw_data) - nrow(cohort)
  ),
  proportion_of_input = c(
    1,
    nrow(cohort) / nrow(raw_data),
    (nrow(raw_data) - nrow(cohort)) / nrow(raw_data)
  )
)

write_table(missingness_audit, "missingness_required_columns.csv")
write_table(cohort_flow, "cohort_flow.csv")

if (nrow(cohort) < 100L) {
  stop("The common cohort contains fewer than 100 rows.")
}
if (uniqueN(cohort$.site_id) < INNER_GROUP_FOLDS) {
  stop("Too few physical sites for grouped inner CV.")
}

# ------------------------------------------------------------------------------
# 6. INPUT QUALITY CONTROL
# ------------------------------------------------------------------------------

if (any(cohort[[TARGET]] < 0 | cohort[[TARGET]] > 99)) {
  stop("BBCH contains values outside 0-99.")
}
if (any(cohort$sowing_DOY < 1 | cohort$sowing_DOY > 366)) {
  stop("sowing_DOY contains values outside 1-366.")
}
if (any(cohort$photoperiod_civil <= 0 | cohort$photoperiod_civil > 24)) {
  stop("photoperiod_civil contains values outside (0, 24] hours.")
}
if (any(cohort$cumDTU_beta < 0)) {
  stop("cumDTU_beta contains negative values.")
}
if (
  "FTSW_0_28_proxy" %in% names(cohort) &&
    any(
      cohort$FTSW_0_28_proxy < -1e-6 |
        cohort$FTSW_0_28_proxy > 1 + 1e-6
    )
) {
  stop(
    "FTSW_0_28_proxy is outside [0, 1]. Fix its scaling upstream; ",
    "do not silently clip it in the benchmark."
  )
}

environment_mapping_qc <- cohort[, .(
  n_sites = uniqueN(.site_id),
  n_years = uniqueN(.season_year)
), by = .env_id]

invalid_environment_ids <- environment_mapping_qc[
  n_sites != 1L | n_years != 1L
]
write_table(environment_mapping_qc, "environment_mapping_qc.csv")
if (nrow(invalid_environment_ids)) {
  write_table(invalid_environment_ids, "invalid_environment_ids.csv")
  stop(
    "At least one environment maps to multiple physical sites or years. ",
    "See invalid_environment_ids.csv."
  )
}

site_coordinate_qc <- cohort[, .(
  n_longitudes = uniqueN(round(.lon, 5)),
  n_latitudes = uniqueN(round(.lat, 5)),
  longitude_range = diff(range(.lon)),
  latitude_range = diff(range(.lat))
), by = .site_id]
write_table(site_coordinate_qc, "site_coordinate_qc.csv")

gdd_qc <- copy(cohort)[order(.env_id, .date), .(
  n_gdd_decreases = if (.N > 1L) {
    sum(diff(cumDTU_beta) < -1e-8)
  } else {
    0L
  },
  n_clock_decreases = if (
    .N > 1L && "cumdevfun_phase" %in% names(.SD)
  ) {
    sum(diff(cumdevfun_phase) < -1e-8)
  } else {
    0L
  }
), by = .env_id]

if (any(gdd_qc$n_gdd_decreases > 0L)) {
  write_table(
    gdd_qc[n_gdd_decreases > 0L],
    "non_monotonic_gdd_trajectories.csv"
  )
  stop("Non-monotonic cumulative GDD trajectories were detected.")
}
if (any(gdd_qc$n_clock_decreases > 0L)) {
  write_table(
    gdd_qc[n_clock_decreases > 0L],
    "non_monotonic_process_clock_trajectories.csv"
  )
  stop("Non-monotonic cumulative process-clock trajectories were detected.")
}

feature_qc <- rbindlist(lapply(required_features, function(column) {
  x <- cohort[[column]]
  data.table(
    feature = column,
    n_unique = uniqueN(x),
    minimum = min(x),
    q01 = as.numeric(quantile(x, 0.01)),
    median = median(x),
    q99 = as.numeric(quantile(x, 0.99)),
    maximum = max(x),
    sd = sd(x)
  )
}))
write_table(feature_qc, "feature_qc_common_cohort.csv")

constant_features <- feature_qc[n_unique < 2L, feature]
if (length(constant_features)) {
  stop(
    "Constant predictors in the common cohort: ",
    paste(constant_features, collapse = ", ")
  )
}

feature_set_definitions <- rbindlist(lapply(
  names(feature_sets),
  function(feature_set) {
    data.table(
      feature_set = feature_set,
      feature = feature_sets[[feature_set]],
      role = if (
        identical(feature_set, PRIMARY_FEATURE_SET)
      ) {
        "primary"
      } else {
        "sensitivity"
      }
    )
  }
))
write_table(feature_set_definitions, "feature_set_definitions.csv")

# ------------------------------------------------------------------------------
# 7. MLR3 TASKS
# ------------------------------------------------------------------------------

make_task <- function(feature_set, features) {
  backend_columns <- unique(c(
    ".obs_id",
    ".env_id",
    ".site_id",
    TARGET,
    ".lon",
    ".lat",
    features
  ))
  backend <- cohort[, ..backend_columns]

  task <- TaskRegrST$new(
    id = feature_set,
    backend = backend,
    target = TARGET,
    coordinate_names = c(".lon", ".lat"),
    crs = "EPSG:4326"
  )

  task$select(features)
  task
}

tasks <- Map(make_task, names(feature_sets), feature_sets)
names(tasks) <- names(feature_sets)

reference_row_ids <- tasks[[PRIMARY_FEATURE_SET]]$row_ids
if (!all(vapply(
  tasks,
  function(task) identical(task$row_ids, reference_row_ids),
  logical(1)
))) {
  stop("Predictor-set tasks do not have identical row IDs.")
}

# ------------------------------------------------------------------------------
# 8. REUSABLE OUTER FOLDS
# ------------------------------------------------------------------------------

fold_vectors <- list()

if ("temporal_block" %in% CV_SCHEMES) {
  fold_vectors$temporal_block <- make_contiguous_year_blocks(
    cohort$.season_year,
    TEMPORAL_BLOCKS
  )
}

if ("spatial" %in% CV_SCHEMES) {
  set.seed(SEED + 20L)
  spatial_resampling <- rsmp(
    "spcv_block",
    range = as.integer(SPATIAL_BLOCK_SIZE_M),
    folds = as.integer(SPATIAL_FOLDS),
    selection = SPATIAL_SELECTION,
    hexagon = FALSE
  )
  spatial_resampling$instantiate(tasks[[PRIMARY_FEATURE_SET]])

  spatial_instance <- as.data.table(spatial_resampling$instance)
  row_id_column <- intersect(
    c("row_id", "row_ids"),
    names(spatial_instance)
  )[1L]
  if (!length(row_id_column) || is.na(row_id_column)) {
    stop("Could not identify row IDs in the spatial resampling instance.")
  }

  fold_vectors$spatial <- spatial_instance$fold[
    match(reference_row_ids, spatial_instance[[row_id_column]])
  ]
  if (anyNA(fold_vectors$spatial)) {
    stop("Spatial CV did not assign every observation to a fold.")
  }
}

if (!identical(sort(names(fold_vectors)), sort(CV_SCHEMES))) {
  stop("One or more requested CV schemes were not constructed.")
}

for (scheme in names(fold_vectors)) {
  environment_check <- data.table(
    env_id = cohort$.env_id,
    fold = as.character(fold_vectors[[scheme]])
  )[, .(n_folds = uniqueN(fold)), by = env_id]

  if (any(environment_check$n_folds != 1L)) {
    stop("An environment was split across outer folds in: ", scheme)
  }
}

if ("spatial" %in% names(fold_vectors)) {
  site_check <- data.table(
    site_id = cohort$.site_id,
    fold = as.character(fold_vectors$spatial)
  )[, .(n_folds = uniqueN(fold)), by = site_id]

  if (any(site_check$n_folds != 1L)) {
    write_table(
      site_check[n_folds != 1L],
      "spatial_sites_split_across_folds.csv"
    )
    stop("Spatial CV split at least one physical site across folds.")
  }
}

fold_manifest <- cohort[, .(
  obs_id = .obs_id,
  source_row_id = .source_row_id,
  env_id = .env_id,
  site_id = .site_id,
  season_year = .season_year,
  Date = .date,
  BBCH = get(TARGET),
  lon = .lon,
  lat = .lat
)]
for (scheme in names(fold_vectors)) {
  fold_manifest[
    ,
    (paste0("fold_", scheme)) := as.character(fold_vectors[[scheme]])
  ]
}
write_table(fold_manifest, "outer_fold_manifest.csv")

fold_sizes <- rbindlist(lapply(names(fold_vectors), function(scheme) {
  answer <- data.table(
    fold = as.character(fold_vectors[[scheme]])
  )[, .N, by = fold]
  answer[, scheme := scheme]
  answer
}))
setcolorder(fold_sizes, c("scheme", "fold", "N"))
write_table(fold_sizes, "outer_fold_sizes.csv")

# ------------------------------------------------------------------------------
# 9. CUSTOM MGCV::gam LEARNER
# ------------------------------------------------------------------------------

LearnerRegrGamMgcv <- R6Class(
  "LearnerRegrGamMgcv",
  inherit = LearnerRegr,

  public = list(
    initialize = function(
      formula,
      family,
      id = "gam"
    ) {
      parameter_set <- ps(
        formula = p_uty(),
        family = p_uty()
      )
      super$initialize(
        id = id,
        feature_types = c("numeric", "integer"),
        predict_types = "response",
        packages = "mgcv",
        param_set = parameter_set
      )
      self$param_set$values <- list(
        formula = formula,
        family = family
      )
    }
  ),

  private = list(
    .train = function(task) {
      values <- self$param_set$values

      mgcv::gam(
        formula = values$formula,
        data = task$data(),
        family = values$family,
        method = "REML",
        select = TRUE
      )
    },

    .predict = function(task) {
      newdata <- task$data(cols = task$feature_names)
      response <- predict(
        self$model,
        newdata = newdata,
        type = "response"
      )
      list(response = as.numeric(response))
    }
  )
)

gam_family <- switch(
  GAM_FAMILY,
  scat = mgcv::scat(),
  gaussian = stats::gaussian(),
  stop("GAM_FAMILY must be 'scat' or 'gaussian'.")
)

# ------------------------------------------------------------------------------
# 10. LEARNERS, NESTED TUNING AND SEARCH SPACES
# ------------------------------------------------------------------------------

make_auto_tuner <- function(
  learner,
  search_space,
  evaluations,
  id
) {
  auto_tuner(
    learner = learner,
    resampling = rsmp("cv", folds = INNER_GROUP_FOLDS),
    measure = msr("regr.mae"),
    tuner = tnr("mbo"),
    search_space = search_space,
    term_evals = evaluations,
    store_tuning_instance = TRUE,
    # The optimization archive is retained below without storing every inner
    # BenchmarkResult, which would make checkpoints unnecessarily large.
    store_benchmark_result = FALSE,
    store_models = FALSE,
    check_values = TRUE,
    id = id
  )
}

make_learners <- function(feature_set, task) {
  p <- length(task$feature_names)

  learner_gam <- LearnerRegrGamMgcv$new(
    formula = gam_formulas[[feature_set]],
    family = gam_family,
    id = "gam"
  )

  learner_rf <- lrn(
    "regr.ranger",
    num.trees = RF_NUM_TREES,
    importance = "none",
    num.threads = 1L,
    seed = SEED,
    predict_type = "response"
  )
  tuned_rf <- make_auto_tuner(
    learner = learner_rf,
    search_space = ps(
      mtry = p_int(1L, p),
      min.node.size = p_int(1L, 100L, logscale = TRUE),
      sample.fraction = p_dbl(0.5, 1),
      max.depth = p_int(0L, 20L)
    ),
    evaluations = TUNING_EVALS_RF,
    id = "random_forest"
  )

  learner_mars <- lrn(
    "regr.earth",
    predict_type = "response"
  )
  tuned_mars <- make_auto_tuner(
    learner = learner_mars,
    search_space = ps(
      degree = p_int(1L, min(3L, p)),
      penalty = p_dbl(0.5, 5, logscale = TRUE),
      nprune = p_int(
        2L,
        max(10L, min(60L, 3L * p + 5L))
      )
    ),
    evaluations = TUNING_EVALS_MARS,
    id = "mars"
  )

  learner_xgb <- lrn(
    "regr.xgboost",
    objective = "reg:squarederror",
    tree_method = "hist",
    booster = "gbtree",
    nthread = 1L,
    seed = SEED,
    predict_type = "response"
  )
  tuned_xgb <- make_auto_tuner(
    learner = learner_xgb,
    search_space = ps(
      nrounds = p_int(200L, 2500L, logscale = TRUE),
      eta = p_dbl(0.01, 0.3, logscale = TRUE),
      gamma = p_dbl(0, 5),
      max_depth = p_int(1L, 10L),
      min_child_weight = p_dbl(0.1, 30, logscale = TRUE),
      subsample = p_dbl(0.5, 1),
      colsample_bytree = p_dbl(0.5, 1),
      lambda = p_dbl(1e-4, 20, logscale = TRUE),
      alpha = p_dbl(1e-4, 20, logscale = TRUE)
    ),
    evaluations = TUNING_EVALS_XGB,
    id = "xgboost"
  )

  list(
    gam = learner_gam,
    mars = tuned_mars,
    random_forest = tuned_rf,
    xgboost = tuned_xgb
  )
}

learners_by_task <- Map(
  make_learners,
  names(tasks),
  tasks
)
names(learners_by_task) <- names(tasks)

# ------------------------------------------------------------------------------
# 11. EXPLICIT NESTED OUTER-CV JOBS
# ------------------------------------------------------------------------------

jobs <- rbindlist(lapply(names(tasks), function(feature_set) {
  rbindlist(lapply(
    names(learners_by_task[[feature_set]]),
    function(learner_id) {
      rbindlist(lapply(names(fold_vectors), function(scheme) {
        fold_ids <- sort(unique(as.character(fold_vectors[[scheme]])))
        data.table(
          feature_set = feature_set,
          learner_id = learner_id,
          scheme = scheme,
          fold_id = fold_ids,
          iteration = seq_along(fold_ids)
        )
      }))
    }
  ))
}))
jobs[, job_id := .I]
write_table(jobs, "benchmark_jobs.csv")

run_outer_job <- function(job) {
  checkpoint_name <- paste0(
    sprintf("%04d", job$job_id), "__",
    safe_id(job$feature_set), "__",
    safe_id(job$learner_id), "__",
    safe_id(job$scheme), "__fold_",
    safe_id(job$fold_id), ".rds"
  )
  checkpoint_file <- file.path(CHECKPOINT_DIR, checkpoint_name)

  if (RESUME_FROM_CHECKPOINTS && file.exists(checkpoint_file)) {
    cached_result <- tryCatch(
      readRDS(checkpoint_file),
      error = function(error) NULL
    )
    if (
      !is.null(cached_result) &&
        identical(
          cached_result$checkpoint_version,
          CHECKPOINT_VERSION
        )
    ) {
      return(cached_result)
    }
  }

  set.seed(SEED + as.integer(job$job_id))

  task <- tasks[[job$feature_set]]$clone(deep = TRUE)
  learner <- learners_by_task[[job$feature_set]][[job$learner_id]]$clone(
    deep = TRUE
  )

  fold_vector <- as.character(fold_vectors[[job$scheme]])
  test_ids <- which(fold_vector == as.character(job$fold_id))
  train_ids <- setdiff(task$row_ids, test_ids)

  if (!length(train_ids) || !length(test_ids)) {
    stop("Empty outer training or test set in job ", job$job_id)
  }
  if (length(intersect(
    unique(cohort$.env_id[train_ids]),
    unique(cohort$.env_id[test_ids])
  ))) {
    stop("Environment leakage in outer job ", job$job_id)
  }
  if (
    identical(job$scheme, "spatial") &&
      length(intersect(
        unique(cohort$.site_id[train_ids]),
        unique(cohort$.site_id[test_ids])
      ))
  ) {
    stop("Physical-site leakage in spatial job ", job$job_id)
  }
  if (uniqueN(cohort$.site_id[train_ids]) < INNER_GROUP_FOLDS) {
    stop(
      "Outer training set has too few sites for inner grouped CV in job ",
      job$job_id
    )
  }

  warning_messages <- character()
  message_messages <- character()
  start_time <- proc.time()[["elapsed"]]

  result <- tryCatch(
    withCallingHandlers({
      learner$train(task, row_ids = train_ids)
      prediction <- learner$predict(task, row_ids = test_ids)
      prediction_table <- prediction_to_table(prediction)

      prediction_table[, feature_set := job$feature_set]
      prediction_table[, learner_id := job$learner_id]
      prediction_table[, scheme := job$scheme]
      prediction_table[, fold_id := as.character(job$fold_id)]
      prediction_table[, iteration := as.integer(job$iteration)]

      prediction_error <- (
        prediction_table$response - prediction_table$truth
      )

      tuned_parameters <- NA_character_
      inner_mae <- NA_real_
      tuning_archive <- NULL

      if (inherits(learner, "AutoTuner")) {
        tuned_parameters <- params_to_text(
          learner$learner$param_set$values
        )
        tuning_result <- learner$tuning_result
        if (
          !is.null(tuning_result) &&
            "regr.mae" %in% names(tuning_result)
        ) {
          inner_mae <- tuning_result[["regr.mae"]][1L]
        }
        if (SAVE_TUNING_ARCHIVES) {
          tuning_archive <- as.data.table(learner$archive$data)
        }
      } else {
        tuned_parameters <- params_to_text(
          learner$param_set$values
        )
      }

      list(
        status = "ok",
        job = copy(job),
        predictions = prediction_table,
        fold_metrics = data.table(
          feature_set = job$feature_set,
          learner_id = job$learner_id,
          scheme = job$scheme,
          fold_id = as.character(job$fold_id),
          iteration = as.integer(job$iteration),
          n_train = length(train_ids),
          n_test = nrow(prediction_table),
          rmse = sqrt(mean(prediction_error^2)),
          mae = mean(abs(prediction_error)),
          bias = mean(prediction_error),
          within_5_bbch = mean(abs(prediction_error) <= 5),
          n_outside_0_99 = sum(
            prediction_table$response < 0 |
              prediction_table$response > 99
          ),
          n_extreme = sum(
            prediction_table$response < EXTREME_PREDICTION_LOWER |
              prediction_table$response > EXTREME_PREDICTION_UPPER
          ),
          min_prediction = min(prediction_table$response),
          max_prediction = max(prediction_table$response),
          inner_mae = inner_mae,
          tuned_parameters = tuned_parameters
        ),
        tuning_archive = tuning_archive,
        warnings = warning_messages,
        messages = message_messages,
        error = NA_character_
      )
    },
    warning = function(warning) {
      warning_messages <<- capture_condition(
        warning_messages,
        warning
      )
      invokeRestart("muffleWarning")
    },
    message = function(message_condition) {
      message_messages <<- capture_condition(
        message_messages,
        message_condition
      )
      invokeRestart("muffleMessage")
    }),
    error = function(error) {
      list(
        status = "error",
        job = copy(job),
        predictions = NULL,
        fold_metrics = data.table(
          feature_set = job$feature_set,
          learner_id = job$learner_id,
          scheme = job$scheme,
          fold_id = as.character(job$fold_id),
          iteration = as.integer(job$iteration),
          n_train = length(train_ids),
          n_test = length(test_ids),
          rmse = NA_real_,
          mae = NA_real_,
          bias = NA_real_,
          within_5_bbch = NA_real_,
          n_outside_0_99 = NA_integer_,
          n_extreme = NA_integer_,
          min_prediction = NA_real_,
          max_prediction = NA_real_,
          inner_mae = NA_real_,
          tuned_parameters = NA_character_
        ),
        tuning_archive = NULL,
        warnings = warning_messages,
        messages = message_messages,
        error = conditionMessage(error)
      )
    }
  )

  result$checkpoint_version <- CHECKPOINT_VERSION
  result$elapsed_seconds <- proc.time()[["elapsed"]] - start_time

  temporary_checkpoint <- paste0(
    checkpoint_file,
    ".tmp_",
    Sys.getpid()
  )
  saveRDS(result, temporary_checkpoint, compress = FALSE)
  if (!file.rename(temporary_checkpoint, checkpoint_file)) {
    stop("Could not atomically write checkpoint: ", checkpoint_file)
  }

  result
}

outer_results_file <- file.path(
  OUTPUT_DIR,
  "outer_results_compact.rds"
)

if (identical(RESULTS_MODE, "results_only")) {
  if (!file.exists(outer_results_file)) {
    stop(
      "Results-only mode requires the completed result object: ",
      outer_results_file
    )
  }

  message(
    "Results-only mode: loading completed outer-fold results; no learner ",
    "fitting or tuning will run."
  )
  outer_results <- readRDS(outer_results_file)

  if (!is.list(outer_results) || length(outer_results) != nrow(jobs)) {
    stop(
      "The saved outer result contains ", length(outer_results),
      " jobs, but the current configuration expects ", nrow(jobs), "."
    )
  }

  saved_job_ids <- vapply(
    outer_results,
    function(result) as.integer(result$job$job_id[[1L]]),
    integer(1)
  )
  saved_checkpoint_versions <- vapply(
    outer_results,
    function(result) as.integer(result$checkpoint_version[[1L]]),
    integer(1)
  )

  if (!identical(unname(saved_job_ids), as.integer(jobs$job_id))) {
    stop(
      "The job identities in outer_results_compact.rds do not match the ",
      "current benchmark configuration."
    )
  }
  if (any(saved_checkpoint_versions != CHECKPOINT_VERSION)) {
    stop(
      "The checkpoint version in outer_results_compact.rds does not match ",
      "CHECKPOINT_VERSION."
    )
  }
} else {
  message(
    "Starting ", nrow(jobs), " outer-fold jobs with ",
    workers, " worker(s)."
  )
  message("Checkpoints: ", CHECKPOINT_DIR)

  options(future.globals.maxSize = 8 * 1024^3)
  future::plan(future::multisession, workers = workers)
  outer_results <- future.apply::future_lapply(
    seq_len(nrow(jobs)),
    function(i) run_outer_job(jobs[i]),
    future.seed = TRUE,
    future.scheduling = 1
  )
  future::plan(future::sequential)

  saveRDS(
    outer_results,
    outer_results_file,
    compress = FALSE
  )
}

# ------------------------------------------------------------------------------
# 12. COMBINE PREDICTIONS, CONDITIONS AND FOLD RESULTS
# ------------------------------------------------------------------------------

job_status <- rbindlist(lapply(outer_results, function(result) {
  data.table(
    job_id = result$job$job_id,
    feature_set = result$job$feature_set,
    learner_id = result$job$learner_id,
    scheme = result$job$scheme,
    fold_id = as.character(result$job$fold_id),
    status = result$status,
    elapsed_seconds = result$elapsed_seconds,
    error = result$error,
    warnings = paste(result$warnings, collapse = " | "),
    messages = paste(result$messages, collapse = " | ")
  )
}), fill = TRUE)
setorder(job_status, job_id)
write_table(job_status, "job_status.csv")

failed_jobs <- job_status[status != "ok"]
write_table(failed_jobs, "failed_jobs.csv")

fold_metrics <- rbindlist(
  lapply(outer_results, `[[`, "fold_metrics"),
  fill = TRUE
)
write_table(fold_metrics, "outer_fold_metrics.csv")

successful_predictions <- lapply(
  outer_results,
  function(result) result$predictions
)
successful_predictions <- successful_predictions[
  !vapply(successful_predictions, is.null, logical(1))
]

if (!length(successful_predictions)) {
  stop("No outer job produced predictions.")
}

oof_predictions <- rbindlist(
  successful_predictions,
  use.names = TRUE,
  fill = TRUE
)

prediction_metadata <- cohort[, .(
  row_ids = .obs_id,
  source_row_id = .source_row_id,
  env_id = .env_id,
  site_id = .site_id,
  season_year = .season_year,
  Date = .date
)]

oof_predictions <- merge(
  oof_predictions,
  prediction_metadata,
  by = "row_ids",
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(oof_predictions$env_id)) {
  stop("Could not match all predictions back to cohort metadata.")
}

oof_predictions[, error := response - truth]
oof_predictions[, absolute_error := abs(error)]
oof_predictions[, squared_error := error^2]
oof_predictions[, outside_0_99 := response < 0 | response > 99]
oof_predictions[
  ,
  extreme_prediction :=
    response < EXTREME_PREDICTION_LOWER |
    response > EXTREME_PREDICTION_UPPER
]

setcolorder(
  oof_predictions,
  c(
    "feature_set", "learner_id", "scheme", "fold_id",
    "iteration", "row_ids", "source_row_id", "env_id",
    "site_id", "season_year", "Date", "truth", "response",
    "error", "absolute_error", "squared_error",
    "outside_0_99", "extreme_prediction"
  )
)
write_table(oof_predictions, "oof_predictions.csv")

if (SAVE_TUNING_ARCHIVES) {
  tuning_archives <- rbindlist(lapply(outer_results, function(result) {
    if (is.null(result$tuning_archive)) return(NULL)
    archive <- copy(result$tuning_archive)
    archive[, job_id := result$job$job_id]
    archive[, feature_set := result$job$feature_set]
    archive[, learner_id := result$job$learner_id]
    archive[, scheme := result$job$scheme]
    archive[, fold_id := as.character(result$job$fold_id)]
    archive
  }), fill = TRUE)

  saveRDS(
    tuning_archives,
    file.path(OUTPUT_DIR, "tuning_archives.rds"),
    compress = FALSE
  )

  # Some mlr3 versions include nested list columns such as x_domain. Retain the
  # complete archive in RDS and write only scalar columns to portable CSV.
  archive_list_columns <- names(tuning_archives)[vapply(
    tuning_archives,
    is.list,
    logical(1)
  )]
  tuning_archives_csv <- tuning_archives[
    ,
    setdiff(names(tuning_archives), archive_list_columns),
    with = FALSE
  ]
  write_table(tuning_archives_csv, "tuning_archives_scalar.csv")
}

# ------------------------------------------------------------------------------
# 13. PERFORMANCE AND STABILITY SUMMARIES
# ------------------------------------------------------------------------------

overall_metrics <- oof_predictions[, .(
  n = .N,
  rmse = sqrt(mean(squared_error)),
  mae = mean(absolute_error),
  bias = mean(error),
  within_5_bbch = mean(absolute_error <= 5),
  n_outside_0_99 = sum(outside_0_99),
  proportion_outside_0_99 = mean(outside_0_99),
  n_extreme = sum(extreme_prediction),
  min_prediction = min(response),
  max_prediction = max(response)
), by = .(feature_set, learner_id, scheme)]

stage_metrics <- oof_predictions[, .(
  n = .N,
  rmse = sqrt(mean(squared_error)),
  mae = mean(absolute_error),
  bias = mean(error),
  within_5_bbch = mean(absolute_error <= 5)
), by = .(
  feature_set,
  learner_id,
  scheme,
  BBCH_stage = truth
)]

# Each environment-by-observed-stage unit receives equal weight.
environment_stage_metrics <- oof_predictions[, .(
  n = .N,
  mse = mean(squared_error),
  mae = mean(absolute_error),
  bias = mean(error)
), by = .(
  feature_set,
  learner_id,
  scheme,
  env_id,
  BBCH_stage = truth
)]
environment_stage_metrics[, rmse := sqrt(mse)]

macro_environment_stage <- environment_stage_metrics[, .(
  n_environment_stage_units = .N,
  macro_mae = mean(mae),
  macro_rmse = sqrt(mean(mse)),
  macro_absolute_bias = mean(abs(bias))
), by = .(feature_set, learner_id, scheme)]

stage_extremes <- stage_metrics[, .(
  n_stages = .N,
  worst_stage_mae = max(mae),
  worst_stage_rmse = max(rmse),
  worst_mae_stage = BBCH_stage[which.max(mae)],
  worst_rmse_stage = BBCH_stage[which.max(rmse)]
), by = .(feature_set, learner_id, scheme)]

performance_summary <- Reduce(
  function(left, right) {
    merge(
      left,
      right,
      by = c("feature_set", "learner_id", "scheme"),
      all = TRUE
    )
  },
  list(
    overall_metrics,
    macro_environment_stage,
    stage_extremes
  )
)

performance_summary[, rank_macro_mae := frank(
  macro_mae,
  ties.method = "average"
), by = .(feature_set, scheme)]
performance_summary[, rank_macro_rmse := frank(
  macro_rmse,
  ties.method = "average"
), by = .(feature_set, scheme)]

setorder(
  performance_summary,
  feature_set,
  scheme,
  rank_macro_mae
)

write_table(overall_metrics, "overall_row_weighted_metrics.csv")
write_table(stage_metrics, "stage_specific_metrics.csv")
write_table(
  environment_stage_metrics,
  "environment_stage_metrics.csv"
)
write_table(performance_summary, "performance_summary.csv")

extreme_folds <- fold_metrics[
  n_extreme > 0L | n_outside_0_99 > 0L
][order(-n_extreme, -n_outside_0_99)]
write_table(extreme_folds, "extreme_prediction_folds.csv")

# Check whether predicted BBCH decreases across successive observation dates.
trajectory_diagnostics <- oof_predictions[
  order(feature_set, learner_id, scheme, env_id, Date),
  .(
    n_observations = .N,
    n_prediction_decreases = if (.N > 1L) {
      sum(diff(response) < -1e-8)
    } else {
      0L
    },
    n_decreases_gt_1 = if (.N > 1L) {
      sum(diff(response) < -1)
    } else {
      0L
    },
    largest_prediction_decrease = if (.N > 1L) {
      min(c(0, diff(response)))
    } else {
      0
    }
  ),
  by = .(feature_set, learner_id, scheme, env_id)
]

trajectory_summary <- trajectory_diagnostics[, .(
  n_environments = .N,
  environments_with_any_decrease = sum(
    n_prediction_decreases > 0L
  ),
  environments_with_decrease_gt_1 = sum(
    n_decreases_gt_1 > 0L
  ),
  proportion_environments_with_decrease_gt_1 = mean(
    n_decreases_gt_1 > 0L
  ),
  total_decreases_gt_1 = sum(n_decreases_gt_1),
  largest_prediction_decrease = min(largest_prediction_decrease)
), by = .(feature_set, learner_id, scheme)]

write_table(
  trajectory_diagnostics,
  "trajectory_diagnostics_by_environment.csv"
)
write_table(
  trajectory_summary,
  "trajectory_diagnostics_summary.csv"
)

# ------------------------------------------------------------------------------
# 14. PRIMARY LEARNER SELECTION
# ------------------------------------------------------------------------------

primary_performance <- performance_summary[
  feature_set == PRIMARY_FEATURE_SET
]

job_completion <- job_status[
  feature_set == PRIMARY_FEATURE_SET,
  .(
    expected_outer_jobs = .N,
    completed_outer_jobs = sum(status == "ok"),
    failed_outer_jobs = sum(status != "ok")
  ),
  by = learner_id
]

primary_stability <- overall_metrics[
  feature_set == PRIMARY_FEATURE_SET,
  .(
    n_extreme_predictions = sum(n_extreme),
    n_outside_0_99 = sum(n_outside_0_99),
    minimum_prediction = min(min_prediction),
    maximum_prediction = max(max_prediction)
  ),
  by = learner_id
]

learner_selection_summary <- primary_performance[, .(
  n_schemes = .N,
  mean_rank_macro_mae = mean(rank_macro_mae),
  worst_rank_macro_mae = max(rank_macro_mae),
  mean_macro_mae = mean(macro_mae),
  worst_macro_mae = max(macro_mae),
  mean_macro_rmse = mean(macro_rmse),
  worst_macro_rmse = max(macro_rmse),
  worst_stage_mae = max(worst_stage_mae)
), by = learner_id]

learner_selection_summary <- merge(
  learner_selection_summary,
  job_completion,
  by = "learner_id",
  all = TRUE
)
learner_selection_summary <- merge(
  learner_selection_summary,
  primary_stability,
  by = "learner_id",
  all = TRUE
)

# Catastrophic predictions or failed folds make a learner ineligible for an
# automatic recommendation. The raw results remain available for diagnosis.
learner_selection_summary[
  ,
  eligible :=
    !is.na(failed_outer_jobs) &
    failed_outer_jobs == 0L &
    !is.na(n_extreme_predictions) &
    n_extreme_predictions == 0L &
    !is.na(n_schemes) &
    n_schemes == length(CV_SCHEMES)
]

setorder(
  learner_selection_summary,
  -eligible,
  mean_rank_macro_mae,
  worst_macro_mae,
  mean_macro_mae
)
learner_selection_summary[, selection_order := .I]

write_table(
  learner_selection_summary,
  "learner_selection_primary_cumulative.csv"
)

# Hybrid is a sensitivity analysis: positive delta means hybrid is worse.
if ("hybrid_process_recent" %in% names(feature_sets)) {
  cumulative_comparison <- macro_environment_stage[
    feature_set == PRIMARY_FEATURE_SET,
    .(
      learner_id,
      scheme,
      cumulative_macro_mae = macro_mae,
      cumulative_macro_rmse = macro_rmse
    )
  ]
  hybrid_comparison <- macro_environment_stage[
    feature_set == "hybrid_process_recent",
    .(
      learner_id,
      scheme,
      hybrid_macro_mae = macro_mae,
      hybrid_macro_rmse = macro_rmse
    )
  ]

  representation_sensitivity <- merge(
    cumulative_comparison,
    hybrid_comparison,
    by = c("learner_id", "scheme"),
    all = TRUE
  )
  representation_sensitivity[
    ,
    delta_mae_hybrid_minus_cumulative :=
      hybrid_macro_mae - cumulative_macro_mae
  ]
  representation_sensitivity[
    ,
    delta_rmse_hybrid_minus_cumulative :=
      hybrid_macro_rmse - cumulative_macro_rmse
  ]
  setorder(
    representation_sensitivity,
    scheme,
    delta_mae_hybrid_minus_cumulative
  )
  write_table(
    representation_sensitivity,
    "hybrid_sensitivity_vs_cumulative.csv"
  )
}

# ------------------------------------------------------------------------------
# 15. PAIRED CLUSTER BOOTSTRAP BETWEEN LEARNERS
# ------------------------------------------------------------------------------

cluster_bootstrap_delta <- function(
  unit_table,
  replicates,
  seed
) {
  environment_summary <- unit_table[, .(
    delta_sum = sum(delta),
    n_units = .N
  ), by = env_id]

  if (nrow(environment_summary) < 2L) {
    return(data.table(
      observed_delta = mean(unit_table$delta),
      ci_low = NA_real_,
      ci_high = NA_real_,
      probability_a_better = NA_real_,
      n_environments = nrow(environment_summary),
      n_units = nrow(unit_table)
    ))
  }

  set.seed(seed)
  bootstrap_delta <- replicate(replicates, {
    sampled_rows <- sample.int(
      nrow(environment_summary),
      nrow(environment_summary),
      replace = TRUE
    )
    sum(environment_summary$delta_sum[sampled_rows]) /
      sum(environment_summary$n_units[sampled_rows])
  })

  data.table(
    observed_delta = mean(unit_table$delta),
    ci_low = as.numeric(quantile(bootstrap_delta, 0.025)),
    ci_high = as.numeric(quantile(bootstrap_delta, 0.975)),
    probability_a_better = mean(bootstrap_delta < 0),
    n_environments = nrow(environment_summary),
    n_units = nrow(unit_table)
  )
}

primary_units <- environment_stage_metrics[
  feature_set == PRIMARY_FEATURE_SET
]
learner_ids <- sort(unique(primary_units$learner_id))

if (length(learner_ids) >= 2L) {
  learner_pairs <- combn(
    learner_ids,
    2L,
    simplify = FALSE
  )

  paired_learner_bootstrap <- rbindlist(lapply(
    seq_along(learner_pairs),
    function(pair_index) {
      pair <- learner_pairs[[pair_index]]

      rbindlist(lapply(
        sort(unique(primary_units$scheme)),
        function(scheme_id) {
          learner_a_units <- primary_units[
            scheme == scheme_id & learner_id == pair[1L],
            .(
              env_id,
              BBCH_stage,
              mae_a = mae
            )
          ]
          learner_b_units <- primary_units[
            scheme == scheme_id & learner_id == pair[2L],
            .(
              env_id,
              BBCH_stage,
              mae_b = mae
            )
          ]

          paired_units <- merge(
            learner_a_units,
            learner_b_units,
            by = c("env_id", "BBCH_stage"),
            all = FALSE
          )
          paired_units[, delta := mae_a - mae_b]

          result <- cluster_bootstrap_delta(
            unit_table = paired_units,
            replicates = BOOTSTRAP_REPLICATES,
            seed = SEED + 1000L * pair_index +
              match(scheme_id, sort(unique(primary_units$scheme)))
          )
          result[, learner_a := pair[1L]]
          result[, learner_b := pair[2L]]
          result[, scheme := scheme_id]
          result
        }
      ))
    }
  ))

  setcolorder(
    paired_learner_bootstrap,
    c(
      "scheme", "learner_a", "learner_b",
      "observed_delta", "ci_low", "ci_high",
      "probability_a_better", "n_environments", "n_units"
    )
  )
  setorder(
    paired_learner_bootstrap,
    scheme,
    learner_a,
    learner_b
  )
  write_table(
    paired_learner_bootstrap,
    "paired_learner_bootstrap_primary_cumulative.csv"
  )
}

# ------------------------------------------------------------------------------
# 16. CONFIGURATION AND PROVENANCE
# ------------------------------------------------------------------------------

configuration <- list(
  run_id = RUN_ID,
  input_file = normalizePath(INPUT_FILE),
  input_md5 = unname(tools::md5sum(INPUT_FILE)),
  seed = SEED,
  target = TARGET,
  site_group_column = site_column,
  environment_column = environment_column,
  year_column = year_column,
  date_column = date_column,
  feature_sets_to_run = FEATURE_SETS_TO_RUN,
  primary_feature_set = PRIMARY_FEATURE_SET,
  feature_sets = feature_sets,
  gam_formulas = lapply(
    gam_formulas,
    function(formula) {
      paste(deparse(formula, width.cutoff = 500L), collapse = " ")
    }
  ),
  gam_engine = "gam",
  gam_family = GAM_FAMILY,
  gam_method = "REML",
  gam_select = TRUE,
  cv_schemes = CV_SCHEMES,
  temporal_blocks = TEMPORAL_BLOCKS,
  spatial_folds = SPATIAL_FOLDS,
  spatial_block_size_m = SPATIAL_BLOCK_SIZE_M,
  spatial_selection = SPATIAL_SELECTION,
  inner_group_folds = INNER_GROUP_FOLDS,
  tuning_measure = "regr.mae",
  tuner = "mbo",
  tuning_evals = list(
    random_forest = TUNING_EVALS_RF,
    mars = TUNING_EVALS_MARS,
    xgboost = TUNING_EVALS_XGB
  ),
  rf_num_trees = RF_NUM_TREES,
  bootstrap_replicates = BOOTSTRAP_REPLICATES,
  extreme_prediction_limits = c(
    EXTREME_PREDICTION_LOWER,
    EXTREME_PREDICTION_UPPER
  ),
  checkpoint_version = CHECKPOINT_VERSION,
  results_mode = RESULTS_MODE,
  workers = workers,
  n_input_rows = nrow(raw_data),
  n_common_cohort_rows = nrow(cohort),
  n_environments = uniqueN(cohort$.env_id),
  n_sites = uniqueN(cohort$.site_id),
  n_years = uniqueN(cohort$.season_year)
)

saveRDS(
  configuration,
  file.path(OUTPUT_DIR, "benchmark_configuration.rds")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "session_info.txt")
)

message("Learner-selection benchmark finished.")
message("Results written to: ", OUTPUT_DIR)
message(
  "Primary learner table: ",
  file.path(
    OUTPUT_DIR,
    "learner_selection_primary_cumulative.csv"
  )
)

# ==============================================================================
# Transfer benchmark for the selected winter-wheat phenology model
# ==============================================================================
# Assess the selected learner and cumulative-weather representation under five
# transfer schemes. Fit a calibrated thermal-time baseline in each outer
# training split. Report scheme-specific errors separately. Preserve fold-level
# fitted settings and out-of-fold predictions for the subsequent TreeSHAP step.
# Macro-RMSE is sqrt(mean(environment-stage MSE)).
# Run from the repository root or set PHENOLOGY_PROJECT_DIR to its path.

# ------------------------------------------------------------------------------
# 1. USER CONFIGURATION
# ------------------------------------------------------------------------------

PROJECT_DIR <- Sys.getenv("PHENOLOGY_PROJECT_DIR", unset = getwd())
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = TRUE)
INPUT_FILE <- Sys.getenv(
  "PHENOLOGY_GENERALISATION_INPUT",
  unset = file.path(PROJECT_DIR, "Data", "pheno_df_CV.csv")
)

# Change RUN_ID or CHECKPOINT_VERSION whenever the cohort, predictors, folds,
# learner definition, tuning space or tuning resampling changes.
SELECTED_LEARNER <- tolower(trimws(Sys.getenv(
  "PHENOLOGY_GENERALISATION_LEARNER",
  unset = "xgboost"
)))

allowed_selected_learners <- c("random_forest", "xgboost")
if (!SELECTED_LEARNER %in% allowed_selected_learners) {
  stop(
    "PHENOLOGY_GENERALISATION_LEARNER must be random_forest or xgboost; got: ",
    SELECTED_LEARNER
  )
}

RUN_ID <- Sys.getenv(
  "PHENOLOGY_GENERALISATION_RUN_ID",
  unset = paste0("multischeme_", SELECTED_LEARNER, "_cumulative_weather_paper")
)

OUTPUT_DIR <- file.path(PROJECT_DIR, "benchmark_results", RUN_ID)
CHECKPOINT_DIR <- file.path(OUTPUT_DIR, "checkpoints")

SEED <- 123L
TARGET <- "BBCH"
FEATURE_SET_ID <- "cumulative_weather"
FEATURES <- c(
  "cumDTU_beta",
  "sowing_DOY",
  "photoperiod_civil",
  "cumRR",
  "cumRad",
  "cumCDD",
  "cumHDD",
  "paw_fc33_0_100_mm"
)

GDD_COLUMN <- "cumDTU_beta"
RUN_GDD_BASELINE <- TRUE

SITE_GROUP_COLUMN <- Sys.getenv("PHENOLOGY_SITE_COLUMN", unset = "Standort")
SITE_COLUMN_CANDIDATES <- unique(c(
  SITE_GROUP_COLUMN, "Standort", "Location"
))

ZONE_COLUMN <- Sys.getenv("PHENOLOGY_ZONE_COLUMN", unset = "env_zone")
ZONE_COLUMN_CANDIDATES <- unique(c(ZONE_COLUMN, "env_zone"))

default_cv_schemes <- c(
  "spatial_block",
  "leave_one_year_out",
  "expanding_origin",
  "environmental_zone",
  "joint_space_time"
)

# Optional comma-separated override, for example:
# PHENOLOGY_GENERALISATION_SCHEMES="spatial_block,expanding_origin"
scheme_override <- trimws(Sys.getenv(
  "PHENOLOGY_GENERALISATION_SCHEMES",
  unset = ""
))
CV_SCHEMES <- if (nzchar(scheme_override)) {
  unique(trimws(strsplit(scheme_override, ",", fixed = TRUE)[[1L]]))
} else {
  default_cv_schemes
}
CV_SCHEMES <- CV_SCHEMES[nzchar(CV_SCHEMES)]

unknown_schemes <- setdiff(CV_SCHEMES, default_cv_schemes)
if (length(unknown_schemes)) {
  stop("Unknown CV scheme(s): ", paste(unknown_schemes, collapse = ", "))
}
if (!length(CV_SCHEMES)) stop("No CV schemes were requested.")

# Outer spatial and joint space-time structure.
SPATIAL_FOLDS <- 5L
SPATIAL_BLOCK_SIZE_M <- 119946.2
SPATIAL_SELECTION <- "random"
JOINT_TEMPORAL_BLOCKS <- 5L

# Expanding-origin starts after this many observed season years. Each later
# observed year is tested once using only earlier years for training.
EXPANDING_MIN_TRAINING_YEARS <- 10L

# Scheme-matched inner tuning. Fallback is grouped site CV if a requested
# scheme cannot form at least two valid inner splits in an outer training set.
INNER_GROUP_FOLDS <- 5L
INNER_TEMPORAL_BLOCKS <- 5L
INNER_ZONE_FOLDS <- 5L
INNER_FORWARD_SPLITS <- 5L
INNER_FORWARD_MIN_TRAINING_YEARS <- 5L
INNER_JOINT_SPLITS <- 5L

# Same search budgets and learner definitions used for learner selection.
TUNING_EVALS_RF <- 40L
TUNING_EVALS_XGB <- 60L
RF_NUM_TREES <- 1000L

WORKERS_REQUESTED <- 10L
RESUME_FROM_CHECKPOINTS <- TRUE
SAVE_TUNING_ARCHIVES <- TRUE
CHECKPOINT_VERSION <- 1L

BOOTSTRAP_REPLICATES <- 2000L

# Predictions are never clipped. These limits only flag numerical failure.
EXTREME_PREDICTION_LOWER <- -100
EXTREME_PREDICTION_UPPER <- 200

# Minimum sizes protect against malformed outer splits. Joint space-time cells
# with no observations are recorded and omitted before jobs are created.
MIN_OUTER_TRAIN_ROWS <- 100L
MIN_OUTER_TEST_ROWS <- 1L

# ------------------------------------------------------------------------------
# 2. PACKAGES AND REPRODUCIBILITY
# ------------------------------------------------------------------------------

required_packages <- c(
  "data.table", "R6", "mlr3", "mlr3learners", "mlr3tuning",
  "mlr3mbo", "paradox", "future", "future.apply", "parallelly",
  "ranger", "xgboost", "mlr3spatiotempcv", "sf", "blockCV"
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
  library(mlr3tuning)
  library(mlr3mbo)
  library(paradox)
  library(future)
  library(future.apply)
  library(parallelly)
  library(mlr3spatiotempcv)
  library(xgboost)
})

RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)
data.table::setDTthreads(1L)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

workers <- max(
  1L,
  min(WORKERS_REQUESTED, parallelly::availableCores(omit = 1L))
)

# ------------------------------------------------------------------------------
# 3. UTILITIES
# ------------------------------------------------------------------------------

safe_id <- function(x) {
  gsub("[^A-Za-z0-9_-]+", "_", as.character(x))
}

parse_boolean <- function(x, default = FALSE) {
  if (length(x) == 0L || is.na(x) || !nzchar(trimws(x))) return(default)
  normalized <- tolower(trimws(x))
  if (normalized %in% c("1", "true", "yes", "y")) return(TRUE)
  if (normalized %in% c("0", "false", "no", "n")) return(FALSE)
  stop("Cannot parse logical value: ", x)
}

# Permit server-side overrides without editing the script.
RUN_GDD_BASELINE <- parse_boolean(
  Sys.getenv("PHENOLOGY_RUN_GDD_BASELINE", unset = as.character(RUN_GDD_BASELINE)),
  default = RUN_GDD_BASELINE
)

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

write_table <- function(x, filename) {
  fwrite(x, file.path(OUTPUT_DIR, filename), na = "NA")
}

atomic_save_rds <- function(object, path, compress = FALSE) {
  temporary_path <- paste0(path, ".tmp_", Sys.getpid())
  saveRDS(object, temporary_path, compress = compress)
  if (!file.rename(temporary_path, path)) {
    unlink(temporary_path)
    stop("Could not atomically write: ", path)
  }
}

capture_condition <- function(store, condition) {
  unique(c(store, conditionMessage(condition)))
}

# mlr3 learner dictionaries are populated through package namespace load hooks.
# The selected learner is constructed inside multisession workers, which start
# as clean R sessions. Explicitly load mlr3learners in the current process before
# asking mlr3::lrn() for a recommended learner; relying on the main session's
# search path is insufficient for side-effect registrations in future workers.
ensure_recommended_learner_registered <- function(learner_key) {
  if (!requireNamespace("mlr3learners", quietly = TRUE)) {
    stop(
      "Package 'mlr3learners' is required to register ", learner_key, "."
    )
  }

  loadNamespace("mlr3learners")

  if (!learner_key %in% mlr3::mlr_learners$keys()) {
    package_version_or_na <- function(package) {
      if (requireNamespace(package, quietly = TRUE)) {
        as.character(utils::packageVersion(package))
      } else {
        "not installed"
      }
    }

    matching_keys <- grep(
      "xgboost|ranger",
      mlr3::mlr_learners$keys(),
      value = TRUE
    )

    stop(
      "Learner '", learner_key,
      "' is still absent after explicitly loading mlr3learners in this R ",
      "process. Installed versions: mlr3=",
      package_version_or_na("mlr3"),
      ", mlr3learners=",
      package_version_or_na("mlr3learners"),
      ", xgboost=",
      package_version_or_na("xgboost"),
      ", ranger=",
      package_version_or_na("ranger"),
      ". Registered tree learner keys: ",
      if (length(matching_keys)) paste(matching_keys, collapse = ", ") else "none",
      ". Restore/update the renv packages before rerunning."
    )
  }

  invisible(TRUE)
}

selected_mlr3_learner_key <- switch(
  SELECTED_LEARNER,
  random_forest = "regr.ranger",
  xgboost = "regr.xgboost"
)
ensure_recommended_learner_registered(selected_mlr3_learner_key)

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
    stop("Prediction does not expose row_ids, truth and response.")
  }
  if (is.data.frame(response)) response <- as.matrix(response)
  if (is.matrix(response) && ncol(response) != 1L) {
    stop("Expected one response column; received ", ncol(response), ".")
  }

  answer <- data.table(
    row_ids = as.integer(row_ids),
    truth = as.numeric(truth),
    response = as.numeric(response)
  )
  if (
    anyNA(answer$row_ids) ||
      any(!is.finite(answer$truth)) ||
      any(!is.finite(answer$response))
  ) {
    stop("Prediction contains missing or non-finite values.")
  }
  answer
}

make_contiguous_year_blocks <- function(year, k) {
  unique_years <- sort(unique(as.integer(year)))
  k <- min(as.integer(k), length(unique_years))
  if (k < 2L) stop("At least two temporal blocks are required.")

  year_block <- ceiling(seq_along(unique_years) * k / length(unique_years))
  year_block <- pmin(year_block, k)
  lookup <- setNames(year_block, unique_years)
  unname(lookup[as.character(as.integer(year))])
}

make_balanced_group_folds <- function(group, k, seed) {
  group <- as.character(group)
  if (anyNA(group) || any(!nzchar(group))) {
    stop("Balanced grouped folds received missing/empty group IDs.")
  }

  group_sizes <- data.table(group = group)[, .N, by = group]
  k <- min(as.integer(k), nrow(group_sizes))
  if (k < 2L) stop("At least two groups are required for grouped folds.")

  set.seed(seed)
  group_sizes[, tie_break := runif(.N)]
  setorder(group_sizes, -N, tie_break)

  fold_load <- numeric(k)
  group_sizes[, fold := NA_integer_]
  for (i in seq_len(nrow(group_sizes))) {
    candidate_folds <- which(fold_load == min(fold_load))
    selected_fold <- candidate_folds[sample.int(length(candidate_folds), 1L)]
    group_sizes$fold[i] <- selected_fold
    fold_load[selected_fold] <- fold_load[selected_fold] + group_sizes$N[i]
  }

  lookup <- setNames(group_sizes$fold, group_sizes$group)
  unname(lookup[group])
}

resampling_from_labels <- function(task, labels) {
  row_ids <- task$row_ids
  labels <- as.character(labels)
  if (length(labels) != length(row_ids) || anyNA(labels)) {
    stop("Inner fold labels do not match the inner task.")
  }

  levels_in_order <- unique(labels)
  train_sets <- lapply(levels_in_order, function(level) {
    row_ids[labels != level]
  })
  test_sets <- lapply(levels_in_order, function(level) {
    row_ids[labels == level]
  })

  custom <- rsmp("custom")
  custom$instantiate(task, train_sets = train_sets, test_sets = test_sets)
  custom
}

resampling_from_sets <- function(task, train_sets, test_sets) {
  if (length(train_sets) != length(test_sets) || length(train_sets) < 2L) {
    stop("Custom inner resampling requires at least two paired splits.")
  }
  custom <- rsmp("custom")
  custom$instantiate(task, train_sets = train_sets, test_sets = test_sets)
  custom
}

summarize_inner_resampling <- function(resampling, strategy, cohort) {
  rbindlist(lapply(seq_len(resampling$iters), function(iteration) {
    train_ids <- resampling$train_set(iteration)
    test_ids <- resampling$test_set(iteration)
    data.table(
      inner_iteration = iteration,
      inner_strategy = strategy,
      n_train = length(train_ids),
      n_test = length(test_ids),
      n_train_sites = uniqueN(cohort$.site_id[train_ids]),
      n_test_sites = uniqueN(cohort$.site_id[test_ids]),
      train_year_min = min(cohort$.season_year[train_ids]),
      train_year_max = max(cohort$.season_year[train_ids]),
      test_year_min = min(cohort$.season_year[test_ids]),
      test_year_max = max(cohort$.season_year[test_ids])
    )
  }))
}

# ------------------------------------------------------------------------------
# 4. LOAD DATA AND FREEZE THE ANALYSIS COHORT
# ------------------------------------------------------------------------------

if (!file.exists(INPUT_FILE)) stop("Input file does not exist: ", INPUT_FILE)

message("Reading: ", INPUT_FILE)
raw_data <- fread(
  INPUT_FILE,
  na.strings = c("", "NA", "NaN", "Inf", "-Inf"),
  showProgress = TRUE
)
raw_data[, .source_row_id := .I]

site_column <- resolve_column(
  raw_data, SITE_COLUMN_CANDIDATES, "physical-site grouping column"
)
environment_column <- resolve_column(
  raw_data, c("env_id", "traj_id"), "site-season environment ID"
)
year_column <- resolve_column(
  raw_data, c("SeasonYear", "season_year", "year"), "season year"
)
date_column <- resolve_column(
  raw_data, c("Date", "phen_date", "date"), "observation date"
)
longitude_column <- resolve_column(
  raw_data, c("lon", "longitude"), "longitude"
)
latitude_column <- resolve_column(
  raw_data, c("lat", "latitude"), "latitude"
)

needs_zone <- "environmental_zone" %in% CV_SCHEMES
zone_column <- if (needs_zone) {
  resolve_column(raw_data, ZONE_COLUMN_CANDIDATES, "environmental-zone column")
} else {
  zone_hits <- ZONE_COLUMN_CANDIDATES[ZONE_COLUMN_CANDIDATES %in% names(raw_data)]
  if (length(zone_hits)) zone_hits[1L] else NA_character_
}

missing_model_columns <- setdiff(c(TARGET, FEATURES), names(raw_data))
if (length(missing_model_columns)) {
  stop(
    "Frozen-pipeline columns are absent: ",
    paste(missing_model_columns, collapse = ", ")
  )
}

non_numeric_features <- FEATURES[
  !vapply(raw_data[, ..FEATURES], is.numeric, logical(1))
]
if (length(non_numeric_features)) {
  stop(
    "All frozen predictors must be numeric. Non-numeric: ",
    paste(non_numeric_features, collapse = ", ")
  )
}

raw_data[, .site_id := trimws(as.character(get(site_column)))]
raw_data[, .env_id := trimws(as.character(get(environment_column)))]
raw_data[, .season_year := as.integer(get(year_column))]
raw_data[, .date := as.IDate(get(date_column))]
raw_data[, .lon := as.numeric(get(longitude_column))]
raw_data[, .lat := as.numeric(get(latitude_column))]
raw_data[, .env_zone := if (!is.na(zone_column)) {
  trimws(as.character(get(zone_column)))
} else {
  NA_character_
}]

raw_data[!nzchar(.site_id), .site_id := NA_character_]
raw_data[!nzchar(.env_id), .env_id := NA_character_]
raw_data[!is.na(.env_zone) & !nzchar(.env_zone), .env_zone := NA_character_]

required_common_columns <- unique(c(
  TARGET, FEATURES, ".site_id", ".env_id", ".season_year", ".date",
  ".lon", ".lat", if (needs_zone) ".env_zone"
))

missingness_audit <- rbindlist(lapply(required_common_columns, function(column) {
  value <- raw_data[[column]]
  data.table(
    column = column,
    n_missing = sum(is.na(value)),
    proportion_missing = mean(is.na(value))
  )
}))

complete_mask <- complete.cases(raw_data[, ..required_common_columns])
finite_columns <- c(TARGET, FEATURES, ".season_year", ".lon", ".lat")
for (column in unique(finite_columns)) {
  complete_mask <- complete_mask & is.finite(raw_data[[column]])
}

cohort <- copy(raw_data[complete_mask])
cohort[, .obs_id := .I]

cohort_flow <- data.table(
  step = c("input", "complete_frozen_pipeline_cohort", "removed"),
  n_rows = c(nrow(raw_data), nrow(cohort), nrow(raw_data) - nrow(cohort)),
  proportion_of_input = c(
    1,
    nrow(cohort) / nrow(raw_data),
    (nrow(raw_data) - nrow(cohort)) / nrow(raw_data)
  )
)

write_table(missingness_audit, "missingness_required_columns.csv")
write_table(cohort_flow, "cohort_flow.csv")

if (nrow(cohort) < MIN_OUTER_TRAIN_ROWS + MIN_OUTER_TEST_ROWS) {
  stop("The frozen-pipeline cohort is too small.")
}
if (uniqueN(cohort$.site_id) < INNER_GROUP_FOLDS) {
  stop("Too few physical sites for grouped inner-CV fallback.")
}

# ------------------------------------------------------------------------------
# 5. INPUT QUALITY CONTROL
# ------------------------------------------------------------------------------

if (any(cohort[[TARGET]] < 0 | cohort[[TARGET]] > 99)) {
  stop("BBCH contains values outside 0-99.")
}
if (any(cohort$sowing_DOY < 1 | cohort$sowing_DOY > 366)) {
  stop("sowing_DOY contains values outside 1-366.")
}
if (any(cohort$Daylength <= 0 | cohort$Daylength > 24)) {
  stop("Daylength contains values outside (0, 24] hours.")
}
if (any(cohort[[GDD_COLUMN]] < 0)) {
  stop(GDD_COLUMN, " contains negative values.")
}

environment_mapping_qc <- cohort[, .(
  n_sites = uniqueN(.site_id),
  n_years = uniqueN(.season_year),
  n_zones = uniqueN(.env_zone, na.rm = TRUE)
), by = .env_id]

invalid_environment_ids <- environment_mapping_qc[
  n_sites != 1L | n_years != 1L | (needs_zone & n_zones != 1L)
]
write_table(environment_mapping_qc, "environment_mapping_qc.csv")
if (nrow(invalid_environment_ids)) {
  write_table(invalid_environment_ids, "invalid_environment_ids.csv")
  stop(
    "At least one environment maps to multiple sites, years or zones. ",
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
    sum(diff(get(GDD_COLUMN)) < -1e-8)
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

feature_qc <- rbindlist(lapply(FEATURES, function(column) {
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
write_table(feature_qc, "feature_qc_frozen_cohort.csv")

constant_features <- feature_qc[n_unique < 2L, feature]
if (length(constant_features)) {
  stop("Constant frozen predictors: ", paste(constant_features, collapse = ", "))
}

write_table(
  data.table(
    feature_set = FEATURE_SET_ID,
    feature = FEATURES,
    role = "frozen_selected_representation"
  ),
  "feature_set_definition.csv"
)

# ------------------------------------------------------------------------------
# 6. MLR3 TASK AND REUSABLE SPATIAL/TEMPORAL BLOCKS
# ------------------------------------------------------------------------------

backend_columns <- unique(c(
  ".obs_id", ".env_id", ".site_id", TARGET, ".lon", ".lat", FEATURES
))
backend <- cohort[, ..backend_columns]

task <- TaskRegrST$new(
  id = FEATURE_SET_ID,
  backend = backend,
  target = TARGET,
  coordinate_names = c(".lon", ".lat"),
  crs = "EPSG:4326"
)
task$select(FEATURES)

needs_spatial_blocks <- any(c("spatial_block", "joint_space_time") %in% CV_SCHEMES)
if (needs_spatial_blocks) {
  set.seed(SEED + 20L)
  spatial_resampling <- rsmp(
    "spcv_block",
    range = as.integer(SPATIAL_BLOCK_SIZE_M),
    folds = as.integer(SPATIAL_FOLDS),
    selection = SPATIAL_SELECTION,
    hexagon = FALSE
  )
  spatial_resampling$instantiate(task)

  spatial_instance <- as.data.table(spatial_resampling$instance)
  row_id_column <- intersect(c("row_id", "row_ids"), names(spatial_instance))[1L]
  if (!length(row_id_column) || is.na(row_id_column)) {
    stop("Could not identify row IDs in the spatial resampling instance.")
  }

  cohort[, .spatial_fold := as.character(spatial_instance$fold[
    match(.obs_id, spatial_instance[[row_id_column]])
  ])]
  if (anyNA(cohort$.spatial_fold)) {
    stop("Spatial blocking did not assign every cohort row.")
  }

  site_spatial_check <- cohort[, .(
    n_spatial_folds = uniqueN(.spatial_fold)
  ), by = .site_id]
  if (any(site_spatial_check$n_spatial_folds != 1L)) {
    write_table(
      site_spatial_check[n_spatial_folds != 1L],
      "spatial_sites_split_across_folds.csv"
    )
    stop("Spatial blocking split at least one physical site.")
  }
} else {
  cohort[, .spatial_fold := NA_character_]
}

if ("joint_space_time" %in% CV_SCHEMES) {
  cohort[, .joint_temporal_block := as.character(make_contiguous_year_blocks(
    .season_year,
    JOINT_TEMPORAL_BLOCKS
  ))]
} else {
  cohort[, .joint_temporal_block := NA_character_]
}

# ------------------------------------------------------------------------------
# 7. EXPLICIT OUTER SPLIT SPECIFICATIONS
# ------------------------------------------------------------------------------

split_parts <- list()

if ("spatial_block" %in% CV_SCHEMES) {
  spatial_levels <- sort(unique(cohort$.spatial_fold))
  split_parts$spatial_block <- data.table(
    scheme = "spatial_block",
    fold_id = paste0("space_", spatial_levels),
    holdout_spatial_fold = spatial_levels,
    holdout_year = NA_integer_,
    holdout_zone = NA_character_,
    holdout_temporal_block = NA_character_
  )
}

if ("leave_one_year_out" %in% CV_SCHEMES) {
  years <- sort(unique(cohort$.season_year))
  if (length(years) < 2L) stop("Leave-one-year-out requires at least two years.")
  split_parts$leave_one_year_out <- data.table(
    scheme = "leave_one_year_out",
    fold_id = paste0("year_", years),
    holdout_spatial_fold = NA_character_,
    holdout_year = years,
    holdout_zone = NA_character_,
    holdout_temporal_block = NA_character_
  )
}

if ("expanding_origin" %in% CV_SCHEMES) {
  years <- sort(unique(cohort$.season_year))
  if (length(years) <= EXPANDING_MIN_TRAINING_YEARS) {
    stop(
      "Expanding-origin CV requires more than ",
      EXPANDING_MIN_TRAINING_YEARS, " observed season years."
    )
  }
  forecast_years <- years[seq_along(years) > EXPANDING_MIN_TRAINING_YEARS]
  split_parts$expanding_origin <- data.table(
    scheme = "expanding_origin",
    fold_id = paste0("forecast_", forecast_years),
    holdout_spatial_fold = NA_character_,
    holdout_year = forecast_years,
    holdout_zone = NA_character_,
    holdout_temporal_block = NA_character_
  )
}

if ("environmental_zone" %in% CV_SCHEMES) {
  zones <- sort(unique(cohort$.env_zone))
  if (length(zones) < 2L) {
    stop("Environmental-zone CV requires at least two zones.")
  }
  split_parts$environmental_zone <- data.table(
    scheme = "environmental_zone",
    fold_id = paste0("zone_", zones),
    holdout_spatial_fold = NA_character_,
    holdout_year = NA_integer_,
    holdout_zone = zones,
    holdout_temporal_block = NA_character_
  )
}

if ("joint_space_time" %in% CV_SCHEMES) {
  joint_grid <- CJ(
    holdout_spatial_fold = sort(unique(cohort$.spatial_fold)),
    holdout_temporal_block = sort(unique(cohort$.joint_temporal_block)),
    unique = TRUE
  )
  joint_grid[, scheme := "joint_space_time"]
  joint_grid[, fold_id := paste0(
    "space_", holdout_spatial_fold,
    "__time_", holdout_temporal_block
  )]
  joint_grid[, holdout_year := NA_integer_]
  joint_grid[, holdout_zone := NA_character_]
  setcolorder(
    joint_grid,
    c(
      "scheme", "fold_id", "holdout_spatial_fold", "holdout_year",
      "holdout_zone", "holdout_temporal_block"
    )
  )
  split_parts$joint_space_time <- joint_grid
}

outer_split_specs <- rbindlist(split_parts, use.names = TRUE, fill = TRUE)
outer_split_specs[, split_id := .I]
setcolorder(outer_split_specs, c("split_id", setdiff(names(outer_split_specs), "split_id")))

resolve_outer_ids <- function(spec) {
  scheme <- as.character(spec$scheme)

  if (identical(scheme, "spatial_block")) {
    test_mask <- cohort$.spatial_fold == as.character(spec$holdout_spatial_fold)
    train_mask <- !test_mask
  } else if (identical(scheme, "leave_one_year_out")) {
    test_mask <- cohort$.season_year == as.integer(spec$holdout_year)
    train_mask <- !test_mask
  } else if (identical(scheme, "expanding_origin")) {
    test_mask <- cohort$.season_year == as.integer(spec$holdout_year)
    train_mask <- cohort$.season_year < as.integer(spec$holdout_year)
  } else if (identical(scheme, "environmental_zone")) {
    test_mask <- cohort$.env_zone == as.character(spec$holdout_zone)
    train_mask <- !test_mask
  } else if (identical(scheme, "joint_space_time")) {
    heldout_space <- as.character(spec$holdout_spatial_fold)
    heldout_time <- as.character(spec$holdout_temporal_block)
    test_mask <- (
      cohort$.spatial_fold == heldout_space &
        cohort$.joint_temporal_block == heldout_time
    )
    train_mask <- (
      cohort$.spatial_fold != heldout_space &
        cohort$.joint_temporal_block != heldout_time
    )
  } else {
    stop("Unknown outer scheme: ", scheme)
  }

  list(
    train_ids = cohort$.obs_id[train_mask],
    test_ids = cohort$.obs_id[test_mask],
    excluded_ids = cohort$.obs_id[!train_mask & !test_mask]
  )
}

outer_split_sizes <- rbindlist(lapply(
  seq_len(nrow(outer_split_specs)),
  function(i) {
    spec <- outer_split_specs[i]
    ids <- resolve_outer_ids(spec)
    data.table(
      split_id = spec$split_id,
      scheme = spec$scheme,
      fold_id = spec$fold_id,
      n_train = length(ids$train_ids),
      n_test = length(ids$test_ids),
      n_excluded = length(ids$excluded_ids),
      n_train_environments = uniqueN(cohort$.env_id[ids$train_ids]),
      n_test_environments = uniqueN(cohort$.env_id[ids$test_ids]),
      n_train_sites = uniqueN(cohort$.site_id[ids$train_ids]),
      n_test_sites = uniqueN(cohort$.site_id[ids$test_ids]),
      train_year_min = if (length(ids$train_ids)) {
        min(cohort$.season_year[ids$train_ids])
      } else NA_integer_,
      train_year_max = if (length(ids$train_ids)) {
        max(cohort$.season_year[ids$train_ids])
      } else NA_integer_,
      test_year_min = if (length(ids$test_ids)) {
        min(cohort$.season_year[ids$test_ids])
      } else NA_integer_,
      test_year_max = if (length(ids$test_ids)) {
        max(cohort$.season_year[ids$test_ids])
      } else NA_integer_
    )
  }
))

valid_split_ids <- outer_split_sizes[
  n_train >= MIN_OUTER_TRAIN_ROWS & n_test >= MIN_OUTER_TEST_ROWS,
  split_id
]
omitted_outer_splits <- outer_split_sizes[!split_id %in% valid_split_ids]
outer_split_specs <- outer_split_specs[split_id %in% valid_split_ids]
outer_split_sizes <- outer_split_sizes[split_id %in% valid_split_ids]

write_table(outer_split_specs, "outer_split_specs.csv")
write_table(outer_split_sizes, "outer_split_sizes.csv")
write_table(omitted_outer_splits, "omitted_outer_splits.csv")

if (!nrow(outer_split_specs)) stop("No valid outer splits remain.")
missing_valid_schemes <- setdiff(CV_SCHEMES, unique(outer_split_specs$scheme))
if (length(missing_valid_schemes)) {
  stop(
    "No valid outer splits for scheme(s): ",
    paste(missing_valid_schemes, collapse = ", ")
  )
}

# Verify all outer leakage constraints before fitting.
outer_leakage_qc <- rbindlist(lapply(seq_len(nrow(outer_split_specs)), function(i) {
  spec <- outer_split_specs[i]
  ids <- resolve_outer_ids(spec)

  environment_overlap <- intersect(
    unique(cohort$.env_id[ids$train_ids]),
    unique(cohort$.env_id[ids$test_ids])
  )
  site_overlap <- intersect(
    unique(cohort$.site_id[ids$train_ids]),
    unique(cohort$.site_id[ids$test_ids])
  )
  year_overlap <- intersect(
    unique(cohort$.season_year[ids$train_ids]),
    unique(cohort$.season_year[ids$test_ids])
  )
  zone_overlap <- intersect(
    unique(cohort$.env_zone[ids$train_ids]),
    unique(cohort$.env_zone[ids$test_ids])
  )

  data.table(
    split_id = spec$split_id,
    scheme = spec$scheme,
    fold_id = spec$fold_id,
    n_environment_overlap = length(environment_overlap),
    n_site_overlap = length(site_overlap),
    n_year_overlap = length(year_overlap),
    n_zone_overlap = length(zone_overlap)
  )
}))
write_table(outer_leakage_qc, "outer_leakage_qc.csv")

if (any(outer_leakage_qc$n_environment_overlap > 0L)) {
  stop("At least one outer split has environment leakage.")
}
if (any(
  outer_leakage_qc$scheme == "spatial_block" &
    outer_leakage_qc$n_site_overlap > 0L
)) {
  stop("Spatial-block CV has physical-site leakage.")
}
if (any(
  outer_leakage_qc$scheme %in% c("leave_one_year_out", "expanding_origin") &
    outer_leakage_qc$n_year_overlap > 0L
)) {
  stop("A temporal outer split has season-year leakage.")
}
if (any(
  outer_leakage_qc$scheme == "environmental_zone" &
    outer_leakage_qc$n_zone_overlap > 0L
)) {
  stop("Environmental-zone CV has zone leakage.")
}
if (any(
  outer_leakage_qc$scheme == "joint_space_time" &
    (
      outer_leakage_qc$n_site_overlap > 0L |
        outer_leakage_qc$n_year_overlap > 0L
    )
)) {
  stop("Joint space-time CV does not fully exclude test sites and years.")
}

# One-row-per-observation manifest for reconstruction and future SHAP. The
# outer_split_specs file defines training rules for expanding-origin and joint
# space-time folds.
fold_manifest <- cohort[, .(
  obs_id = .obs_id,
  source_row_id = .source_row_id,
  env_id = .env_id,
  site_id = .site_id,
  season_year = .season_year,
  env_zone = .env_zone,
  Date = .date,
  BBCH = get(TARGET),
  lon = .lon,
  lat = .lat,
  spatial_block = .spatial_fold,
  joint_temporal_block = .joint_temporal_block
)]

if ("spatial_block" %in% CV_SCHEMES) {
  fold_manifest[, fold_spatial_block := paste0("space_", spatial_block)]
}
if ("leave_one_year_out" %in% CV_SCHEMES) {
  fold_manifest[, fold_leave_one_year_out := paste0("year_", season_year)]
}
if ("expanding_origin" %in% CV_SCHEMES) {
  first_forecast_year <- min(
    outer_split_specs[scheme == "expanding_origin", holdout_year]
  )
  fold_manifest[, fold_expanding_origin := ifelse(
    season_year >= first_forecast_year,
    paste0("forecast_", season_year),
    NA_character_
  )]
}
if ("environmental_zone" %in% CV_SCHEMES) {
  fold_manifest[, fold_environmental_zone := paste0("zone_", env_zone)]
}
if ("joint_space_time" %in% CV_SCHEMES) {
  fold_manifest[, fold_joint_space_time := paste0(
    "space_", spatial_block, "__time_", joint_temporal_block
  )]
}
write_table(fold_manifest, "outer_fold_manifest.csv")

split_coverage <- rbindlist(lapply(unique(outer_split_specs$scheme), function(scheme_id) {
  scheme_specs <- outer_split_specs[scheme == scheme_id]
  tested_ids <- unique(unlist(lapply(seq_len(nrow(scheme_specs)), function(i) {
    resolve_outer_ids(scheme_specs[i])$test_ids
  }), use.names = FALSE))
  data.table(
    scheme = scheme_id,
    n_outer_splits = nrow(scheme_specs),
    n_unique_test_rows = length(tested_ids),
    proportion_cohort_tested = length(tested_ids) / nrow(cohort),
    n_unique_test_environments = uniqueN(cohort$.env_id[tested_ids]),
    n_unique_test_sites = uniqueN(cohort$.site_id[tested_ids]),
    n_unique_test_years = uniqueN(cohort$.season_year[tested_ids])
  )
}))
write_table(split_coverage, "outer_split_coverage.csv")

# ------------------------------------------------------------------------------
# 8. TRAINING-ONLY GDD BASELINE
# ------------------------------------------------------------------------------

LearnerRegrGDDThreshold <- R6Class(
  "LearnerRegrGDDThreshold",
  inherit = LearnerRegr,
  public = list(
    gdd_column = NULL,

    initialize = function(gdd_column = GDD_COLUMN, id = "gdd_threshold") {
      self$gdd_column <- gdd_column
      super$initialize(
        id = id,
        feature_types = c("numeric", "integer"),
        predict_types = "response",
        packages = c("stats", "data.table")
      )
    }
  ),
  private = list(
    .train = function(task) {
      target_column <- task$target_names
      data <- task$data(cols = c(target_column, self$gdd_column))
      setnames(data, c(target_column, self$gdd_column), c("stage", "gdd"))

      thresholds <- data[, .(gdd = median(gdd)), by = stage][order(stage)]
      if (nrow(thresholds) < 2L) {
        stop("GDD calibration requires at least two training stages.")
      }

      thresholds[, gdd := stats::isoreg(stage, gdd)$yf]
      thresholds <- thresholds[, .(stage = max(stage)), by = gdd][order(gdd)]
      if (nrow(thresholds) < 2L) {
        stop("Monotonic GDD calibration produced fewer than two thresholds.")
      }

      list(
        gdd = thresholds$gdd,
        stage = thresholds$stage,
        thresholds = thresholds
      )
    },

    .predict = function(task) {
      new_gdd <- task$data(cols = self$gdd_column)[[self$gdd_column]]
      response <- approx(
        x = self$model$gdd,
        y = self$model$stage,
        xout = new_gdd,
        rule = 2,
        ties = "ordered"
      )$y
      list(response = as.numeric(response))
    }
  )
)

# ------------------------------------------------------------------------------
# 9. SCHEME-MATCHED INNER RESAMPLING
# ------------------------------------------------------------------------------

build_site_grouped_inner <- function(inner_task, seed) {
  ids <- inner_task$row_ids
  labels <- make_balanced_group_folds(
    cohort$.site_id[ids],
    INNER_GROUP_FOLDS,
    seed
  )
  list(
    resampling = resampling_from_labels(inner_task, labels),
    strategy = "site_grouped_fallback"
  )
}

build_inner_resampling <- function(spec, inner_task, seed) {
  ids <- inner_task$row_ids
  scheme <- as.character(spec$scheme)

  result <- tryCatch({
    if (identical(scheme, "spatial_block")) {
      labels <- cohort$.spatial_fold[ids]
      if (uniqueN(labels) < 2L) stop("Too few remaining spatial blocks.")
      list(
        resampling = resampling_from_labels(inner_task, labels),
        strategy = "remaining_outer_spatial_blocks"
      )
    } else if (identical(scheme, "leave_one_year_out")) {
      labels <- make_contiguous_year_blocks(
        cohort$.season_year[ids],
        INNER_TEMPORAL_BLOCKS
      )
      list(
        resampling = resampling_from_labels(inner_task, labels),
        strategy = "contiguous_temporal_blocks"
      )
    } else if (identical(scheme, "environmental_zone")) {
      zone <- cohort$.env_zone[ids]
      labels <- make_balanced_group_folds(zone, INNER_ZONE_FOLDS, seed)
      list(
        resampling = resampling_from_labels(inner_task, labels),
        strategy = "environmental_zone_grouped"
      )
    } else if (identical(scheme, "expanding_origin")) {
      years <- sort(unique(cohort$.season_year[ids]))
      candidate_years <- years[
        seq_along(years) > INNER_FORWARD_MIN_TRAINING_YEARS
      ]
      if (length(candidate_years) < 2L) {
        stop("Too few years for forward inner tuning.")
      }

      n_splits <- min(INNER_FORWARD_SPLITS, length(candidate_years))
      positions <- unique(as.integer(round(seq(
        1,
        length(candidate_years),
        length.out = n_splits
      ))))
      validation_years <- candidate_years[positions]

      train_sets <- lapply(validation_years, function(year) {
        ids[cohort$.season_year[ids] < year]
      })
      test_sets <- lapply(validation_years, function(year) {
        ids[cohort$.season_year[ids] == year]
      })

      valid <- vapply(train_sets, length, integer(1)) > 0L &
        vapply(test_sets, length, integer(1)) > 0L
      train_sets <- train_sets[valid]
      test_sets <- test_sets[valid]

      list(
        resampling = resampling_from_sets(inner_task, train_sets, test_sets),
        strategy = "forward_chaining_years"
      )
    } else if (identical(scheme, "joint_space_time")) {
      row_space <- cohort$.spatial_fold[ids]
      row_time <- cohort$.joint_temporal_block[ids]
      combinations <- CJ(
        heldout_space = sort(unique(row_space)),
        heldout_time = sort(unique(row_time)),
        unique = TRUE
      )
      combinations[, n_test := vapply(seq_len(.N), function(i) {
        sum(
          row_space == combinations$heldout_space[i] &
            row_time == combinations$heldout_time[i]
        )
      }, integer(1))]
      combinations[, n_train := vapply(seq_len(.N), function(i) {
        sum(
          row_space != combinations$heldout_space[i] &
            row_time != combinations$heldout_time[i]
        )
      }, integer(1))]
      combinations <- combinations[n_test > 0L & n_train >= MIN_OUTER_TRAIN_ROWS]
      if (nrow(combinations) < 2L) {
        stop("Too few valid joint inner split combinations.")
      }

      set.seed(seed)
      combinations[, tie_break := runif(.N)]
      setorder(combinations, -n_test, tie_break)
      combinations <- combinations[seq_len(min(INNER_JOINT_SPLITS, .N))]

      train_sets <- lapply(seq_len(nrow(combinations)), function(i) {
        ids[
          row_space != combinations$heldout_space[i] &
            row_time != combinations$heldout_time[i]
        ]
      })
      test_sets <- lapply(seq_len(nrow(combinations)), function(i) {
        ids[
          row_space == combinations$heldout_space[i] &
            row_time == combinations$heldout_time[i]
        ]
      })

      list(
        resampling = resampling_from_sets(inner_task, train_sets, test_sets),
        strategy = "joint_space_time_margin_exclusion"
      )
    } else {
      stop("Unknown inner-resampling scheme: ", scheme)
    }
  }, error = function(error) {
    fallback <- build_site_grouped_inner(inner_task, seed + 100000L)
    fallback$strategy <- paste0(
      fallback$strategy,
      "__reason_",
      safe_id(conditionMessage(error))
    )
    fallback
  })

  result$manifest <- summarize_inner_resampling(
    result$resampling,
    result$strategy,
    cohort
  )
  result
}

# ------------------------------------------------------------------------------
# 10. SELECTED TREE LEARNER AND FROZEN TUNING SPACES
# ------------------------------------------------------------------------------

make_auto_tuner <- function(
  learner,
  search_space,
  evaluations,
  inner_resampling,
  id
) {
  auto_tuner(
    learner = learner,
    resampling = inner_resampling,
    measure = msr("regr.rmse"),
    tuner = tnr("mbo"),
    search_space = search_space,
    term_evals = evaluations,
    store_tuning_instance = TRUE,
    store_benchmark_result = FALSE,
    store_models = FALSE,
    check_values = TRUE,
    id = id
  )
}

make_selected_learner <- function(inner_resampling) {
  p <- length(FEATURES)

  if (identical(SELECTED_LEARNER, "random_forest")) {
    ensure_recommended_learner_registered("regr.ranger")
    base_learner <- mlr3::lrn(
      "regr.ranger",
      num.trees = RF_NUM_TREES,
      importance = "none",
      num.threads = 1L,
      seed = SEED,
      predict_type = "response"
    )
    return(make_auto_tuner(
      learner = base_learner,
      search_space = ps(
        mtry = p_int(1L, p),
        min.node.size = p_int(1L, 100L, logscale = TRUE),
        sample.fraction = p_dbl(0.5, 1),
        max.depth = p_int(0L, 20L)
      ),
      evaluations = TUNING_EVALS_RF,
      inner_resampling = inner_resampling,
      id = "random_forest"
    ))
  }

  ensure_recommended_learner_registered("regr.xgboost")
  base_learner <- mlr3::lrn(
    "regr.xgboost",
    objective = "reg:squarederror",
    tree_method = "hist",
    booster = "gbtree",
    nthread = 1L,
    seed = SEED,
    predict_type = "response"
  )
  make_auto_tuner(
    learner = base_learner,
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
    inner_resampling = inner_resampling,
    id = "xgboost"
  )
}

# ------------------------------------------------------------------------------
# 11. CHECKPOINTED OUTER JOBS
# ------------------------------------------------------------------------------

model_ids <- c(
  SELECTED_LEARNER,
  if (RUN_GDD_BASELINE) "gdd_threshold"
)

jobs <- outer_split_specs[
  rep(seq_len(nrow(outer_split_specs)), each = length(model_ids))
]
jobs[, model_id := rep(model_ids, times = nrow(outer_split_specs))]
jobs[, model_order := match(model_id, model_ids)]
setorder(jobs, scheme, split_id, model_order)
jobs[, model_order := NULL]
jobs[, job_id := .I]
setcolorder(jobs, c("job_id", "model_id", setdiff(names(jobs), c("job_id", "model_id"))))
write_table(jobs, "benchmark_jobs.csv")

input_md5 <- unname(tools::md5sum(INPUT_FILE))

run_outer_job <- function(job) {
  checkpoint_name <- paste0(
    sprintf("%04d", job$job_id), "__",
    safe_id(job$model_id), "__",
    safe_id(job$scheme), "__",
    safe_id(job$fold_id), ".rds"
  )
  checkpoint_file <- file.path(CHECKPOINT_DIR, checkpoint_name)

  if (RESUME_FROM_CHECKPOINTS && file.exists(checkpoint_file)) {
    cached <- tryCatch(readRDS(checkpoint_file), error = function(error) NULL)
    if (
      !is.null(cached) &&
        identical(cached$status, "ok") &&
        identical(cached$checkpoint_version, CHECKPOINT_VERSION) &&
        identical(cached$input_md5, input_md5) &&
        identical(as.character(cached$job$model_id), as.character(job$model_id)) &&
        identical(as.integer(cached$job$split_id), as.integer(job$split_id)) &&
        identical(as.character(cached$job$scheme), as.character(job$scheme)) &&
        identical(as.character(cached$job$fold_id), as.character(job$fold_id))
    ) {
      return(cached)
    }
  }

  set.seed(SEED + as.integer(job$job_id))
  ids <- resolve_outer_ids(job)
  train_ids <- ids$train_ids
  test_ids <- ids$test_ids

  full_task <- task$clone(deep = TRUE)
  inner_task <- task$clone(deep = TRUE)
  inner_task$filter(train_ids)

  warning_messages <- character()
  message_messages <- character()
  start_time <- proc.time()[["elapsed"]]

  result <- tryCatch(
    withCallingHandlers({
      inner_strategy <- NA_character_
      inner_manifest <- NULL

      if (identical(job$model_id, "gdd_threshold")) {
        learner <- LearnerRegrGDDThreshold$new(
          gdd_column = GDD_COLUMN,
          id = "gdd_threshold"
        )
      } else {
        inner <- build_inner_resampling(
          spec = job,
          inner_task = inner_task,
          seed = SEED + 1000L + as.integer(job$split_id)
        )
        learner <- make_selected_learner(inner$resampling)
        inner_strategy <- inner$strategy
        inner_manifest <- inner$manifest
      }

      learner$train(inner_task)
      prediction <- learner$predict(full_task, row_ids = test_ids)
      prediction_table <- prediction_to_table(prediction)

      prediction_table[, feature_set := FEATURE_SET_ID]
      prediction_table[, learner_id := as.character(job$model_id)]
      prediction_table[, scheme := as.character(job$scheme)]
      prediction_table[, split_id := as.integer(job$split_id)]
      prediction_table[, fold_id := as.character(job$fold_id)]

      prediction_error <- prediction_table$response - prediction_table$truth

      tuned_parameter_values <- NULL
      tuned_parameters <- NA_character_
      inner_rmse <- NA_real_
      tuning_archive <- NULL
      gdd_thresholds <- NULL

      if (inherits(learner, "AutoTuner")) {
        tuned_parameter_values <- learner$learner$param_set$values
        tuned_parameters <- params_to_text(tuned_parameter_values)
        tuning_result <- learner$tuning_result
        if (!is.null(tuning_result) && "regr.rmse" %in% names(tuning_result)) {
          inner_rmse <- tuning_result[["regr.rmse"]][1L]
        }
        if (SAVE_TUNING_ARCHIVES) {
          tuning_archive <- as.data.table(learner$archive$data)
        }
      } else if (identical(job$model_id, "gdd_threshold")) {
        gdd_thresholds <- copy(learner$model$thresholds)
        tuned_parameters <- "training_fold_median_stage_thresholds_isotonic"
      }

      list(
        status = "ok",
        job = copy(job),
        predictions = prediction_table,
        fold_metrics = data.table(
          feature_set = FEATURE_SET_ID,
          learner_id = as.character(job$model_id),
          scheme = as.character(job$scheme),
          split_id = as.integer(job$split_id),
          fold_id = as.character(job$fold_id),
          n_train = length(train_ids),
          n_test = nrow(prediction_table),
          n_excluded = length(ids$excluded_ids),
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
          inner_strategy = inner_strategy,
          inner_rmse = inner_rmse,
          tuned_parameters = tuned_parameters
        ),
        tuned_parameter_values = tuned_parameter_values,
        gdd_thresholds = gdd_thresholds,
        inner_manifest = inner_manifest,
        tuning_archive = tuning_archive,
        warnings = warning_messages,
        messages = message_messages,
        error = NA_character_
      )
    },
    warning = function(warning) {
      warning_messages <<- capture_condition(warning_messages, warning)
      invokeRestart("muffleWarning")
    },
    message = function(message_condition) {
      message_messages <<- capture_condition(message_messages, message_condition)
      invokeRestart("muffleMessage")
    }),
    error = function(error) {
      list(
        status = "error",
        job = copy(job),
        predictions = NULL,
        fold_metrics = data.table(
          feature_set = FEATURE_SET_ID,
          learner_id = as.character(job$model_id),
          scheme = as.character(job$scheme),
          split_id = as.integer(job$split_id),
          fold_id = as.character(job$fold_id),
          n_train = length(train_ids),
          n_test = length(test_ids),
          n_excluded = length(ids$excluded_ids),
          rmse = NA_real_,
          mae = NA_real_,
          bias = NA_real_,
          within_5_bbch = NA_real_,
          n_outside_0_99 = NA_integer_,
          n_extreme = NA_integer_,
          min_prediction = NA_real_,
          max_prediction = NA_real_,
          inner_strategy = NA_character_,
          inner_rmse = NA_real_,
          tuned_parameters = NA_character_
        ),
        tuned_parameter_values = NULL,
        gdd_thresholds = NULL,
        inner_manifest = NULL,
        tuning_archive = NULL,
        warnings = warning_messages,
        messages = message_messages,
        error = conditionMessage(error)
      )
    }
  )

  result$checkpoint_version <- CHECKPOINT_VERSION
  result$input_md5 <- input_md5
  result$elapsed_seconds <- proc.time()[["elapsed"]] - start_time
  atomic_save_rds(result, checkpoint_file, compress = FALSE)
  result
}

message(
  "Starting ", nrow(jobs), " outer jobs with ", workers, " worker(s)."
)
message("Selected learner: ", SELECTED_LEARNER)
message("Schemes: ", paste(CV_SCHEMES, collapse = ", "))
message("Checkpoints: ", CHECKPOINT_DIR)

options(future.globals.maxSize = 8 * 1024^3)
future::plan(future::multisession, workers = workers)
outer_results <- future.apply::future_lapply(
  seq_len(nrow(jobs)),
  function(i) run_outer_job(jobs[i]),
  future.packages = c(
    "data.table", "R6", "mlr3", "mlr3learners", "mlr3tuning",
    "mlr3mbo", "paradox", "mlr3spatiotempcv", "ranger", "xgboost"
  ),
  future.seed = TRUE,
  future.scheduling = 1
)
future::plan(future::sequential)

atomic_save_rds(
  outer_results,
  file.path(OUTPUT_DIR, "outer_results_compact.rds"),
  compress = FALSE
)

# ------------------------------------------------------------------------------
# 12. COMBINE RESULTS AND RETAIN REFIT SPECIFICATIONS FOR SHAP
# ------------------------------------------------------------------------------

job_status <- rbindlist(lapply(outer_results, function(result) {
  data.table(
    job_id = result$job$job_id,
    model_id = result$job$model_id,
    scheme = result$job$scheme,
    split_id = result$job$split_id,
    fold_id = result$job$fold_id,
    status = result$status,
    elapsed_seconds = result$elapsed_seconds,
    error = result$error,
    warnings = paste(result$warnings, collapse = " | "),
    messages = paste(result$messages, collapse = " | ")
  )
}), fill = TRUE)
setorder(job_status, job_id)
write_table(job_status, "job_status.csv")
write_table(job_status[status != "ok"], "failed_jobs.csv")

fold_metrics <- rbindlist(
  lapply(outer_results, `[[`, "fold_metrics"),
  fill = TRUE
)
write_table(fold_metrics, "outer_fold_metrics.csv")

successful_predictions <- lapply(outer_results, `[[`, "predictions")
successful_predictions <- successful_predictions[
  !vapply(successful_predictions, is.null, logical(1))
]
if (!length(successful_predictions)) stop("No outer job produced predictions.")

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
  env_zone = .env_zone,
  Date = .date,
  spatial_block = .spatial_fold,
  joint_temporal_block = .joint_temporal_block
)]

oof_predictions <- merge(
  oof_predictions,
  prediction_metadata,
  by = "row_ids",
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(oof_predictions$env_id)) {
  stop("Could not match all predictions to cohort metadata.")
}

oof_predictions[, error := response - truth]
oof_predictions[, absolute_error := abs(error)]
oof_predictions[, squared_error := error^2]
oof_predictions[, outside_0_99 := response < 0 | response > 99]
oof_predictions[, extreme_prediction :=
  response < EXTREME_PREDICTION_LOWER |
    response > EXTREME_PREDICTION_UPPER]

oof_predictions[, bootstrap_cluster := fcase(
  scheme == "spatial_block", site_id,
  scheme == "leave_one_year_out", as.character(season_year),
  scheme == "expanding_origin", as.character(season_year),
  scheme == "environmental_zone", env_zone,
  scheme == "joint_space_time", fold_id,
  default = env_id
)]

setcolorder(
  oof_predictions,
  c(
    "feature_set", "learner_id", "scheme", "split_id", "fold_id",
    "row_ids", "source_row_id", "env_id", "site_id", "season_year",
    "env_zone", "Date", "truth", "response", "error", "absolute_error",
    "squared_error", "outside_0_99", "extreme_prediction",
    "bootstrap_cluster", "spatial_block", "joint_temporal_block"
  )
)
write_table(oof_predictions, "oof_predictions.csv")

# Exact fold refit information for later OOF SHAP. Training rules are in
# outer_split_specs.csv and predictor membership is in benchmark_configuration.rds.
fold_model_specifications <- lapply(outer_results, function(result) {
  list(
    status = result$status,
    job = result$job,
    tuned_parameter_values = result$tuned_parameter_values,
    gdd_thresholds = result$gdd_thresholds,
    inner_manifest = result$inner_manifest,
    checkpoint_version = result$checkpoint_version,
    input_md5 = result$input_md5,
    error = result$error
  )
})
atomic_save_rds(
  fold_model_specifications,
  file.path(OUTPUT_DIR, "fold_model_specifications.rds"),
  compress = FALSE
)

inner_resampling_manifests <- rbindlist(lapply(outer_results, function(result) {
  if (is.null(result$inner_manifest)) return(NULL)
  answer <- copy(result$inner_manifest)
  answer[, job_id := result$job$job_id]
  answer[, learner_id := result$job$model_id]
  answer[, scheme := result$job$scheme]
  answer[, split_id := result$job$split_id]
  answer[, fold_id := result$job$fold_id]
  answer
}), fill = TRUE)
write_table(inner_resampling_manifests, "inner_resampling_manifests.csv")

if (SAVE_TUNING_ARCHIVES) {
  tuning_archives <- rbindlist(lapply(outer_results, function(result) {
    if (is.null(result$tuning_archive)) return(NULL)
    archive <- copy(result$tuning_archive)
    archive[, job_id := result$job$job_id]
    archive[, learner_id := result$job$model_id]
    archive[, scheme := result$job$scheme]
    archive[, split_id := result$job$split_id]
    archive[, fold_id := result$job$fold_id]
    archive
  }), fill = TRUE)

  atomic_save_rds(
    tuning_archives,
    file.path(OUTPUT_DIR, "tuning_archives.rds"),
    compress = FALSE
  )

  archive_list_columns <- names(tuning_archives)[vapply(
    tuning_archives, is.list, logical(1)
  )]
  tuning_archives_csv <- tuning_archives[
    , setdiff(names(tuning_archives), archive_list_columns), with = FALSE
  ]
  write_table(tuning_archives_csv, "tuning_archives_scalar.csv")
}

# ------------------------------------------------------------------------------
# 13. PERFORMANCE, STABILITY AND TRAJECTORY SUMMARIES
# ------------------------------------------------------------------------------

job_completion <- job_status[, .(
  expected_outer_jobs = .N,
  completed_outer_jobs = sum(status == "ok"),
  failed_outer_jobs = sum(status != "ok"),
  complete = all(status == "ok")
), by = .(learner_id = model_id, scheme)]

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
  feature_set, learner_id, scheme, BBCH_stage = truth
)]

# Primary units: one row per environment x observed BBCH stage. MSE is retained
# so macro-RMSE can be calculated as sqrt(mean(unit MSE)).
environment_stage_metrics <- oof_predictions[, .(
  n = .N,
  mae = mean(absolute_error),
  mse = mean(squared_error),
  rmse = sqrt(mean(squared_error)),
  bias = mean(error),
  bootstrap_cluster = bootstrap_cluster[1L]
), by = .(
  feature_set, learner_id, scheme, env_id, site_id, season_year,
  env_zone, BBCH_stage = truth
)]

macro_environment_stage <- environment_stage_metrics[, .(
  n_environment_stage_units = .N,
  macro_mae = mean(mae),
  macro_rmse = sqrt(mean(mse)),
  mean_environment_stage_rmse = mean(rmse),
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
  list(overall_metrics, macro_environment_stage, stage_extremes)
)
performance_summary <- merge(
  performance_summary,
  job_completion,
  by = c("learner_id", "scheme"),
  all.x = TRUE
)
setorder(performance_summary, scheme, learner_id)

write_table(overall_metrics, "overall_row_weighted_metrics.csv")
write_table(stage_metrics, "stage_specific_metrics.csv")
write_table(environment_stage_metrics, "environment_stage_metrics.csv")
write_table(performance_summary, "generalisation_performance_by_scheme.csv")

extreme_folds <- fold_metrics[
  n_extreme > 0L | n_outside_0_99 > 0L
][order(-n_extreme, -n_outside_0_99)]
write_table(extreme_folds, "extreme_prediction_folds.csv")

trajectory_diagnostics <- oof_predictions[
  order(learner_id, scheme, env_id, Date),
  .(
    n_observations = .N,
    n_prediction_decreases = if (.N > 1L) {
      sum(diff(response) < -1e-8)
    } else 0L,
    n_decreases_gt_1 = if (.N > 1L) {
      sum(diff(response) < -1)
    } else 0L,
    largest_prediction_decrease = if (.N > 1L) {
      min(c(0, diff(response)))
    } else 0
  ),
  by = .(feature_set, learner_id, scheme, env_id)
]

trajectory_summary <- trajectory_diagnostics[, .(
  n_environments = .N,
  environments_with_any_decrease = sum(n_prediction_decreases > 0L),
  environments_with_decrease_gt_1 = sum(n_decreases_gt_1 > 0L),
  proportion_environments_with_decrease_gt_1 = mean(n_decreases_gt_1 > 0L),
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
# 14. SCHEME-APPROPRIATE CLUSTER BOOTSTRAP
# ------------------------------------------------------------------------------

bootstrap_model_metrics <- function(unit_table, replicates, seed) {
  cluster_summary <- unit_table[, .(
    sum_mae = sum(mae),
    sum_mse = sum(mse),
    n_units = .N
  ), by = bootstrap_cluster]

  observed_mae <- mean(unit_table$mae)
  observed_rmse <- sqrt(mean(unit_table$mse))

  if (nrow(cluster_summary) < 2L) {
    return(data.table(
      macro_mae = observed_mae,
      macro_mae_ci_low = NA_real_,
      macro_mae_ci_high = NA_real_,
      macro_rmse = observed_rmse,
      macro_rmse_ci_low = NA_real_,
      macro_rmse_ci_high = NA_real_,
      n_bootstrap_clusters = nrow(cluster_summary),
      n_environment_stage_units = nrow(unit_table)
    ))
  }

  set.seed(seed)
  bootstrap_values <- replicate(replicates, {
    sampled <- sample.int(
      nrow(cluster_summary),
      nrow(cluster_summary),
      replace = TRUE
    )
    denominator <- sum(cluster_summary$n_units[sampled])
    c(
      mae = sum(cluster_summary$sum_mae[sampled]) / denominator,
      rmse = sqrt(sum(cluster_summary$sum_mse[sampled]) / denominator)
    )
  })

  data.table(
    macro_mae = observed_mae,
    macro_mae_ci_low = as.numeric(quantile(bootstrap_values["mae", ], 0.025)),
    macro_mae_ci_high = as.numeric(quantile(bootstrap_values["mae", ], 0.975)),
    macro_rmse = observed_rmse,
    macro_rmse_ci_low = as.numeric(quantile(bootstrap_values["rmse", ], 0.025)),
    macro_rmse_ci_high = as.numeric(quantile(bootstrap_values["rmse", ], 0.975)),
    n_bootstrap_clusters = nrow(cluster_summary),
    n_environment_stage_units = nrow(unit_table)
  )
}

bootstrap_groups <- unique(environment_stage_metrics[, .(learner_id, scheme)])
performance_bootstrap <- rbindlist(lapply(seq_len(nrow(bootstrap_groups)), function(i) {
  group <- bootstrap_groups[i]
  units <- environment_stage_metrics[
    learner_id == group$learner_id & scheme == group$scheme
  ]
  result <- bootstrap_model_metrics(
    units,
    BOOTSTRAP_REPLICATES,
    SEED + 10000L + i
  )
  result[, learner_id := group$learner_id]
  result[, scheme := group$scheme]
  result
}))
setcolorder(
  performance_bootstrap,
  c("learner_id", "scheme", setdiff(names(performance_bootstrap), c("learner_id", "scheme")))
)
setorder(performance_bootstrap, scheme, learner_id)
write_table(
  performance_bootstrap,
  "cluster_bootstrap_performance_by_scheme.csv"
)

bootstrap_paired_difference <- function(paired_units, replicates, seed) {
  cluster_summary <- paired_units[, .(
    sum_delta_mae = sum(mae_selected - mae_baseline),
    sum_mse_selected = sum(mse_selected),
    sum_mse_baseline = sum(mse_baseline),
    n_units = .N
  ), by = bootstrap_cluster]

  observed_delta_mae <- mean(
    paired_units$mae_selected - paired_units$mae_baseline
  )
  observed_delta_rmse <-
    sqrt(mean(paired_units$mse_selected)) -
    sqrt(mean(paired_units$mse_baseline))

  if (nrow(cluster_summary) < 2L) {
    return(data.table(
      delta_macro_mae_selected_minus_gdd = observed_delta_mae,
      delta_macro_mae_ci_low = NA_real_,
      delta_macro_mae_ci_high = NA_real_,
      probability_selected_better_mae = NA_real_,
      delta_macro_rmse_selected_minus_gdd = observed_delta_rmse,
      delta_macro_rmse_ci_low = NA_real_,
      delta_macro_rmse_ci_high = NA_real_,
      probability_selected_better_rmse = NA_real_,
      n_bootstrap_clusters = nrow(cluster_summary),
      n_paired_units = nrow(paired_units)
    ))
  }

  set.seed(seed)
  bootstrap_values <- replicate(replicates, {
    sampled <- sample.int(
      nrow(cluster_summary),
      nrow(cluster_summary),
      replace = TRUE
    )
    denominator <- sum(cluster_summary$n_units[sampled])
    delta_mae <- sum(cluster_summary$sum_delta_mae[sampled]) / denominator
    selected_rmse <- sqrt(
      sum(cluster_summary$sum_mse_selected[sampled]) / denominator
    )
    baseline_rmse <- sqrt(
      sum(cluster_summary$sum_mse_baseline[sampled]) / denominator
    )
    c(mae = delta_mae, rmse = selected_rmse - baseline_rmse)
  })

  data.table(
    delta_macro_mae_selected_minus_gdd = observed_delta_mae,
    delta_macro_mae_ci_low = as.numeric(quantile(bootstrap_values["mae", ], 0.025)),
    delta_macro_mae_ci_high = as.numeric(quantile(bootstrap_values["mae", ], 0.975)),
    probability_selected_better_mae = mean(bootstrap_values["mae", ] < 0),
    delta_macro_rmse_selected_minus_gdd = observed_delta_rmse,
    delta_macro_rmse_ci_low = as.numeric(quantile(bootstrap_values["rmse", ], 0.025)),
    delta_macro_rmse_ci_high = as.numeric(quantile(bootstrap_values["rmse", ], 0.975)),
    probability_selected_better_rmse = mean(bootstrap_values["rmse", ] < 0),
    n_bootstrap_clusters = nrow(cluster_summary),
    n_paired_units = nrow(paired_units)
  )
}

if (RUN_GDD_BASELINE) {
  selected_units <- environment_stage_metrics[
    learner_id == SELECTED_LEARNER,
    .(
      scheme, env_id, BBCH_stage, bootstrap_cluster,
      mae_selected = mae,
      mse_selected = mse
    )
  ]
  baseline_units <- environment_stage_metrics[
    learner_id == "gdd_threshold",
    .(
      scheme, env_id, BBCH_stage,
      bootstrap_cluster_baseline = bootstrap_cluster,
      mae_baseline = mae,
      mse_baseline = mse
    )
  ]
  paired_units <- merge(
    selected_units,
    baseline_units,
    by = c("scheme", "env_id", "BBCH_stage"),
    all = FALSE
  )
  if (any(
    paired_units$bootstrap_cluster != paired_units$bootstrap_cluster_baseline
  )) {
    stop("Selected learner and GDD baseline have mismatched bootstrap clusters.")
  }

  paired_schemes <- sort(unique(paired_units$scheme))
  selected_vs_gdd <- rbindlist(lapply(seq_along(paired_schemes), function(i) {
    scheme_id <- paired_schemes[i]
    result <- bootstrap_paired_difference(
      paired_units[scheme == scheme_id],
      BOOTSTRAP_REPLICATES,
      SEED + 20000L + i
    )
    result[, selected_learner := SELECTED_LEARNER]
    result[, scheme := scheme_id]
    result
  }))
  setcolorder(
    selected_vs_gdd,
    c("selected_learner", "scheme", setdiff(names(selected_vs_gdd), c("selected_learner", "scheme")))
  )
  write_table(
    selected_vs_gdd,
    "selected_learner_vs_gdd_paired_bootstrap.csv"
  )
}

# ------------------------------------------------------------------------------
# 15. CONFIGURATION AND PROVENANCE
# ------------------------------------------------------------------------------

configuration <- list(
  run_id = RUN_ID,
  scientific_stage = "frozen_pipeline_multischeme_generalisation",
  input_file = normalizePath(INPUT_FILE),
  input_md5 = input_md5,
  seed = SEED,
  target = TARGET,
  include_bbch_zero = TRUE,
  feature_set = FEATURE_SET_ID,
  features = FEATURES,
  selected_learner = SELECTED_LEARNER,
  comparator_models = if (RUN_GDD_BASELINE) "gdd_threshold" else character(),
  site_group_column = site_column,
  environment_column = environment_column,
  year_column = year_column,
  date_column = date_column,
  zone_column = zone_column,
  cv_schemes = CV_SCHEMES,
  scheme_claims = list(
    spatial_block = "new physical locations/spatial blocks",
    leave_one_year_out = "unseen climatic season with past and future years permitted in training",
    expanding_origin = "future-season forecasting from earlier seasons only",
    environmental_zone = "new climatic-edaphic environmental zone",
    joint_space_time = "new spatial block in a held-out temporal block"
  ),
  spatial_folds = SPATIAL_FOLDS,
  spatial_block_size_m = SPATIAL_BLOCK_SIZE_M,
  spatial_selection = SPATIAL_SELECTION,
  joint_temporal_blocks = JOINT_TEMPORAL_BLOCKS,
  joint_outer_training_rule = "spatial_fold != heldout_space AND temporal_block != heldout_time",
  joint_outer_test_rule = "spatial_fold == heldout_space AND temporal_block == heldout_time",
  expanding_min_training_years = EXPANDING_MIN_TRAINING_YEARS,
  inner_group_folds = INNER_GROUP_FOLDS,
  inner_temporal_blocks = INNER_TEMPORAL_BLOCKS,
  inner_zone_folds = INNER_ZONE_FOLDS,
  inner_forward_splits = INNER_FORWARD_SPLITS,
  inner_forward_min_training_years = INNER_FORWARD_MIN_TRAINING_YEARS,
  inner_joint_splits = INNER_JOINT_SPLITS,
  tuning_measure = "regr.rmse",
  tuner = "mbo",
  tuning_evals = list(
    random_forest = TUNING_EVALS_RF,
    xgboost = TUNING_EVALS_XGB
  ),
  package_versions = list(
    mlr3 = as.character(utils::packageVersion("mlr3")),
    mlr3learners = as.character(utils::packageVersion("mlr3learners")),
    mlr3tuning = as.character(utils::packageVersion("mlr3tuning")),
    mlr3mbo = as.character(utils::packageVersion("mlr3mbo")),
    ranger = as.character(utils::packageVersion("ranger")),
    xgboost = as.character(utils::packageVersion("xgboost"))
  ),
  rf_num_trees = RF_NUM_TREES,
  primary_metric = "environment_stage_macro_rmse",
  macro_rmse_definition = "sqrt(mean(environment_stage_mse))",
  bootstrap_replicates = BOOTSTRAP_REPLICATES,
  bootstrap_clusters = list(
    spatial_block = "physical_site",
    leave_one_year_out = "season_year",
    expanding_origin = "season_year",
    environmental_zone = "environmental_zone",
    joint_space_time = "joint_outer_space_time_cell"
  ),
  predictions_clipped = FALSE,
  extreme_prediction_limits = c(
    EXTREME_PREDICTION_LOWER,
    EXTREME_PREDICTION_UPPER
  ),
  checkpoint_version = CHECKPOINT_VERSION,
  workers = workers,
  n_input_rows = nrow(raw_data),
  n_cohort_rows = nrow(cohort),
  n_environments = uniqueN(cohort$.env_id),
  n_sites = uniqueN(cohort$.site_id),
  n_years = uniqueN(cohort$.season_year),
  n_environmental_zones = uniqueN(cohort$.env_zone, na.rm = TRUE),
  n_outer_splits = nrow(outer_split_specs),
  n_jobs = nrow(jobs)
)

atomic_save_rds(
  configuration,
  file.path(OUTPUT_DIR, "benchmark_configuration.rds")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "session_info.txt")
)

message("Multi-scheme generalisation benchmark finished.")
message("Results written to: ", OUTPUT_DIR)
message(
  "Primary performance table: ",
  file.path(OUTPUT_DIR, "generalisation_performance_by_scheme.csv")
)

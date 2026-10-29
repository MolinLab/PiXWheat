# ==============================================================================
# Independent BRIWECS validation of XGBoost phenology predictions
# ==============================================================================
# Reconstruct the German-Austrian development cohort, use its matching frozen
# model cache (or explicitly refit on development data), and predict BRIWECS
# environment-level milestone targets and daily threshold-crossing dates.
# A targeted two-predictor GAM comparison tests the thermal-time by radiation
# interaction. BRIWECS outcomes are not used to fit or tune either model.
# Run from the repository root or set PHENOLOGY_PROJECT_DIR to its path.

# ------------------------------------------------------------------------------
# 1. USER CONFIGURATION
# ------------------------------------------------------------------------------
PROJECT_DIR <- Sys.getenv("PHENOLOGY_PROJECT_DIR", unset = getwd())
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = TRUE)

BENCHMARK_RUN_ID <- Sys.getenv(
  "PHENOLOGY_GENERALISATION_RUN_ID",
  unset = "multischeme_xgboost_paper_hybrid"
)
BENCHMARK_DIR <- Sys.getenv(
  "PHENOLOGY_GENERALISATION_BENCHMARK_DIR",
  unset = file.path(PROJECT_DIR, "benchmark_results", BENCHMARK_RUN_ID)
)

BRIWECS_DIR <- Sys.getenv(
  "PHENOLOGY_BRIWECS_DIR",
  unset = file.path(PROJECT_DIR, "Data", "BRIWECS")
)
BRIWECS_ENVIRONMENT_TARGETS_FILE <- Sys.getenv(
  "PHENOLOGY_BRIWECS_ENVIRONMENT_TARGETS_FILE",
  unset = file.path(
    Sys.getenv(
      "PHENOLOGY_BRIWECS_TARGET_OUTPUT_DIR",
      unset = file.path(PROJECT_DIR, "external_validation", "briwecs_environment_targets")
    ),
    "briwecs_environment_milestone_targets_paper_style.csv"
  )
)
REFERENCE_TREATMENT <- trimws(Sys.getenv(
  "PHENOLOGY_BRIWECS_REFERENCE_TREATMENT", unset = "HN_WF_RF"
))
TARGET_DEFINITION <- "paper_style_environment_v1"
BRIWECS_LOCATION_FILE <- Sys.getenv(
  "PHENOLOGY_BRIWECS_LOCATION_FILE",
  unset = file.path(BRIWECS_DIR, "BRIWECS_location_coordinates.csv")
)
BRIWECS_WEATHER_DIR <- Sys.getenv(
  "PHENOLOGY_BRIWECS_WEATHER_DIR",
  unset = BRIWECS_DIR
)
BRIWECS_WEATHER_PATTERN <- Sys.getenv(
  "PHENOLOGY_BRIWECS_WEATHER_PATTERN",
  unset = "^DWD_HYRAS_2014_2022_all_locations.*[.]csv$"
)
# Accept the current ERA5-Land extraction directory and its legacy spelling.
era5_dir_override <- trimws(Sys.getenv(
  "PHENOLOGY_BRIWECS_ERA5_DIR",
  unset = ""
))
era5_dir_candidates <- unique(c(
  era5_dir_override,
  file.path(
    PROJECT_DIR, "Data", "Era5",
    "era5land_yearly_BRIWECS", "temp_extract"
  ),
  file.path(
    PROJECT_DIR, "Data", "Era5",
    "era5lan_yearly_BRIWECS", "temp_extract"
  )
))
era5_dir_candidates <- era5_dir_candidates[nzchar(era5_dir_candidates)]
existing_era5_dirs <- era5_dir_candidates[dir.exists(era5_dir_candidates)]
BRIWECS_ERA5_DIR <- if (length(existing_era5_dirs)) {
  existing_era5_dirs[1L]
} else {
  era5_dir_candidates[1L]
}

# The preferred file has one row per BRIWECS location, a PAW value (or the
# field-capacity layers needed to derive it), and a wilting-point storage value
# (or theta_wp_l1--l3). These fields are required to reproduce FTSW.
soil_override <- trimws(Sys.getenv("PHENOLOGY_BRIWECS_SOIL_FILE", unset = ""))
soil_candidates <- unique(c(
  soil_override,
  file.path(BRIWECS_DIR, "BRIWECS_soil_hydraulic.csv"),
  file.path(PROJECT_DIR, "Data", "Soil_hydraulic_ERA5_0_100.csv")
))
soil_candidates <- soil_candidates[nzchar(soil_candidates)]
existing_soil_candidates <- soil_candidates[file.exists(soil_candidates)]
BRIWECS_SOIL_FILE <- if (length(existing_soil_candidates)) {
  existing_soil_candidates[1L]
} else {
  soil_candidates[1L]
}

# Select exactly one deployment representation per run. The daily BRIWECS table
# still calculates the union of both feature sets, which makes switching a
# configuration choice rather than a feature-engineering code change. Examples:
#   PHENOLOGY_BRIWECS_FEATURE_SET=cumulative Rscript <this_script.R>
#   PHENOLOGY_BRIWECS_FEATURE_SET=hybrid     Rscript <this_script.R>
FEATURE_SET_REQUEST <- tolower(trimws(Sys.getenv(
  "PHENOLOGY_BRIWECS_FEATURE_SET",
  unset = "hybrid"
)))
FEATURE_SET_ALIASES <- c(
  hybrid = "hybrid_process_recent_0_100",
  hybrid_process_recent_0_100 = "hybrid_process_recent_0_100",
  cumulative = "cumulative_weather_0_100",
  cumulative_weather_0_100 = "cumulative_weather_0_100"
)
if (!FEATURE_SET_REQUEST %in% names(FEATURE_SET_ALIASES)) {
  stop(
    "PHENOLOGY_BRIWECS_FEATURE_SET must be one of: hybrid, cumulative, ",
    "hybrid_process_recent_0_100, cumulative_weather_0_100."
  )
}
FEATURE_SET <- unname(FEATURE_SET_ALIASES[[FEATURE_SET_REQUEST]])
FEATURE_SET_SHORT <- if (
  identical(FEATURE_SET, "hybrid_process_recent_0_100")
) "hybrid" else "cumulative"

VALIDATION_RUN_ID <- Sys.getenv(
  "PHENOLOGY_BRIWECS_RUN_ID",
  unset = paste0("briwecs_xgboost_", FEATURE_SET_SHORT, "_environment_targets_v1")
)
OUTPUT_DIR <- file.path(
  PROJECT_DIR,
  "external_validation",
  VALIDATION_RUN_ID
)
MODEL_DIR <- file.path(OUTPUT_DIR, "model")
PLOT_DIR <- file.path(OUTPUT_DIR, "plots")
# Reuse development-only models from the previous run without mixing its
# plot-based outputs into the environment-target run. Override for a custom run.
MODEL_SOURCE_DIR <- Sys.getenv(
  "PHENOLOGY_BRIWECS_MODEL_SOURCE_DIR",
  unset = file.path(
    PROJECT_DIR, "external_validation",
    paste0("briwecs_xgboost_", FEATURE_SET_SHORT, "_external_v2.2.1"), "model"
  )
)

# Optional path to an already completed cumulative-feature BRIWECS run. It is
# used only when the selected representation is hybrid. When its environment-
# unit table exists, the script writes a paired hybrid-versus-cumulative
# comparison without changing or recalibrating either model.
CUMULATIVE_VALIDATION_DIR <- Sys.getenv(
  "PHENOLOGY_BRIWECS_CUMULATIVE_RESULTS_DIR",
  unset = file.path(
    PROJECT_DIR,
    "external_validation",
    "briwecs_xgboost_cumulative_environment_targets_v1"
  )
)

SEED <- 123L
FINAL_TUNING_SCHEME <- trimws(Sys.getenv(
  "PHENOLOGY_BRIWECS_TUNING_SCHEME",
  unset = "joint_space_time"
))
ALLOWED_TUNING_SCHEMES <- c("joint_space_time", "spatial_block")
if (!FINAL_TUNING_SCHEME %in% ALLOWED_TUNING_SCHEMES) {
  stop(
    "PHENOLOGY_BRIWECS_TUNING_SCHEME must be one of: ",
    paste(ALLOWED_TUNING_SCHEMES, collapse = ", ")
  )
}

read_env_integer <- function(name, default, minimum = 1L) {
  raw <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(raw)) return(as.integer(default))
  value <- suppressWarnings(as.numeric(raw))
  parsed <- as.integer(value)
  if (!is.finite(value) || value != parsed || parsed < minimum) {
    stop(name, " must be an integer >= ", minimum, ".")
  }
  parsed
}

read_env_boolean <- function(name, default = FALSE) {
  raw <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(raw)) return(isTRUE(default))
  normalized <- tolower(raw)
  if (normalized %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (normalized %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(name, " must be TRUE or FALSE.")
}

# Updating observed targets must not silently retune or refit the predictor.
# Enable only for an intentional first fit on DE+AT, never to optimise BRIWECS.
ALLOW_MODEL_REFIT <- read_env_boolean("PHENOLOGY_BRIWECS_ALLOW_MODEL_REFIT", FALSE)
# A new interaction pair needs new GAMs, but must not force XGBoost retuning.
# Retain the existing all-model opt-in as the default for backwards compatibility.
ALLOW_GAM_ABLATION_REFIT <- read_env_boolean(
  "PHENOLOGY_BRIWECS_ALLOW_GAM_REFIT", ALLOW_MODEL_REFIT
)

TUNING_EVALUATIONS <- read_env_integer(
  "PHENOLOGY_BRIWECS_TUNING_EVALS",
  60L
)
FINAL_TUNING_FOLDS <- read_env_integer(
  "PHENOLOGY_BRIWECS_TUNING_FOLDS",
  if (FINAL_TUNING_SCHEME == "joint_space_time") 5L else 5L,
  minimum = 2L
)
XGBOOST_THREADS <- read_env_integer(
  "PHENOLOGY_BRIWECS_XGBOOST_THREADS",
  3L
)
BOOTSTRAP_REPLICATES <- read_env_integer(
  "PHENOLOGY_BRIWECS_BOOTSTRAP_REPLICATES",
  2000L,
  minimum = 100L
)
MAX_DAYS_AFTER_SOWING <- read_env_integer(
  "PHENOLOGY_BRIWECS_MAX_DAYS_AFTER_SOWING",
  365L,
  minimum = 200L
)
SOIL_MAX_DISTANCE_M <- read_env_integer(
  "PHENOLOGY_BRIWECS_SOIL_MAX_DISTANCE_M",
  1000L
)
ERA5_MAX_DISTANCE_M <- read_env_integer(
  "PHENOLOGY_BRIWECS_ERA5_MAX_DISTANCE_M",
  20000L
)

RUN_GAM_INTERACTION_ABLATION_REQUESTED <- read_env_boolean(
  "PHENOLOGY_BRIWECS_RUN_GAM_INTERACTION_ABLATION",
  TRUE
)
RUN_GAM_INTERACTION_ABLATION <- RUN_GAM_INTERACTION_ABLATION_REQUESTED
GAM_MAIN_K <- read_env_integer(
  "PHENOLOGY_BRIWECS_GAM_MAIN_K",
  10L,
  minimum = 4L
)
GAM_INTERACTION_K <- read_env_integer(
  "PHENOLOGY_BRIWECS_GAM_INTERACTION_K",
  5L,
  minimum = 3L
)
GAM_THREADS <- read_env_integer(
  "PHENOLOGY_BRIWECS_GAM_THREADS",
  XGBOOST_THREADS
)

INCLUDE_GDD_BASELINE <- TRUE
INCLUDE_MONOTONE_SENSITIVITY <- TRUE
MODEL_CACHE_VERSION <- 4L
RECENT_LOOKBACK_DAYS <- 30L
XGBOOST_TUNING_MEASURE_ID <- "regr.rmse"

# Exact feature-engineering constants from calculate predictors scripts.
T_BASE <- 0
T_OPTIMUM <- 27.5
T_CEILING <- 40
COLD_THRESHOLD <- 0
HEAT_THRESHOLD <- 30

# Exact phase-aware Soltani/Sinclair constants used by calculate predictors scripts.
# German and Austrian feature scripts.
CPP <- 21
PPSEN <- 0.004
SUN_ANGLE_BELOW_HORIZON <- 6
TBV <- -1
TP1V <- 0
TP2V <- 8
TCV <- 12
VDSAT <- 50
VSEN <- 0.030
BD_SOW_EMR <- 4.5
BD_PHOTO_START <- 4.5
BD_PHOTO_END <- 17.8
BD_VERN_START <- 4.5
BD_VERN_END <- 17.8
DEVERN_CUMVER_LIMIT <- 10
DEVERN_TMAX_LIMIT <- 30
DEVERN_RATE <- 0.5

# Native ERA5-Land layer thicknesses used by the training feature scripts.
SWVL1_DEPTH_MM <- 70
SWVL2_DEPTH_MM <- 210
SWVL3_DEPTH_MM <- 720
SW_DEPTH_0_28_MM <- SWVL1_DEPTH_MM + SWVL2_DEPTH_MM
SW_DEPTH_0_100_MM <- SW_DEPTH_0_28_MM + SWVL3_DEPTH_MM

# ------------------------------------------------------------------------------
# 2. PACKAGES AND BASIC UTILITIES
# ------------------------------------------------------------------------------

required_packages <- c(
  "data.table",
  "R6",
  "mlr3",
  "mlr3learners",
  "mlr3spatiotempcv",
  "mlr3tuning",
  "mlr3mbo",
  "paradox",
  "bbotk",
  "xgboost",
  "mgcv",
  "ggplot2",
  "geosphere"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install these packages before running BRIWECS validation: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(R6)
  library(mlr3)
  library(mlr3learners)
  library(mlr3spatiotempcv)
  library(mlr3tuning)
  library(mlr3mbo)
  library(paradox)
  library(bbotk)
  library(xgboost)
  library(mgcv)
  library(ggplot2)
  library(geosphere)
})

RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)
data.table::setDTthreads(1L)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MODEL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
message(
  "Selected BRIWECS representation: ", FEATURE_SET,
  " (set PHENOLOGY_BRIWECS_FEATURE_SET=hybrid or cumulative to switch)."
)

assert_file <- function(path, label) {
  if (!file.exists(path)) stop(label, " does not exist: ", path)
  invisible(path)
}

resolve_column <- function(data, candidates, label, required = TRUE) {
  hits <- candidates[candidates %in% names(data)]
  if (!length(hits)) {
    if (!required) return(NA_character_)
    stop(
      "Could not resolve ", label, ". Tried: ",
      paste(candidates, collapse = ", ")
    )
  }
  hits[1L]
}

numeric_locale <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  suppressWarnings(as.numeric(gsub(",", ".", trimws(as.character(x)), fixed = TRUE)))
}

parse_target_logical <- function(x, label) {
  value <- tolower(trimws(as.character(x)))
  result <- rep(NA, length(value))
  result[value %in% c("true", "t", "1")] <- TRUE
  result[value %in% c("false", "f", "0")] <- FALSE
  invalid <- !is.na(value) & !value %in% c("", "na", "true", "t", "1", "false", "f", "0")
  if (any(invalid)) stop(label, " contains invalid logical values.")
  result
}

parse_target_date <- function(x, label) {
  if (inherits(x, "Date")) return(as.IDate(x))
  text <- trimws(as.character(x))
  text[text %in% c("", "NA")] <- NA_character_
  result <- as.IDate(text, format = "%Y-%m-%d")
  if (any(!is.na(text) & (is.na(result) |
      format(result, "%Y-%m-%d") != text))) {
    stop(label, " must contain ISO calendar dates (YYYY-MM-DD), not DOY values.")
  }
  result
}

prepare_environment_targets <- function(data, reference_treatment, max_days) {
  targets <- copy(as.data.table(data))
  requested_treatment <- trimws(as.character(reference_treatment))
  if (length(requested_treatment) != 1L || is.na(requested_treatment) || !nzchar(requested_treatment)) {
    stop("Select exactly one nonempty reference treatment.")
  }
  required <- c(
    "environment_id", "Location", "SeasonYear", "milestone", "BBCH",
    "reference_treatment", "reference_sowing_date", "target_available",
    "observed_environment_milestone_date", "observed_environment_days_after_sowing",
    "target_source", "n_plots", "n_genotypes", "target_standard_error_days",
    "target_lower_95_das", "target_upper_95_das"
  )
  missing <- setdiff(required, names(targets))
  if (length(missing)) {
    stop("Paper-style target CSV is missing: ", paste(missing, collapse = ", "),
         ". Regenerate it with calculate_briwecs_environment_targets_paper_style.R.")
  }
  if (!nrow(targets)) stop("The paper-style target CSV is empty.")
  for (column in c("environment_id", "Location", "milestone", "reference_treatment", "target_source")) {
    set(targets, j = column, value = trimws(as.character(targets[[column]])))
  }
  targets[, target_row_id := .I]
  targets[, target_available := parse_target_logical(target_available, "target_available")]
  if (anyNA(targets$target_available)) stop("target_available must be TRUE or FALSE for every input row.")
  input_audit <- targets[, .(n_input_targets = .N), by = .(reference_treatment, target_available)]
  input_audit[, selected_for_validation :=
    !is.na(reference_treatment) & reference_treatment == requested_treatment & target_available]
  # Compute the mask outside data.table's column scope.
  selected_rows <- which(
    targets$target_available & !is.na(targets$reference_treatment) &
      targets$reference_treatment == requested_treatment
  )
  targets <- targets[selected_rows]
  if (!nrow(targets)) {
    stop("No available environment targets for reference treatment ", reference_treatment, ".")
  }
  for (column in c("environment_id", "Location", "milestone", "target_source")) {
    if (anyNA(targets[[column]]) || any(!nzchar(targets[[column]]))) {
      stop("Available targets have missing/empty ", column, " values.")
    }
  }
  numeric_columns <- c(
    "SeasonYear", "BBCH", "observed_environment_days_after_sowing", "n_plots",
    "n_genotypes", "target_standard_error_days", "target_lower_95_das", "target_upper_95_das"
  )
  for (column in numeric_columns) set(targets, j = column, value = numeric_locale(targets[[column]]))
  for (column in c("SeasonYear", "BBCH", "n_plots", "n_genotypes")) {
    values <- targets[[column]]
    if (any(!is.finite(values) | values < 1 | values != floor(values))) {
      stop(column, " must contain positive integers for available targets.")
    }
    set(targets, j = column, value = as.integer(values))
  }
  if (targets[, any(n_genotypes > n_plots)]) stop("Target genotype counts exceed plot counts.")
  if (targets[, any(!milestone %in% c("BBCH10", "BBCH59", "BBCH87") |
      milestone != paste0("BBCH", BBCH))]) {
    stop("Target milestone and BBCH must agree and be BBCH10, BBCH59 or BBCH87.")
  }
  for (column in c("reference_sowing_date", "observed_environment_milestone_date")) {
    set(targets, j = column, value = parse_target_date(targets[[column]], column))
    if (anyNA(targets[[column]])) stop("Available targets contain missing ", column, ".")
  }
  if (!"adjustment_singular" %in% names(targets)) targets[, adjustment_singular := NA]
  targets[, adjustment_singular := parse_target_logical(adjustment_singular, "adjustment_singular")]
  invalid_das <- targets[!is.finite(observed_environment_days_after_sowing) |
    observed_environment_days_after_sowing < 0 | observed_environment_days_after_sowing > max_days]
  if (nrow(invalid_das)) stop("Target days after sowing are missing or outside the fixed prediction horizon.")
  if (targets[, any(as.numeric(observed_environment_milestone_date - reference_sowing_date) !=
      round(observed_environment_days_after_sowing))]) {
    stop("Target calendar dates do not equal reference_sowing_date + round(target DAS). Regenerate the target CSV.")
  }
  if (targets[, any(is.finite(target_standard_error_days) & target_standard_error_days < 0)] ||
      targets[, any(is.finite(target_lower_95_das) & is.finite(target_upper_95_das) &
        target_lower_95_das > target_upper_95_das)]) {
    stop("Target uncertainty fields contain a negative SE or reversed interval.")
  }
  if (targets[, anyDuplicated(.SD), .SDcols = c("Location", "SeasonYear", "milestone")]) {
    stop("Targets are not unique by Location x SeasonYear x milestone. Do not average duplicate targets.")
  }
  environment_map <- unique(targets[, .(Location, SeasonYear, environment_id)])
  if (anyDuplicated(environment_map$environment_id) ||
      environment_map[, anyDuplicated(.SD), .SDcols = c("Location", "SeasonYear")]) {
    stop("environment_id must map one-to-one to Location x SeasonYear.")
  }
  sowing_audit <- targets[, .(
    n_reference_sowing_dates = uniqueN(reference_sowing_date),
    reference_sowing_date = min(reference_sowing_date),
    n_milestone_targets = .N
  ), by = .(environment_id, Location, SeasonYear)]
  if (sowing_audit[, any(n_reference_sowing_dates != 1L)]) {
    stop("Multiple reference sowing dates within an environment. Regenerate the paper-style targets with one environment-level reference.")
  }
  setorder(targets, Location, SeasonYear, BBCH)
  chronology <- targets[, .(
    invalid_stage_order = any(diff(observed_environment_days_after_sowing) < 0)
  ), by = environment_id]
  if (chronology[, any(invalid_stage_order)]) stop("Environment target milestone chronology is reversed.")
  targets[, `:=`(
    environment_milestone_id = paste(environment_id, milestone, sep = "__"),
    sowing_date = reference_sowing_date,
    sowing_DOY = as.integer(format(reference_sowing_date, "%j")),
    observed_date = observed_environment_milestone_date,
    BBCH_obs = as.numeric(BBCH),
    target_date_rounding_days = as.numeric(observed_environment_milestone_date - reference_sowing_date) -
      observed_environment_days_after_sowing
  )]
  list(targets = targets, input_audit = input_audit, sowing_audit = sowing_audit)
}

make_environment_stage_units <- function(predictions) {
  units <- copy(as.data.table(predictions))
  keys <- c("learner_id", "environment_id", "milestone")
  if (units[, anyDuplicated(.SD), .SDcols = keys]) {
    stop("Expected exactly one prediction per learner x environment x milestone.")
  }
  if (any(!is.finite(units$prediction))) stop("Non-finite environment target prediction.")
  # Retain the established mean_* export names for downstream figure scripts.
  # Each is now a SINGLE target prediction/loss, never a mean of plot losses.
  units[, `:=`(
    n_target_rows = 1L, truth = BBCH_obs,
    mean_prediction = prediction,
    mean_error = prediction - BBCH_obs,
    mean_abs_error = abs(prediction - BBCH_obs),
    mean_squared_error = (prediction - BBCH_obs)^2,
    proportion_within_3 = as.numeric(abs(prediction - BBCH_obs) <= 3),
    proportion_within_5 = as.numeric(abs(prediction - BBCH_obs) <= 5),
    proportion_within_10 = as.numeric(abs(prediction - BBCH_obs) <= 10),
    proportion_outside_0_99 = as.numeric(prediction < 0 | prediction > 99)
  )]
  units
}

read_compatible_model_cache <- function(
  paths, signature, allow_refit, label,
  refit_environment_variable = "PHENOLOGY_BRIWECS_ALLOW_MODEL_REFIT"
) {
  for (path in unique(paths)) {
    if (!file.exists(path)) next
    cached <- tryCatch(readRDS(path), error = function(error) NULL)
    if (is.list(cached) && identical(cached$signature, signature)) {
      return(list(object = cached, path = normalizePath(path)))
    }
  }
  if (!allow_refit) {
    stop("No compatible frozen ", label, " cache found. Checked: ",
         paste(unique(paths), collapse = ", "),
         ". Set PHENOLOGY_BRIWECS_MODEL_SOURCE_DIR to the existing model folder. ",
         "Only for an intentional development-only refit set ",
         refit_environment_variable, "=TRUE.")
  }
  list(object = NULL, path = NA_character_)
}

make_environment_date_offsets <- function(targets, threshold_dates) {
  keys <- c("trajectory_id", "Location", "SeasonYear", "BBCH_obs")
  if (threshold_dates[, anyDuplicated(.SD), .SDcols = c(keys, "learner_id")]) {
    stop("Predicted crossing dates are duplicated for an environment, stage and model.")
  }
  if (targets[, anyDuplicated(.SD), .SDcols = keys]) {
    stop("Observed milestone targets are duplicated in the date-offset join.")
  }
  result <- merge(targets, threshold_dates, by = keys, all.x = TRUE, sort = FALSE)
  n_models <- uniqueN(threshold_dates$learner_id)
  if (!n_models || nrow(result) != nrow(targets) * n_models || anyNA(result$learner_id) ||
      result[, anyDuplicated(.SD), .SDcols = c("learner_id", "environment_id", "milestone")]) {
    stop("Crossing-date coverage does not contain one row per target and model.")
  }
  result[, `:=`(
    threshold_reached = !is.na(predicted_date),
    predicted_days_after_sowing = as.numeric(predicted_date - reference_sowing_date),
    offset_days_rounded_date = as.numeric(predicted_date - observed_date)
  )]
  # Preserve fractional BLUE/mean DAS. The supplied rounded date is used only
  # for the daily BBCH lookup and the explicitly named rounded-date sensitivity.
  result[, offset_days := predicted_days_after_sowing - observed_environment_days_after_sowing]
  result[, `:=`(abs_offset_days = abs(offset_days), squared_offset_days = offset_days^2)]
  result
}

cumulative_target_compatibility <- function(current, reference) {
  keys <- c("environment_id", "Location", "SeasonYear", "milestone")
  metadata <- c(
    "target_definition", "target_input_md5", "reference_treatment", "reference_sowing_date",
    "observed_date", "observed_environment_days_after_sowing", "target_source", "truth", "n_target_rows"
  )
  measures <- c("mean_abs_error", "mean_squared_error", "mean_error")
  required <- c(keys, metadata, measures)
  if (length(setdiff(required, names(reference)))) {
    return(paste("Cumulative output is not an environment-target run; missing:",
                 paste(setdiff(required, names(reference)), collapse = ", ")))
  }
  if (!nrow(current) || !nrow(reference)) return("One representation has no XGBoost target units.")
  if (reference[, anyDuplicated(.SD), .SDcols = keys] ||
      current[, anyDuplicated(.SD), .SDcols = keys]) return("Duplicated paired target keys.")
  left <- copy(current[, ..required])
  right <- copy(reference[, ..required])
  setorderv(left, keys)
  setorderv(right, keys)
  if (nrow(left) != nrow(right)) return("Representations have different environment-milestone coverage.")
  for (column in c(keys, metadata)) {
    a <- left[[column]]
    b <- right[[column]]
    if (anyNA(a) || anyNA(b)) return(paste("Missing target-comparison metadata:", column))
    equal <- if (column %in% c("observed_environment_days_after_sowing", "truth", "n_target_rows")) {
      isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-10))
    } else {
      identical(as.character(a), as.character(b))
    }
    if (!equal) return(paste("Representations differ in target keys/definition:", column))
  }
  if (any(left$n_target_rows != 1L) || any(right$n_target_rows != 1L)) {
    return("Paired outputs contain aggregated plot losses rather than one target loss.")
  }
  if (any(!is.finite(as.matrix(right[, ..measures])))) return("Cumulative target losses are non-finite.")
  NULL
}

atomic_save_rds <- function(object, path, compress = TRUE) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = compress)
  if (!file.rename(temporary, path)) {
    stop("Could not atomically save: ", path)
  }
  invisible(path)
}

write_table <- function(data, filename) {
  path <- file.path(OUTPUT_DIR, filename)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  fwrite(as.data.table(data), temporary)
  if (!file.rename(temporary, path)) stop("Could not save: ", path)
  invisible(path)
}

format_list_cell_for_csv <- function(value) {
  if (is.null(value) || !length(value)) return(NA_character_)
  paste(capture.output(dput(value)), collapse = " ")
}

flatten_list_columns_for_csv <- function(data) {
  result <- copy(as.data.table(data))
  list_columns <- names(result)[vapply(result, is.list, logical(1))]
  for (column in list_columns) {
    set(
      result,
      j = column,
      value = vapply(
        result[[column]],
        format_list_cell_for_csv,
        character(1)
      )
    )
  }
  result
}

write_model_table <- function(data, filename) {
  path <- file.path(MODEL_DIR, filename)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  fwrite(flatten_list_columns_for_csv(data), temporary)
  if (!file.rename(temporary, path)) stop("Could not save: ", path)
  invisible(path)
}

package_version_or_na <- function(package) {
  if (requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else {
    NA_character_
  }
}

ensure_xgboost_registered <- function() {
  suppressPackageStartupMessages(requireNamespace("mlr3learners", quietly = TRUE))
  if (!"regr.xgboost" %in% mlr3::mlr_learners$keys()) {
    stop(
      "regr.xgboost is not registered. Restore the benchmark renv with ",
      "mlr3learners and xgboost installed."
    )
  }
  invisible(TRUE)
}

quantile_or_na <- function(x, probability) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(quantile(x, probability, names = FALSE))
}

# Process states require complete daily histories. A missing required daily
# value invalidates the state from that date onward rather than being treated
# as zero.
cumsum_complete_to_date <- function(x) {
  running_sum <- cumsum(fcoalesce(x, 0))
  running_missing <- cumsum(is.na(x))
  fifelse(running_missing == 0L, running_sum, NA_real_)
}

calc_civil_photoperiod <- function(
    latitude,
    date,
    sun_angle_below_horizon = SUN_ANGLE_BELOW_HORIZON
) {
  doy <- as.integer(format(as.Date(date), "%j"))
  radians <- pi / 180
  solar_mean_anomaly <- 0.9856 * doy - 3.251
  solar_longitude <- solar_mean_anomaly +
    1.916 * sin(solar_mean_anomaly * radians) +
    0.020 * sin(2 * solar_mean_anomaly * radians) +
    282.565
  declination <- asin(0.39779 * sin(solar_longitude * radians))
  latitude_radians <- as.numeric(latitude) * radians
  solar_altitude <- -(sun_angle_below_horizon * radians)
  cos_hour_angle <- (
    sin(solar_altitude) -
      sin(latitude_radians) * sin(declination)
  ) / (cos(latitude_radians) * cos(declination))
  cos_hour_angle <- pmin(1, pmax(-1, cos_hour_angle))
  24 * acos(cos_hour_angle) / pi
}

accumulate_vernalization <- function(verday, tmax) {
  result <- rep(NA_real_, length(verday))
  state <- 0
  complete_history <- TRUE
  for (i in seq_along(verday)) {
    if (!complete_history || is.na(verday[i]) || is.na(tmax[i])) {
      complete_history <- FALSE
      next
    }
    state <- state + verday[i]
    if (state < DEVERN_CUMVER_LIMIT && tmax[i] > DEVERN_TMAX_LIMIT) {
      state <- state - DEVERN_RATE * (tmax[i] - DEVERN_TMAX_LIMIT)
    }
    state <- max(0, state)
    result[i] <- state
  }
  result
}

calc_soltani_phase_state <- function(tempfun, ppfun, verfun) {
  n <- length(tempfun)
  ppfun_phase <- rep(NA_real_, n)
  verfun_phase <- rep(NA_real_, n)
  devfun_phase <- rep(NA_real_, n)
  cumulative_biological_days <- 0
  complete_history <- TRUE

  for (i in seq_len(n)) {
    if (
      !complete_history ||
        is.na(tempfun[i]) || is.na(ppfun[i]) || is.na(verfun[i])
    ) {
      complete_history <- FALSE
      next
    }
    photo_active <- cumulative_biological_days >= BD_PHOTO_START &&
      cumulative_biological_days <= BD_PHOTO_END
    vern_active <- cumulative_biological_days >= BD_VERN_START &&
      cumulative_biological_days <= BD_VERN_END
    ppfun_phase[i] <- if (photo_active) ppfun[i] else 1
    verfun_phase[i] <- if (vern_active) verfun[i] else 1
    devfun_phase[i] <- tempfun[i] * ppfun_phase[i] * verfun_phase[i]
    cumulative_biological_days <-
      cumulative_biological_days + devfun_phase[i]
  }
  list(
    ppfun_phase = ppfun_phase,
    verfun_phase = verfun_phase,
    devfun_phase = devfun_phase
  )
}

rolling_sum_available <- function(x, width) {
  running_sum <- cumsum(fcoalesce(x, 0))
  running_sum - shift(running_sum, n = width, fill = 0)
}

rolling_n_available <- function(x, width) {
  running_n <- cumsum(!is.na(x))
  running_n - shift(running_n, n = width, fill = 0)
}

rolling_calendar_complete <- function(date, width) {
  date_integer <- as.integer(date)
  seq_along(date_integer) >= width &
    date_integer - shift(date_integer, n = width - 1L) == width - 1L
}

rolling_sum_complete <- function(x, date, width) {
  complete <- rolling_calendar_complete(date, width) &
    rolling_n_available(x, width) == width
  result <- rep(NA_real_, length(x))
  result[complete] <- rolling_sum_available(x, width)[complete]
  result
}

rolling_mean_complete <- function(x, date, width) {
  rolling_sum_complete(x, date, width) / width
}

first_existing <- function(paths) {
  paths <- paths[nzchar(paths)]
  hits <- paths[file.exists(paths)]
  if (length(hits)) hits[1L] else NA_character_
}

save_plot <- function(plot, filename, width, height) {
  ggsave(
    file.path(PLOT_DIR, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
  )
  ggsave(
    file.path(PLOT_DIR, paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 450,
    limitsize = FALSE
  )
}

theme_pheno <- function() {
  theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9.5, margin = margin(b = 7)),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey94", color = "grey75"),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
}

get_gam_ablation_spec <- function(feature_set) {
  pairs <- list(
    hybrid_process_recent_0_100 = c("cumdevfun_phase", "sumRad_30"),
    cumulative_weather_0_100 = c("cumDTU_beta", "cumRad")
  )
  if (length(feature_set) != 1L || is.na(feature_set) ||
      !feature_set %in% names(pairs)) {
    stop("No GAM ablation specification exists for the selected feature set.")
  }
  features <- pairs[[feature_set]]
  is_hybrid <- identical(feature_set, "hybrid_process_recent_0_100")
  list(
    feature_set = feature_set,
    features = features,
    driver = features[1L],
    radiation = features[2L],
    pair_id = paste(features, collapse = "__X__"),
    label = if (is_hybrid) "Development x recent radiation" else "Thermal time x cumulative radiation",
    cache_filename = paste0("final_gam_", feature_set, "_ablation.rds"),
    # Only hybrid may reuse the original hybrid-only cache, after signature checks.
    legacy_cache_filename = if (is_hybrid) "final_gam_development_radiation_ablation.rds" else character(),
    figure_stem = if (is_hybrid) {
      "figure_gam_development_radiation_interaction_ablation"
    } else {
      "figure_gam_thermal_time_cumulative_radiation_interaction_ablation"
    }
  )
}

annotate_gam_ablation <- function(
  data, specification, reference_treatment, external_target_definition,
  external_target_md5 = NA_character_
) {
  result <- copy(as.data.table(data))
  if (!"dataset" %in% names(result) || anyNA(result$dataset) ||
      any(!result$dataset %in% c("Internal development CV", "External BRIWECS"))) {
    stop("GAM ablation output must identify internal development or external BRIWECS rows.")
  }
  external <- result$dataset == "External BRIWECS"
  # Assign explicitly outside data.table's column scope to avoid name collisions.
  metadata <- list(
    feature_set = rep(specification$feature_set, nrow(result)),
    interaction_pair = rep(specification$pair_id, nrow(result)),
    interaction_driver = rep(specification$driver, nrow(result)),
    interaction_radiation = rep(specification$radiation, nrow(result)),
    evaluation_population = ifelse(external, "BRIWECS_reference_treatment", "DE_AT_development_cohort"),
    reference_treatment = ifelse(external, reference_treatment, NA_character_),
    target_definition = ifelse(external, external_target_definition, "development_observations"),
    target_input_md5 = ifelse(external, external_target_md5, NA_character_)
  )
  for (column in names(metadata)) set(result, j = column, value = metadata[[column]])
  result
}

gam_effective_k <- function(x, requested_k, minimum_k = 3L) {
  n_unique <- data.table::uniqueN(x)
  if (n_unique <= minimum_k) {
    stop(
      "A GAM ablation predictor has only ", n_unique,
      " unique values; at least ", minimum_k + 1L, " are required."
    )
  }
  as.integer(max(minimum_k, min(requested_k, n_unique - 1L)))
}

make_gam_ablation_formula <- function(
  data, target, include_interaction, features = GAM_ABLATION_FEATURES
) {
  if (length(features) != 2L || anyNA(features) || anyDuplicated(features) ||
      !all(features %in% c("cumdevfun_phase", "sumRad_30", "cumDTU_beta", "cumRad"))) {
    stop("The GAM ablation requires two distinct supported predictor names.")
  }
  missing <- setdiff(c(target, features), names(data))
  if (length(missing)) {
    stop("GAM ablation data are missing: ", paste(missing, collapse = ", "), ".")
  }
  driver <- features[1L]
  radiation <- features[2L]
  driver_k <- gam_effective_k(
    data[[driver]],
    GAM_MAIN_K,
    minimum_k = 4L
  )
  radiation_k <- gam_effective_k(
    data[[radiation]],
    GAM_MAIN_K,
    minimum_k = 4L
  )
  interaction_driver_k <- gam_effective_k(
    data[[driver]],
    GAM_INTERACTION_K,
    minimum_k = 3L
  )
  interaction_radiation_k <- gam_effective_k(
    data[[radiation]],
    GAM_INTERACTION_K,
    minimum_k = 3L
  )

  terms <- c(
    sprintf("s(%s, k = %d, bs = 'cr')", driver, driver_k),
    sprintf("s(%s, k = %d, bs = 'cr')", radiation, radiation_k)
  )
  if (isTRUE(include_interaction)) {
    terms <- c(
      terms,
      sprintf(
        paste0(
          "ti(%s, %s, k = c(%d, %d), ",
          "bs = c('cr', 'cr'))"
        ),
        driver,
        radiation,
        interaction_driver_k,
        interaction_radiation_k
      )
    )
  }

  stats::as.formula(paste0(
    "`", target, "` ~ ",
    paste(terms, collapse = " + ")
  ))
}

fit_gam_ablation_model <- function(
  data, target, include_interaction, features = GAM_ABLATION_FEATURES
) {
  model_columns <- c(target, features)
  model_data <- as.data.frame(data[, ..model_columns])
  formula <- make_gam_ablation_formula(
    model_data,
    target,
    include_interaction,
    features = features
  )
  mgcv::bam(
    formula = formula,
    data = model_data,
    family = stats::gaussian(),
    method = "fREML",
    discrete = TRUE,
    select = TRUE,
    nthreads = GAM_THREADS,
    na.action = stats::na.fail
  )
}

predict_gam_ablation <- function(model, newdata, label) {
  prediction <- as.numeric(stats::predict(
    model,
    newdata = as.data.frame(newdata),
    type = "response"
  ))
  if (length(prediction) != nrow(newdata) || any(!is.finite(prediction))) {
    stop(label, " produced missing or non-finite predictions.")
  }
  prediction
}

# ------------------------------------------------------------------------------
# 3. RECONSTRUCT THE FROZEN DEVELOPMENT COHORT
# ------------------------------------------------------------------------------

assert_file(BRIWECS_ENVIRONMENT_TARGETS_FILE, "Paper-style environment-target CSV")
# Do not leave obsolete plot/treatment summaries beside the new primary results.
# The default run ID is new; old output folders are preserved, never cleaned here.
legacy_outputs <- c(
  "briwecs_replicate_weighted_performance_supplement.csv",
  "briwecs_performance_by_treatment.csv", "briwecs_performance_by_genotype.csv"
)
if (any(file.exists(file.path(OUTPUT_DIR, legacy_outputs)))) {
  stop("OUTPUT_DIR contains plot-based results. Use a new PHENOLOGY_BRIWECS_RUN_ID; old outputs have not been removed.")
}

configuration_path <- file.path(BENCHMARK_DIR, "benchmark_configuration.rds")
fold_manifest_path <- file.path(BENCHMARK_DIR, "outer_fold_manifest.csv")
split_specs_path <- file.path(BENCHMARK_DIR, "outer_split_specs.csv")

assert_file(configuration_path, "Benchmark configuration")
assert_file(fold_manifest_path, "Outer-fold manifest")
assert_file(split_specs_path, "Outer-split specifications")

benchmark_configuration <- readRDS(configuration_path)
if (!identical(
  tolower(as.character(benchmark_configuration$selected_learner)),
  "xgboost"
)) {
  stop("The selected generalisation learner is not XGBoost.")
}

TARGET <- as.character(benchmark_configuration$target)
FEATURE_SETS <- list(
  hybrid_process_recent_0_100 = c(
    "cumdevfun_phase",
    "sowing_DOY",
    "paw_fc33_0_100_mm",
    "FTSW_0_100_proxy",
    "sumRR_30",
    "sumRad_30",
    "sumCDD_30",
    "sumHDD_30",
    "sumPEvap_30",
    "climaticWaterBalance_30",
    "meanSoilWater_0_100_30"
  ),
  cumulative_weather_0_100 = c(
    "cumDTU_beta",
    "sowing_DOY",
    "photoperiod_civil",
    "cumRR",
    "cumRad",
    "cumCDD",
    "cumHDD",
    "paw_fc33_0_100_mm"
  )
)
if (!FEATURE_SET %in% names(FEATURE_SETS)) {
  stop("No feature definition exists for selected set: ", FEATURE_SET)
}
FEATURES <- FEATURE_SETS[[FEATURE_SET]]
benchmark_features <- as.character(benchmark_configuration$features)
if (length(benchmark_features) != length(FEATURES) ||
    !setequal(benchmark_features, FEATURES)) {
  stop(
    "The selected BRIWECS feature set does not match the frozen benchmark. ",
    "Benchmark features: ", paste(benchmark_features, collapse = ", "),
    "; requested features: ", paste(FEATURES, collapse = ", "), "."
  )
}
ALL_EXTERNAL_FEATURES <- unique(unlist(FEATURE_SETS, use.names = FALSE))
GAM_ABLATION_SPEC <- get_gam_ablation_spec(FEATURE_SET)
GAM_ABLATION_FEATURES <- GAM_ABLATION_SPEC$features
if (
  RUN_GAM_INTERACTION_ABLATION &&
    !all(GAM_ABLATION_FEATURES %in% FEATURES)
) {
  stop("The GAM ablation features are not part of the selected frozen feature set: ", FEATURE_SET, ".")
}
if (RUN_GAM_INTERACTION_ABLATION) {
  message(
    "GAM interaction comparison enabled for ", FEATURE_SET_SHORT, ": ",
    paste(GAM_ABLATION_FEATURES, collapse = " x "),
    "; GAM refitting permitted: ", ALLOW_GAM_ABLATION_REFIT,
    "; XGBoost refitting permitted: ", ALLOW_MODEL_REFIT, "."
  )
}

supported_external_features <- ALL_EXTERNAL_FEATURES
unsupported_features <- setdiff(FEATURES, supported_external_features)
if (length(unsupported_features)) {
  stop(
    "The BRIWECS builder does not yet implement frozen predictor(s): ",
    paste(unsupported_features, collapse = ", ")
  )
}

# Archived benchmark configurations may contain an obsolete server path.
# The repository copy is acceptable only if its checksum matches below.
recorded_training_file <- as.character(benchmark_configuration$input_file)
if (length(recorded_training_file) != 1L || !nzchar(recorded_training_file)) {
  stop("The benchmark configuration does not specify one development input.")
}
repo_training_file <- file.path(PROJECT_DIR, "Data", "pheno_df_CV.csv")
default_training_file <- if (file.exists(recorded_training_file)) {
  recorded_training_file
} else {
  repo_training_file
}
training_input_file <- Sys.getenv(
  "PHENOLOGY_INPUT_FILE", unset = default_training_file
)
assert_file(training_input_file, "Frozen development input")
training_input_md5 <- unname(tools::md5sum(training_input_file))
if (
  !is.null(benchmark_configuration$input_md5) &&
    !identical(training_input_md5, benchmark_configuration$input_md5)
) {
  stop("The development input MD5 differs from the generalisation benchmark.")
}

message("Reconstructing frozen development cohort.")
training_raw <- fread(
  training_input_file,
  na.strings = c("", "NA", "NaN", "Inf", "-Inf"),
  showProgress = TRUE
)
training_raw[, source_row_id_internal := .I]
fold_manifest <- fread(fold_manifest_path)
outer_split_specs <- fread(split_specs_path)

required_manifest_columns <- c(
  "obs_id", "source_row_id", "env_id", "site_id", "season_year",
  "BBCH", "lon", "lat", "spatial_block", "joint_temporal_block"
)
missing_manifest_columns <- setdiff(required_manifest_columns, names(fold_manifest))
if (length(missing_manifest_columns)) {
  stop(
    "Outer-fold manifest is missing: ",
    paste(missing_manifest_columns, collapse = ", ")
  )
}
if (anyDuplicated(fold_manifest$source_row_id) || anyDuplicated(fold_manifest$obs_id)) {
  stop("The fold manifest contains duplicated source_row_id or obs_id values.")
}

training_index <- match(
  fold_manifest$source_row_id,
  training_raw$source_row_id_internal
)
if (anyNA(training_index)) stop("Could not reconstruct all benchmark rows.")
training <- copy(training_raw[training_index])
training[, `.obs_id` := as.integer(fold_manifest$obs_id)]
training[, `.env_id` := as.character(fold_manifest$env_id)]
training[, `.site_id` := as.character(fold_manifest$site_id)]
training[, `.season_year` := as.integer(fold_manifest$season_year)]
training[, `.lon` := as.numeric(fold_manifest$lon)]
training[, `.lat` := as.numeric(fold_manifest$lat)]
training[, `.spatial_block` := as.character(fold_manifest$spatial_block)]
training[, `.joint_temporal_block` := as.character(
  fold_manifest$joint_temporal_block
)]

missing_training_columns <- setdiff(c(TARGET, FEATURES), names(training))
if (length(missing_training_columns)) {
  stop(
    "Frozen development columns are absent: ",
    paste(missing_training_columns, collapse = ", ")
  )
}
if (!"cumDTU_beta" %in% names(training)) {
  stop(
    "The development table needs cumDTU_beta for the unchanged thermal-time ",
    "baseline, although it is not necessarily a selected XGBoost predictor."
  )
}
if (anyNA(training[, c(TARGET, FEATURES), with = FALSE])) {
  stop("The frozen development cohort contains missing model values.")
}
if (any(!vapply(
  training[, c(TARGET, FEATURES), with = FALSE],
  function(x) is.numeric(x) && all(is.finite(x)),
  logical(1)
))) {
  stop("All frozen development model columns must be finite numeric values.")
}

backend_columns <- unique(c(
  ".obs_id", ".env_id", ".site_id", TARGET, ".lon", ".lat", FEATURES
))
training_backend <- training[, ..backend_columns]
training_task <- TaskRegrST$new(
  id = paste0(FEATURE_SET, "_final_deployment"),
  backend = training_backend,
  target = TARGET,
  coordinate_names = c(".lon", ".lat")
)
training_task$select(FEATURES)
MODEL_FEATURE_ORDER <- as.character(training_task$feature_names)
if (!setequal(MODEL_FEATURE_ORDER, FEATURES)) {
  stop("The reconstructed task does not contain the frozen predictor set.")
}
if (!identical(
  as.integer(training_task$row_ids),
  as.integer(training$.obs_id)
)) {
  stop("The reconstructed task row IDs differ from the benchmark manifest.")
}

training_feature_support <- rbindlist(lapply(FEATURES, function(feature) {
  x <- training[[feature]]
  data.table(
    feature = feature,
    train_min = min(x),
    train_q01 = quantile_or_na(x, 0.01),
    train_median = median(x),
    train_q99 = quantile_or_na(x, 0.99),
    train_max = max(x)
  )
}))
write_table(training_feature_support, "training_feature_support.csv")

# ------------------------------------------------------------------------------
# 4. DEVELOPMENT-ONLY FINAL TUNING AND FULL-COHORT REFIT
# ------------------------------------------------------------------------------

resolve_internal_split <- function(spec) {
  scheme <- as.character(spec$scheme)
  if (scheme == "spatial_block") {
    test_mask <- training$.spatial_block == as.character(
      spec$holdout_spatial_fold
    )
    train_mask <- !test_mask
  } else if (scheme == "joint_space_time") {
    heldout_space <- as.character(spec$holdout_spatial_fold)
    heldout_time <- as.character(spec$holdout_temporal_block)
    test_mask <- (
      training$.spatial_block == heldout_space &
        training$.joint_temporal_block == heldout_time
    )
    train_mask <- (
      training$.spatial_block != heldout_space &
        training$.joint_temporal_block != heldout_time
    )
  } else {
    stop("Unsupported final-tuning scheme: ", scheme)
  }
  train_mask[is.na(train_mask)] <- FALSE
  test_mask[is.na(test_mask)] <- FALSE
  list(
    train_ids = training$.obs_id[train_mask],
    test_ids = training$.obs_id[test_mask]
  )
}

tuning_specs <- copy(outer_split_specs[scheme == FINAL_TUNING_SCHEME])
if (!nrow(tuning_specs)) {
  stop(
    "The benchmark has no outer splits for final tuning scheme: ",
    FINAL_TUNING_SCHEME
  )
}
tuning_specs[, `:=`(
  n_train_internal = vapply(seq_len(.N), function(i) {
    length(resolve_internal_split(tuning_specs[i])$train_ids)
  }, integer(1)),
  n_test_internal = vapply(seq_len(.N), function(i) {
    length(resolve_internal_split(tuning_specs[i])$test_ids)
  }, integer(1))
)]
tuning_specs <- tuning_specs[n_train_internal > 0L & n_test_internal > 0L]
if (nrow(tuning_specs) < 2L) stop("Fewer than two valid final-tuning splits.")

if (nrow(tuning_specs) > FINAL_TUNING_FOLDS) {
  set.seed(SEED + 910L)
  tuning_specs[, tie_break__ := runif(.N)]
  setorder(tuning_specs, -n_test_internal, tie_break__)
  tuning_specs <- tuning_specs[seq_len(FINAL_TUNING_FOLDS)]
  tuning_specs[, tie_break__ := NULL]
}

train_sets <- lapply(seq_len(nrow(tuning_specs)), function(i) {
  resolve_internal_split(tuning_specs[i])$train_ids
})
test_sets <- lapply(seq_len(nrow(tuning_specs)), function(i) {
  resolve_internal_split(tuning_specs[i])$test_ids
})

internal_leakage_audit <- rbindlist(lapply(seq_along(train_sets), function(i) {
  train_ids <- train_sets[[i]]
  test_ids <- test_sets[[i]]
  data.table(
    tuning_iteration = i,
    split_id = tuning_specs$split_id[i],
    fold_id = tuning_specs$fold_id[i],
    n_train = length(train_ids),
    n_test = length(test_ids),
    n_environment_overlap = length(intersect(
      unique(training$.env_id[train_ids]),
      unique(training$.env_id[test_ids])
    )),
    n_site_overlap = length(intersect(
      unique(training$.site_id[train_ids]),
      unique(training$.site_id[test_ids])
    )),
    n_year_overlap = length(intersect(
      unique(training$.season_year[train_ids]),
      unique(training$.season_year[test_ids])
    ))
  )
}))
if (any(internal_leakage_audit$n_environment_overlap > 0L)) {
  stop("Environment leakage was detected in final tuning.")
}
if (FINAL_TUNING_SCHEME == "spatial_block" &&
    any(internal_leakage_audit$n_site_overlap > 0L)) {
  stop("Site leakage was detected in spatial final tuning.")
}
if (FINAL_TUNING_SCHEME == "joint_space_time" &&
    any(
      internal_leakage_audit$n_site_overlap > 0L |
        internal_leakage_audit$n_year_overlap > 0L
    )) {
  stop("Site/year leakage was detected in joint space-time final tuning.")
}
write_model_table(internal_leakage_audit, "final_tuning_split_audit.csv")

final_resampling <- rsmp("custom")
final_resampling$instantiate(
  training_task,
  train_sets = train_sets,
  test_sets = test_sets
)

ensure_xgboost_registered()
reconstruction_versions <- c("mlr3", "mlr3learners", "xgboost")
current_versions <- vapply(
  reconstruction_versions,
  package_version_or_na,
  character(1)
)
benchmark_versions <- vapply(reconstruction_versions, function(package) {
  value <- benchmark_configuration$package_versions[[package]]
  if (is.null(value)) NA_character_ else as.character(value)
}, character(1))
version_audit <- data.table(
  package = reconstruction_versions,
  benchmark_version = unname(benchmark_versions),
  current_version = unname(current_versions),
  exact_match = is.na(benchmark_versions) | benchmark_versions == current_versions
)
write_model_table(version_audit, "final_model_package_version_audit.csv")
if (any(!version_audit$exact_match)) {
  stop(
    "Restore the benchmark renv before fitting the final model. Version ",
    "mismatch: ",
    paste(
      version_audit[!exact_match, paste0(
        package, " benchmark=", benchmark_version,
        " current=", current_version
      )],
      collapse = "; "
    )
  )
}

model_signature <- paste(
  MODEL_CACHE_VERSION,
  training_input_md5,
  FEATURE_SET,
  paste(FEATURES, collapse = ","),
  FINAL_TUNING_SCHEME,
  paste(tuning_specs$split_id, collapse = ","),
  TUNING_EVALUATIONS,
  XGBOOST_TUNING_MEASURE_ID,
  XGBOOST_THREADS,
  SEED,
  paste(current_versions, collapse = ","),
  sep = "__"
)
model_cache_path <- file.path(MODEL_DIR, "final_xgboost_model.rds")
model_cache_result <- read_compatible_model_cache(
  c(model_cache_path, file.path(MODEL_SOURCE_DIR, "final_xgboost_model.rds")),
  model_signature, ALLOW_MODEL_REFIT, "XGBoost"
)
cached_model <- model_cache_result$object
model_loaded_from <- model_cache_result$path

if (!is.null(cached_model) && identical(cached_model$signature, model_signature)) {
  message("Reusing the completed full-cohort XGBoost deployment model.")
  final_learner <- cached_model$learner
  tuned_parameter_values <- cached_model$tuned_parameter_values
  tuning_archive <- cached_model$tuning_archive
  if (!identical(model_loaded_from, normalizePath(model_cache_path, mustWork = FALSE))) {
    atomic_save_rds(cached_model, model_cache_path, compress = FALSE)
  }
} else {
  message(
    "Tuning the final XGBoost model on development data only (",
    FINAL_TUNING_SCHEME, "; ", TUNING_EVALUATIONS, " evaluations)."
  )
  base_learner <- lrn(
    "regr.xgboost",
    objective = "reg:squarederror",
    tree_method = "hist",
    booster = "gbtree",
    nthread = XGBOOST_THREADS,
    seed = SEED,
    predict_type = "response"
  )
  search_space <- ps(
    nrounds = p_int(200L, 2500L, logscale = TRUE),
    eta = p_dbl(0.01, 0.3, logscale = TRUE),
    gamma = p_dbl(0, 5),
    max_depth = p_int(1L, 10L),
    min_child_weight = p_dbl(0.1, 30, logscale = TRUE),
    subsample = p_dbl(0.5, 1),
    colsample_bytree = p_dbl(0.5, 1),
    lambda = p_dbl(1e-4, 20, logscale = TRUE),
    alpha = p_dbl(1e-4, 20, logscale = TRUE)
  )
  final_auto_tuner <- auto_tuner(
    learner = base_learner,
    resampling = final_resampling,
    measure = msr(XGBOOST_TUNING_MEASURE_ID),
    tuner = tnr("mbo"),
    search_space = search_space,
    term_evals = TUNING_EVALUATIONS,
    store_tuning_instance = TRUE,
    store_benchmark_result = FALSE,
    store_models = FALSE,
    check_values = TRUE,
    id = "xgboost_final_briwecs"
  )
  final_auto_tuner$train(training_task)
  tuned_parameter_values <- final_auto_tuner$learner$param_set$values
  tuning_archive <- tryCatch(
    as.data.table(final_auto_tuner$tuning_instance$archive$data),
    error = function(error) data.table()
  )

  # AutoTuner fits its parametrized inner learner on the complete training task
  # after optimization. Store that plain, fully trained XGBoost learner.
  final_learner <- final_auto_tuner$learner$clone(deep = TRUE)

  cached_model <- list(
    cache_version = MODEL_CACHE_VERSION,
    signature = model_signature,
    learner = final_learner,
    tuned_parameter_values = tuned_parameter_values,
    tuning_archive = tuning_archive,
    feature_set = FEATURE_SET,
    configured_features = FEATURES,
    model_feature_order = MODEL_FEATURE_ORDER,
    tuning_scheme = FINAL_TUNING_SCHEME,
    tuning_split_ids = tuning_specs$split_id,
    tuning_evaluations = TUNING_EVALUATIONS,
    tuning_measure_id = XGBOOST_TUNING_MEASURE_ID,
    development_input_md5 = training_input_md5,
    fitted_at = format(Sys.time(), tz = "UTC")
  )
  atomic_save_rds(cached_model, model_cache_path, compress = FALSE)
}

if (nrow(tuning_archive)) {
  atomic_save_rds(
    tuning_archive,
    file.path(MODEL_DIR, "final_xgboost_tuning_archive.rds"),
    compress = FALSE
  )
  write_model_table(tuning_archive, "final_xgboost_tuning_archive.csv")
}
tuned_parameter_table <- data.table(
  parameter = names(tuned_parameter_values),
  value = vapply(tuned_parameter_values, function(x) {
    paste(as.character(x), collapse = ",")
  }, character(1))
)
write_model_table(tuned_parameter_table, "final_xgboost_parameters.csv")

# Targeted interaction comparison (H9 in the figure script): the two GAMs use
# the selected representation's same two main effects and differ only by the
# tensor-product interaction. This is not a feature/full-XGBoost ablation.
# All internal performance values are derived from held-out development
# environments; BRIWECS is not used here.
gam_additive_full <- NULL
gam_interaction_full <- NULL
internal_gam_ablation_units <- data.table()
gam_loaded_from <- NA_character_

if (RUN_GAM_INTERACTION_ABLATION) {
  # The ordered predictor pair is part of this signature. Keep the established
  # signature format so compatible frozen hybrid GAMs remain reusable; cumulative
  # uses a different pair. The XGBoost
  # cache version/signature are intentionally unchanged by this extension.
  gam_cache_signature <- paste(
    MODEL_CACHE_VERSION,
    training_input_md5,
    TARGET,
    paste(GAM_ABLATION_FEATURES, collapse = ","),
    FINAL_TUNING_SCHEME,
    paste(tuning_specs$split_id, collapse = ","),
    GAM_MAIN_K,
    GAM_INTERACTION_K,
    GAM_THREADS,
    SEED,
    package_version_or_na("mgcv"),
    sep = "__"
  )
  gam_cache_path <- file.path(
    MODEL_DIR,
    GAM_ABLATION_SPEC$cache_filename
  )
  gam_cache_filenames <- c(
    GAM_ABLATION_SPEC$cache_filename, GAM_ABLATION_SPEC$legacy_cache_filename
  )
  gam_cache_result <- read_compatible_model_cache(
    c(file.path(MODEL_DIR, gam_cache_filenames),
      file.path(MODEL_SOURCE_DIR, gam_cache_filenames)),
    gam_cache_signature, ALLOW_GAM_ABLATION_REFIT,
    paste0("GAM interaction-ablation (", GAM_ABLATION_SPEC$pair_id, ")"),
    refit_environment_variable = "PHENOLOGY_BRIWECS_ALLOW_GAM_REFIT"
  )
  cached_gam_ablation <- gam_cache_result$object
  gam_loaded_from <- gam_cache_result$path

  if (!is.null(cached_gam_ablation) &&
      identical(cached_gam_ablation$signature, gam_cache_signature)) {
    message("Reusing the completed GAM interaction-ablation models.")
    gam_additive_full <- cached_gam_ablation$gam_additive_full
    gam_interaction_full <- cached_gam_ablation$gam_interaction_full
    internal_gam_ablation_units <- as.data.table(
      cached_gam_ablation$internal_units
    )
    if (!identical(gam_loaded_from, normalizePath(gam_cache_path, mustWork = FALSE))) {
      atomic_save_rds(cached_gam_ablation, gam_cache_path, compress = FALSE)
    }
  } else {
    message(
      "Fitting additive and interaction GAMs in development-only grouped ",
      "folds (", FINAL_TUNING_SCHEME, ")."
    )
    internal_gam_ablation_units <- rbindlist(lapply(
      seq_along(train_sets),
      function(i) {
        development_train <- training[.obs_id %in% train_sets[[i]]]
        development_test <- training[.obs_id %in% test_sets[[i]]]
        if (!nrow(development_train) || !nrow(development_test)) {
          stop("An interaction-ablation fold is empty: ", i)
        }

        additive_model <- fit_gam_ablation_model(
          development_train,
          TARGET,
          include_interaction = FALSE
        )
        interaction_model <- fit_gam_ablation_model(
          development_train,
          TARGET,
          include_interaction = TRUE
        )
        prediction_data <- development_test[, ..GAM_ABLATION_FEATURES]
        fold_predictions <- data.table(
          model = rep(
            c("gam_additive", "gam_interaction"),
            each = nrow(development_test)
          ),
          environment_id = rep(development_test$.env_id, 2L),
          site_id = rep(development_test$.site_id, 2L),
          season_year = rep(development_test$.season_year, 2L),
          truth = rep(development_test[[TARGET]], 2L),
          prediction = c(
            predict_gam_ablation(
              additive_model,
              prediction_data,
              "Additive GAM"
            ),
            predict_gam_ablation(
              interaction_model,
              prediction_data,
              "Interaction GAM"
            )
          )
        )
        fold_predictions[, `:=`(
          error = prediction - truth,
          abs_error = abs(prediction - truth),
          squared_error = (prediction - truth)^2,
          milestone = paste0("BBCH", format(truth, trim = TRUE))
        )]
        fold_predictions[, .(
          n_records = .N,
          truth = first(truth),
          mean_prediction = mean(prediction),
          mean_error = mean(error),
          mean_abs_error = mean(abs_error),
          mean_squared_error = mean(squared_error)
        ), by = .(
          model,
          environment_id,
          site_id,
          season_year,
          milestone
        )][, `:=`(
          dataset = "Internal development CV",
          scheme = FINAL_TUNING_SCHEME,
          split_id = as.character(tuning_specs$split_id[i]),
          fold_id = as.character(tuning_specs$fold_id[i])
        )]
      }
    ), use.names = TRUE, fill = TRUE)

    message("Fitting both GAM variants on the complete development cohort.")
    gam_additive_full <- fit_gam_ablation_model(
      training,
      TARGET,
      include_interaction = FALSE
    )
    gam_interaction_full <- fit_gam_ablation_model(
      training,
      TARGET,
      include_interaction = TRUE
    )
    cached_gam_ablation <- list(
      cache_version = MODEL_CACHE_VERSION,
      signature = gam_cache_signature,
      feature_set = FEATURE_SET,
      ablation_features = GAM_ABLATION_FEATURES,
      interaction_pair = GAM_ABLATION_SPEC$pair_id,
      gam_additive_full = gam_additive_full,
      gam_interaction_full = gam_interaction_full,
      internal_units = internal_gam_ablation_units,
      reference_model = "gam_additive",
      alternative_model = "gam_interaction",
      formulas = c(
        additive = paste(deparse(stats::formula(gam_additive_full)), collapse = " "),
        interaction = paste(
          deparse(stats::formula(gam_interaction_full)),
          collapse = " "
        )
      ),
      fitted_at = format(Sys.time(), tz = "UTC")
    )
    atomic_save_rds(cached_gam_ablation, gam_cache_path, compress = FALSE)
  }

  if (!nrow(internal_gam_ablation_units)) {
    stop("The internal GAM interaction ablation produced no held-out units.")
  }
  for (model in list(gam_additive_full, gam_interaction_full)) {
    if (!setequal(all.vars(stats::formula(model)), c(TARGET, GAM_ABLATION_FEATURES))) {
      stop("A GAM model formula does not match the selected interaction pair and target.")
    }
  }
  internal_gam_ablation_units <- annotate_gam_ablation(
    internal_gam_ablation_units, GAM_ABLATION_SPEC,
    REFERENCE_TREATMENT, TARGET_DEFINITION
  )
  write_table(
    internal_gam_ablation_units,
    "gam_interaction_ablation_internal_environment_stage_units.csv"
  )
  write_model_table(
    data.table(
      model = c("gam_additive", "gam_interaction"),
      role = c("reference", "alternative"),
      feature_set = FEATURE_SET,
      interaction_pair = GAM_ABLATION_SPEC$pair_id,
      formula = c(
        paste(deparse(stats::formula(gam_additive_full)), collapse = " "),
        paste(deparse(stats::formula(gam_interaction_full)), collapse = " ")
      ),
      fitting_data = "complete frozen development cohort",
      external_outcomes_used = FALSE
    ),
    "gam_interaction_ablation_formulas.csv"
  )

  concurvity_result <- tryCatch(
    mgcv::concurvity(gam_interaction_full, full = TRUE),
    error = function(error) error
  )
  if (inherits(concurvity_result, "error")) {
    write_model_table(
      data.table(
        feature_set = FEATURE_SET,
        interaction_pair = GAM_ABLATION_SPEC$pair_id,
        status = "failed",
        message = conditionMessage(concurvity_result)
      ),
      "gam_interaction_ablation_concurvity.csv"
    )
  } else {
    concurvity_table <- as.data.table(
      as.data.frame(concurvity_result),
      keep.rownames = "diagnostic"
    )
    concurvity_table[, `:=`(
      status = "ok", feature_set = FEATURE_SET,
      interaction_pair = GAM_ABLATION_SPEC$pair_id
    )]
    write_model_table(
      concurvity_table,
      "gam_interaction_ablation_concurvity.csv"
    )
  }
}

# Full-development thermal-time comparator. This is calibrated before BRIWECS
# is opened and is never recalibrated on the external outcomes.
fit_thermal_baseline <- function(data, target, thermal_column) {
  thresholds <- data[, .(
    thermal_time = median(get(thermal_column))
  ), by = .(stage = get(target))][order(stage)]
  if (nrow(thresholds) < 2L) stop("Thermal baseline needs >= 2 stages.")
  thresholds[, thermal_time := stats::isoreg(stage, thermal_time)$yf]
  thresholds <- thresholds[
    , .(stage = max(stage)),
    by = thermal_time
  ][order(thermal_time)]
  if (nrow(thresholds) < 2L) {
    stop("Thermal baseline collapsed to fewer than two thresholds.")
  }
  thresholds
}

predict_thermal_baseline <- function(thermal_time, thresholds) {
  as.numeric(approx(
    x = thresholds$thermal_time,
    y = thresholds$stage,
    xout = thermal_time,
    rule = 2,
    ties = "ordered"
  )$y)
}

thermal_baseline_thresholds <- fit_thermal_baseline(
  training,
  TARGET,
  "cumDTU_beta"
)
write_model_table(
  thermal_baseline_thresholds,
  "full_development_thermal_baseline_thresholds.csv"
)

# ------------------------------------------------------------------------------
# 5. LOAD AND STANDARDISE THE INDEPENDENT BRIWECS OUTCOMES
# ------------------------------------------------------------------------------

assert_file(BRIWECS_ENVIRONMENT_TARGETS_FILE, "BRIWECS environment-target file")
assert_file(BRIWECS_LOCATION_FILE, "BRIWECS location file")
assert_file(BRIWECS_SOIL_FILE, "BRIWECS soil file")
if (!dir.exists(BRIWECS_ERA5_DIR)) {
  stop("BRIWECS ERA5-Land yearly directory does not exist: ", BRIWECS_ERA5_DIR)
}

message("Loading paper-style environment targets after freezing the development models.")
prepared_targets <- prepare_environment_targets(
  fread(BRIWECS_ENVIRONMENT_TARGETS_FILE, na.strings = c("", "NA", "NaN")),
  reference_treatment = REFERENCE_TREATMENT,
  max_days = MAX_DAYS_AFTER_SOWING
)
environment_targets <- prepared_targets$targets
target_input_md5 <- unname(tools::md5sum(BRIWECS_ENVIRONMENT_TARGETS_FILE))
environment_targets[, `:=`(
  feature_set = FEATURE_SET,
  evaluation_population = "BRIWECS_reference_treatment",
  target_definition = TARGET_DEFINITION,
  target_input_md5 = ..target_input_md5
)]
message(
  "External metrics use ", nrow(environment_targets), " available paper-style ",
  "environment-milestone targets for ", REFERENCE_TREATMENT,
  " only (including Overall summaries). Internal DE+AT CV is a separate population."
)
write_table(prepared_targets$input_audit, "briwecs_environment_target_input_audit.csv")
write_table(prepared_targets$sowing_audit, "briwecs_sowing_date_variation.csv")
date_qc <- environment_targets[, .(
  environment_milestone_id, environment_id, Location, SeasonYear, milestone,
  reference_sowing_date, observed_date, observed_environment_days_after_sowing,
  target_date_rounding_days, target_source,
  invalid_before_sowing = observed_date < reference_sowing_date,
  outside_fixed_horizon = observed_date > reference_sowing_date + MAX_DAYS_AFTER_SOWING
)]
write_table(date_qc, "briwecs_milestone_date_qc.csv")

location_data <- fread(
  BRIWECS_LOCATION_FILE,
  na.strings = c("", "NA", "NaN"),
  encoding = "UTF-8"
)
location_name_col <- resolve_column(
  location_data,
  c("Location", "location", "Standort", "standort"),
  "location-table site name"
)
location_lat_col <- resolve_column(
  location_data,
  c("Lat", "lat", "latitude", "Latitude"),
  "location latitude"
)
location_lon_col <- resolve_column(
  location_data,
  c("Long","Lon", "lon", "longitude", "Longitude"),
  "location longitude"
)
locations <- location_data[, .(
  Location = trimws(as.character(get(location_name_col))),
  lat = numeric_locale(get(location_lat_col)),
  lon = numeric_locale(get(location_lon_col))
)]
locations <- unique(locations)
if (anyDuplicated(locations$Location)) {
  duplicate_locations <- locations[, .N, by = Location][N > 1L, Location]
  stop(
    "Location coordinate table has non-unique locations: ",
    paste(duplicate_locations, collapse = ", ")
  )
}
# Restrict environmental extraction to sites represented by the selected targets.
missing_locations <- setdiff(environment_targets$Location, locations$Location)
if (length(missing_locations)) {
  stop("Environment targets lack coordinates for: ", paste(missing_locations, collapse = ", "))
}
locations <- locations[Location %in% environment_targets$Location]
if (locations[, any(!is.finite(lat) | !is.finite(lon) | abs(lat) > 90 | abs(lon) > 180)]) {
  stop("BRIWECS target coordinates are missing or invalid.")
}
environment_targets <- merge(
  environment_targets, locations, by = "Location", all.x = TRUE, sort = FALSE
)

# ------------------------------------------------------------------------------
# 6. SOIL CONTRACT AND DAILY WEATHER FEATURE ENGINEERING
# ------------------------------------------------------------------------------

soil_raw <- fread(
  BRIWECS_SOIL_FILE,
  na.strings = c("", "NA", "NaN", "Inf", "-Inf"),
  encoding = "UTF-8"
)
wp_layer_columns <- c("theta_wp_l1", "theta_wp_l2", "theta_wp_l3")
fc_layer_columns <- c("theta_fc33_l1", "theta_fc33_l2", "theta_fc33_l3")
wp_storage_column <- resolve_column(
  soil_raw,
  c(
    "wilting_point_0_100_mm",
    "wp_0_100_mm",
    "total_wilting_point_0_100_mm"
  ),
  "0-100 cm wilting-point water storage",
  required = FALSE
)

if (is.na(wp_storage_column)) {
  missing_wp_layers <- setdiff(wp_layer_columns, names(soil_raw))
  if (length(missing_wp_layers)) {
    stop(
      "FTSW_0_100_proxy requires wilting_point_0_100_mm or these columns: ",
      paste(wp_layer_columns, collapse = ", ")
    )
  }
  soil_raw[, wilting_point_0_100_mm :=
    SWVL1_DEPTH_MM * numeric_locale(theta_wp_l1) +
      SWVL2_DEPTH_MM * numeric_locale(theta_wp_l2) +
      SWVL3_DEPTH_MM * numeric_locale(theta_wp_l3)]
} else {
  soil_raw[, wilting_point_0_100_mm := numeric_locale(
    get(wp_storage_column)
  )]
}

if (!"paw_fc33_0_100_mm" %in% names(soil_raw)) {
  missing_fc_layers <- setdiff(fc_layer_columns, names(soil_raw))
  if (length(missing_fc_layers)) {
    stop(
      "Soil table needs paw_fc33_0_100_mm or these field-capacity columns: ",
      paste(fc_layer_columns, collapse = ", ")
    )
  }
  soil_raw[, field_capacity_0_100_mm :=
    SWVL1_DEPTH_MM * numeric_locale(theta_fc33_l1) +
      SWVL2_DEPTH_MM * numeric_locale(theta_fc33_l2) +
      SWVL3_DEPTH_MM * numeric_locale(theta_fc33_l3)]
  soil_raw[, paw_fc33_0_100_mm :=
    field_capacity_0_100_mm - wilting_point_0_100_mm]
}
soil_raw[, paw_fc33_0_100_mm := numeric_locale(paw_fc33_0_100_mm)]

soil_site_col <- resolve_column(
  soil_raw,
  c("Location", "location", "Standort", "standort"),
  "soil location",
  required = FALSE
)
soil_lat_col <- resolve_column(
  soil_raw,
  c("lat", "Lat", "latitude", "Latitude"),
  "soil latitude",
  required = FALSE
)
soil_lon_col <- resolve_column(
  soil_raw,
  c("Long","lon", "Lon", "longitude", "Longitude"),
  "soil longitude",
  required = FALSE
)
library(dplyr)
soil_raw <-soil_raw %>% filter(standort %in% c("KIE","GGE", "HAN", "KAL", "QLB", "RHH"))
soil_join_audit <- NULL
if (!is.na(soil_site_col)) {
  soil_by_site <- soil_raw[, .(
    Location = trimws(as.character(get(soil_site_col))),
    paw_fc33_0_100_mm,
    wilting_point_0_100_mm
  )]
  soil_by_site <- unique(soil_by_site)
  if (anyDuplicated(soil_by_site$Location)) {
    stop("Soil file has duplicated location names with different records.")
  }
  location_soil <- merge(
    locations,
    soil_by_site,
    by = "Location",
    all.x = TRUE,
    sort = FALSE
  )
    soil_join_audit <- location_soil[, .(
      Location,
      soil_match_method = "site_name",
      soil_distance_m = 0,
      paw_fc33_0_100_mm,
      wilting_point_0_100_mm
    )]
} else if (!is.na(soil_lat_col) && !is.na(soil_lon_col)) {
  soil_points <- soil_raw[, .(
    soil_lat = numeric_locale(get(soil_lat_col)),
    soil_lon = numeric_locale(get(soil_lon_col)),
    paw_fc33_0_100_mm,
    wilting_point_0_100_mm
  )]
  if (soil_points[, anyNA(soil_lat) | anyNA(soil_lon)]) {
    stop("Soil coordinates contain missing values.")
  }
  soil_join_audit <- rbindlist(lapply(seq_len(nrow(locations)), function(i) {
    distances <- geosphere::distHaversine(
      cbind(soil_points$soil_lon, soil_points$soil_lat),
      c(locations$lon[i], locations$lat[i])
    )
    nearest <- which.min(distances)
    data.table(
      Location = locations$Location[i],
      soil_match_method = "nearest_coordinate",
      soil_distance_m = distances[nearest],
      paw_fc33_0_100_mm = soil_points$paw_fc33_0_100_mm[nearest],
      wilting_point_0_100_mm =
        soil_points$wilting_point_0_100_mm[nearest]
    )
  }))
} else {
  stop("Soil file needs either a site-name column or latitude/longitude.")
}

# If a general soil table contains site names that do not match BRIWECS but
# also contains coordinates, use an explicitly audited nearest-point fallback.
if (
  soil_join_audit[, anyNA(paw_fc33_0_100_mm) |
    anyNA(wilting_point_0_100_mm)] &&
    !is.na(soil_lat_col) && !is.na(soil_lon_col)
) {
  soil_points <- soil_raw[, .(
    soil_lat = numeric_locale(get(soil_lat_col)),
    soil_lon = numeric_locale(get(soil_lon_col)),
    paw_fc33_0_100_mm,
    wilting_point_0_100_mm
  )][
    is.finite(soil_lat) & is.finite(soil_lon) &
      is.finite(paw_fc33_0_100_mm) &
      is.finite(wilting_point_0_100_mm)
  ]
  if (!nrow(soil_points)) stop("No finite coordinate-based soil records remain.")
  missing_locations <- soil_join_audit[
    is.na(paw_fc33_0_100_mm) | is.na(wilting_point_0_100_mm),
    Location
  ]
  nearest_fallback <- rbindlist(lapply(missing_locations, function(site_name) {
    site <- locations[Location == site_name]
    distances <- geosphere::distHaversine(
      cbind(soil_points$soil_lon, soil_points$soil_lat),
      c(site$lon, site$lat)
    )
    nearest <- which.min(distances)
    data.table(
      Location = site_name,
      soil_match_method = "nearest_coordinate_fallback",
      soil_distance_m = distances[nearest],
      paw_fc33_0_100_mm = soil_points$paw_fc33_0_100_mm[nearest],
      wilting_point_0_100_mm =
        soil_points$wilting_point_0_100_mm[nearest]
    )
  }))
  soil_join_audit <- rbind(
    soil_join_audit[!Location %in% missing_locations],
    nearest_fallback,
    use.names = TRUE
  )
}
setorder(soil_join_audit, Location)
write_table(soil_join_audit, "briwecs_soil_join_audit.csv")

if (soil_join_audit[, anyNA(paw_fc33_0_100_mm) |
    anyNA(wilting_point_0_100_mm) |
    any(!is.finite(paw_fc33_0_100_mm)) |
    any(!is.finite(wilting_point_0_100_mm))]) {
  stop(
    "Missing PAW or wilting-point storage for BRIWECS location(s): ",
    paste(
      soil_join_audit[
        !is.finite(paw_fc33_0_100_mm) |
          !is.finite(wilting_point_0_100_mm),
        Location
      ],
      collapse = ", "
    )
  )
}
if (soil_join_audit[, any(paw_fc33_0_100_mm <= 0)]) {
  stop("BRIWECS plant-available water capacity must be positive.")
}
if (soil_join_audit[, any(
  wilting_point_0_100_mm < 0 | wilting_point_0_100_mm > SW_DEPTH_0_100_MM
)]) {
  stop("BRIWECS 0-100 cm wilting-point storage must lie within 0-1000 mm.")
}
if (soil_join_audit[, any(soil_distance_m > SOIL_MAX_DISTANCE_M, na.rm = TRUE)]) {
  stop(
    "Nearest soil match exceeds ", SOIL_MAX_DISTANCE_M,
    " m. Inspect briwecs_soil_join_audit.csv or increase the explicit limit."
  )
}

weather_files <- list.files(
  BRIWECS_WEATHER_DIR,
  pattern = BRIWECS_WEATHER_PATTERN,
  full.names = TRUE
)
if (!length(weather_files)) {
  stop(
    "No BRIWECS weather files matched ", BRIWECS_WEATHER_PATTERN,
    " in ", BRIWECS_WEATHER_DIR
  )
}
weather_files <- sort(weather_files)
message("Reading ", length(weather_files), " BRIWECS weather files.")
weather_raw <- rbindlist(lapply(weather_files, function(path) {
  data <- fread(
    path,
    sep = ";",
    na.strings = c("", "NA", "NaN"),
    encoding = "UTF-8"
  )
  data[, source_weather_file := basename(path)]
  data
}), use.names = TRUE, fill = TRUE)

weather_location_col <- resolve_column(
  weather_raw,
  c("Location", "location", "Standort", "standort"),
  "weather location"
)
weather_date_col <- resolve_column(
  weather_raw,
  c("Date", "date"),
  "weather date"
)
weather_tmean_col <- resolve_column(
  weather_raw,
  c("Temperature.Air.Mean_C", "Tmean", "tmean", "TMP"),
  "daily mean temperature",
  required = FALSE
)
weather_tmin_col <- resolve_column(
  weather_raw,
  c("Temperature.Air.Min_C", "TN", "tmin"),
  "daily minimum temperature"
)
weather_tmax_col <- resolve_column(
  weather_raw,
  c("Temperature.Air.Max_C", "TX", "tmax"),
  "daily maximum temperature"
)
weather_rr_col <- resolve_column(
  weather_raw,
  c("Precipitation_mm_d", "RR", "prec", "precipitation"),
  "daily precipitation"
)
weather_rad_col <- resolve_column(
  weather_raw,
  c("Global.Radiation_MJ_m2_d", "Rad", "rad", "radiation"),
  "daily global radiation"
)

weather <- weather_raw[, .(
  Location = trimws(as.character(get(weather_location_col))),
  Date = as.IDate(get(weather_date_col)),
  Tmean_reported = if (is.na(weather_tmean_col)) NA_real_ else
    numeric_locale(get(weather_tmean_col)),
  TN = numeric_locale(get(weather_tmin_col)),
  TX = numeric_locale(get(weather_tmax_col)),
  RR = numeric_locale(get(weather_rr_col)),
  Rad_DWD = numeric_locale(get(weather_rad_col)),
  source_weather_file
)]
# Match the development scripts exactly: daily mean temperature is the mean of
# daily minimum and maximum, even when the source supplies another daily mean.
weather[, Tmean := (TN + TX) / 2]
weather_temperature_audit <- if (weather[, any(!is.na(Tmean_reported))]) {
  weather[!is.na(Tmean_reported), .(
    n = .N,
    mean_reported_minus_minmax = mean(Tmean_reported - Tmean),
    rmse_reported_vs_minmax = sqrt(mean((Tmean_reported - Tmean)^2)),
    maximum_absolute_difference = max(abs(Tmean_reported - Tmean))
  )]
} else {
  data.table(
    n = 0L,
    mean_reported_minus_minmax = NA_real_,
    rmse_reported_vs_minmax = NA_real_,
    maximum_absolute_difference = NA_real_
  )
}
write_table(
  weather_temperature_audit,
  "briwecs_reported_vs_minmax_mean_temperature.csv"
)
if (weather[, any(TN > TX, na.rm = TRUE)]) {
  stop("BRIWECS weather contains daily minimum temperature above maximum.")
}
if (weather[, any(RR < 0, na.rm = TRUE)]) {
  stop("BRIWECS weather contains negative precipitation.")
}
if (weather[, any(Rad_DWD < 0, na.rm = TRUE)]) {
  stop("BRIWECS weather contains negative global radiation.")
}
if (weather[, any(Rad_DWD > 60, na.rm = TRUE)]) {
  stop(
    "The DWD BRIWECS daily radiation field exceeds 60 MJ m-2 d-1. Verify ",
    "its units. The sumRad_30 predictor is calculated from ",
    "ERA5-Land ssrd, matching the training feature pipeline."
  )
}

weather_duplicate_qc <- weather[, .N, by = .(Location, Date)][N > 1L]
if (nrow(weather_duplicate_qc)) {
  write_table(weather_duplicate_qc, "briwecs_duplicate_weather_dates.csv")
  stop("BRIWECS weather has duplicate Location-Date rows.")
}
if (weather[, anyNA(Location) | anyNA(Date)]) {
  stop("BRIWECS weather contains missing location/date values.")
}

# ERA5-Land inputs are stored as one file per variable and calendar year,
# matching the corrected development-feature pipeline (for example
# pev_2025.csv). Read and standardize every required variable independently,
# then retain exactly one source point per BRIWECS location.
REQUIRED_ERA5_VARIABLES <- c(
  "ssrd", "pev", "swvl1", "swvl2", "swvl3", "swvl4"
)

standardize_era5_variable <- function(variable, directory, locations) {
  pattern <- paste0("^", variable, "_[0-9]{4}[.]csv$")
  files <- sort(list.files(directory, pattern = pattern, full.names = TRUE))
  if (!length(files)) {
    stop(
      "No ERA5-Land files found for ", variable, " in ", directory,
      ". Expected names such as ", variable, "_2025.csv."
    )
  }

  raw <- rbindlist(lapply(files, function(path) {
    value <- fread(
      path,
      na.strings = c("", "NA", "NaN", "Inf", "-Inf"),
      showProgress = FALSE
    )
    value[, source_era5_file := basename(path)]
    value
  }), use.names = TRUE, fill = TRUE)

  site_col <- resolve_column(
    raw,
    c(
      "Location", "location", "Standort", "standort",
      "Stations_id", "station_id", "Station"
    ),
    paste(variable, "ERA5-Land site identifier"),
    required = FALSE
  )
  date_col <- resolve_column(
    raw,
    c("Date", "date"),
    paste(variable, "ERA5-Land date")
  )
  lat_col <- resolve_column(
    raw,
    c("lat_orig", "lat", "Lat", "latitude", "Latitude"),
    paste(variable, "ERA5-Land latitude"),
    required = FALSE
  )
  lon_col <- resolve_column(
    raw,
    c("lon_orig", "lon", "Lon", "longitude", "Longitude"),
    paste(variable, "ERA5-Land longitude"),
    required = FALSE
  )
  if (is.na(site_col) && (is.na(lat_col) || is.na(lon_col))) {
    stop(
      variable,
      " files need a site identifier or lat/lon columns for spatial matching."
    )
  }

  standardized <- raw[, .(
    source_site = if (is.na(site_col)) {
      NA_character_
    } else {
      trimws(as.character(get(site_col)))
    },
    source_lat = if (is.na(lat_col)) NA_real_ else
      numeric_locale(get(lat_col)),
    source_lon = if (is.na(lon_col)) NA_real_ else
      numeric_locale(get(lon_col)),
    Date = as.IDate(get(date_col)),
    source_era5_file
  )]

  if (variable == "pev") {
    pev_col <- resolve_column(
      raw,
      c("pev", "PEV"),
      "native ERA5-Land potential evaporation",
      required = FALSE
    )
    pev_mm_col <- resolve_column(
      raw,
      c("potential_evaporation_mm", "PEvap_mm", "pev_mm"),
      "positive potential evaporation in mm",
      required = FALSE
    )
    if (is.na(pev_col) && is.na(pev_mm_col)) {
      stop(
        "The pev files need native pev or positive potential_evaporation_mm."
      )
    }
    standardized[, value := if (!is.na(pev_mm_col)) {
      numeric_locale(raw[[pev_mm_col]])
    } else {
      -1000 * numeric_locale(raw[[pev_col]])
    }]
    output_name <- "potential_evaporation_mm"
    value_transform <- if (!is.na(pev_mm_col)) {
      "potential_evaporation_mm retained"
    } else {
      "-1000 * native pev"
    }
  } else if (variable == "ssrd") {
    value_col <- resolve_column(
      raw,
      c("ssrd", "SSRD"),
      "ERA5-Land surface solar radiation"
    )
    standardized[, value := numeric_locale(raw[[value_col]]) / 1e6]
    output_name <- "Rad"
    value_transform <- "ssrd / 1e6"
  } else {
    value_col <- resolve_column(
      raw,
      c(variable, toupper(variable)),
      paste("ERA5-Land", variable)
    )
    standardized[, value := numeric_locale(raw[[value_col]])]
    output_name <- variable
    value_transform <- "retained in native volumetric units"
  }

  standardized[, source_key := paste(
    fcoalesce(source_site, "<no-site>"),
    fifelse(
      is.finite(source_lat),
      sprintf("%.7f", source_lat),
      "<no-lat>"
    ),
    fifelse(
      is.finite(source_lon),
      sprintf("%.7f", source_lon),
      "<no-lon>"
    ),
    sep = "__"
  )]
  source_sites <- unique(standardized[, .(
    source_key, source_site, source_lat, source_lon
  )])
  if (!nrow(source_sites)) stop("No ", variable, " ERA5-Land sites remain.")

  site_mapping <- rbindlist(lapply(seq_len(nrow(locations)), function(i) {
    location_name <- locations$Location[i]
    exact <- source_sites[
      !is.na(source_site) & source_site == location_name
    ]
    candidates <- if (nrow(exact)) exact else source_sites

    if (nrow(candidates) > 1L) {
      finite_coordinates <- candidates[
        is.finite(source_lat) & is.finite(source_lon)
      ]
      if (!nrow(finite_coordinates)) {
        stop(
          "Multiple ", variable, " source records match ", location_name,
          " but their coordinates are unavailable."
        )
      }
      distances <- geosphere::distHaversine(
        cbind(
          finite_coordinates$source_lon,
          finite_coordinates$source_lat
        ),
        c(locations$lon[i], locations$lat[i])
      )
      selected <- finite_coordinates[which.min(distances)]
      distance <- min(distances)
    } else {
      selected <- candidates[1L]
      distance <- if (
        is.finite(selected$source_lat) &&
          is.finite(selected$source_lon)
      ) {
        geosphere::distHaversine(
          c(selected$source_lon, selected$source_lat),
          c(locations$lon[i], locations$lat[i])
        )
      } else {
        0
      }
    }

    match_method <- if (
      !is.na(selected$source_site) &&
        selected$source_site == location_name
    ) {
      "site_identifier"
    } else {
      "nearest_coordinate"
    }
    data.table(
      variable = variable,
      Location = location_name,
      source_key = selected$source_key,
      source_site = selected$source_site,
      source_lat = selected$source_lat,
      source_lon = selected$source_lon,
      era5_match_method = match_method,
      era5_distance_m = as.numeric(distance)
    )
  }))

  if (site_mapping[, any(
    era5_match_method == "nearest_coordinate" &
      era5_distance_m > ERA5_MAX_DISTANCE_M
  )]) {
    stop(
      "Nearest ", variable, " ERA5-Land match exceeds ",
      ERA5_MAX_DISTANCE_M, " m. Inspect the ERA5 join audit or increase ",
      "PHENOLOGY_BRIWECS_ERA5_MAX_DISTANCE_M explicitly."
    )
  }

  output <- merge(
    standardized,
    site_mapping[, .(Location, source_key)],
    by = "source_key",
    allow.cartesian = TRUE,
    sort = FALSE
  )[, .(Location, Date, value)]
  if (output[, anyNA(Location) | anyNA(Date)]) {
    stop(variable, " contains missing standardized location/date values.")
  }
  duplicate_rows <- output[, .N, by = .(Location, Date)][N > 1L]
  if (nrow(duplicate_rows)) {
    write_table(
      duplicate_rows,
      paste0("briwecs_duplicate_era5_", variable, "_dates.csv")
    )
    stop(
      "ERA5-Land ", variable,
      " has duplicate Location-Date rows after matching."
    )
  }
  setnames(output, "value", output_name)

  file_audit <- rbindlist(lapply(files, function(path) {
    data.table(
      variable = variable,
      year = as.integer(sub(
        paste0("^", variable, "_([0-9]{4})[.]csv$"),
        "\\1",
        basename(path)
      )),
      file = normalizePath(path),
      md5 = unname(tools::md5sum(path)),
      value_transform = value_transform
    )
  }))
  list(data = output, site_audit = site_mapping, file_audit = file_audit)
}

message(
  "Reading yearly BRIWECS ERA5-Land variables from: ",
  BRIWECS_ERA5_DIR
)
era5_results <- setNames(lapply(REQUIRED_ERA5_VARIABLES, function(variable) {
  standardize_era5_variable(variable, BRIWECS_ERA5_DIR, locations)
}), REQUIRED_ERA5_VARIABLES)

era5_file_audit <- rbindlist(lapply(
  era5_results,
  function(result) result$file_audit
))
required_era5_years <- seq.int(
  as.integer(format(
    min(as.Date(environment_targets$reference_sowing_date)) - (RECENT_LOOKBACK_DAYS - 1L),
    "%Y"
  )),
  as.integer(format(
    max(as.Date(environment_targets$reference_sowing_date)) + MAX_DAYS_AFTER_SOWING,
    "%Y"
  ))
)
required_era5_files <- CJ(
  variable = REQUIRED_ERA5_VARIABLES,
  year = required_era5_years,
  unique = TRUE
)
missing_era5_files <- required_era5_files[
  !era5_file_audit,
  on = .(variable, year)
]
if (nrow(missing_era5_files)) {
  write_table(missing_era5_files, "briwecs_missing_era5_year_files.csv")
  stop(
    "ERA5-Land variable-year files are missing for the required BRIWECS ",
    "lookback/prediction period. Inspect briwecs_missing_era5_year_files.csv."
  )
}
era5_join_audit <- rbindlist(lapply(
  era5_results,
  function(result) result$site_audit
))
setorder(era5_file_audit, variable, year)
setorder(era5_join_audit, variable, Location)
write_table(era5_file_audit, "briwecs_era5_file_audit.csv")
write_table(era5_join_audit, "briwecs_era5_join_audit.csv")

era5 <- Reduce(
  function(left, right) {
    merge(
      left,
      right,
      by = c("Location", "Date"),
      all = TRUE,
      sort = FALSE
    )
  },
  lapply(era5_results, function(result) result$data)
)
setorder(era5, Location, Date)
era5_duplicate_qc <- era5[, .N, by = .(Location, Date)][N > 1L]
if (nrow(era5_duplicate_qc)) {
  write_table(era5_duplicate_qc, "briwecs_duplicate_era5_dates.csv")
  stop("The joined ERA5-Land table has duplicate Location-Date rows.")
}
pev_values <- era5$potential_evaporation_mm
finite_pev <- pev_values[is.finite(pev_values)]
pev_value_audit <- data.table(
  n = length(pev_values),
  n_missing_or_nonfinite = sum(!is.finite(pev_values)),
  n_negative_after_conversion = sum(finite_pev < 0),
  proportion_negative_after_conversion = if (length(finite_pev)) {
    mean(finite_pev < 0)
  } else NA_real_,
  minimum_mm = if (length(finite_pev)) min(finite_pev) else NA_real_,
  q01_mm = quantile_or_na(finite_pev, 0.01),
  median_mm = if (length(finite_pev)) median(finite_pev) else NA_real_,
  q99_mm = quantile_or_na(finite_pev, 0.99),
  maximum_mm = if (length(finite_pev)) max(finite_pev) else NA_real_
)
write_table(pev_value_audit, "briwecs_era5_pev_value_audit.csv")
if (pev_value_audit$n_negative_after_conversion > 0L) {
  warning(
    pev_value_audit$n_negative_after_conversion,
    " daily potential-evaporation values are negative after conversion. ",
    "They are retained because the corrected German and Austrian training ",
    "feature scripts use potential_evaporation_mm = -1000 * pev without ",
    "clipping. Inspect briwecs_era5_pev_value_audit.csv."
  )
}
if (era5[, any(abs(potential_evaporation_mm) > 30, na.rm = TRUE)]) {
  stop(
    "Absolute daily potential evaporation exceeds 30 mm after conversion; ",
    "verify the source units."
  )
}
for (column in c("swvl1", "swvl2", "swvl3", "swvl4")) {
  if (era5[, any(
    get(column) < 0 | get(column) > 1,
    na.rm = TRUE
  )]) {
    stop(column, " must be a volumetric fraction within [0, 1].")
  }
}
trajectories <- unique(environment_targets[, .(
  environment_id, Location, SeasonYear,
  sowing_date = reference_sowing_date, sowing_DOY, lat, lon
)])
if (trajectories[, anyDuplicated(.SD), .SDcols = c("Location", "SeasonYear")]) {
  stop("Expected exactly one sowing reference and trajectory per Location x SeasonYear.")
}
setorder(trajectories, Location, SeasonYear)
trajectories[, trajectory_id := sprintf("bri_traj_%05d", .I)]
environment_targets <- merge(
  environment_targets,
  trajectories[, .(environment_id, trajectory_id)],
  by = "environment_id", all.x = TRUE, sort = FALSE
)
if (anyNA(environment_targets$trajectory_id)) stop("A target has no environment trajectory.")
write_table(environment_targets, "briwecs_environment_targets_used.csv")
write_table(trajectories, "briwecs_environment_trajectories.csv")

daily <- merge(
  trajectories,
  weather,
  by = "Location",
  allow.cartesian = TRUE,
  sort = FALSE
)
daily <- daily[
  Date >= sowing_date - (RECENT_LOOKBACK_DAYS - 1L) &
    Date <= sowing_date + MAX_DAYS_AFTER_SOWING
]
daily <- merge(
  daily,
  era5,
  by = c("Location", "Date"),
  all.x = TRUE,
  sort = FALSE
)
daily <- merge(
  daily,
  soil_join_audit[, .(
    Location,
    paw_fc33_0_100_mm,
    wilting_point_0_100_mm
  )],
  by = "Location",
  all.x = TRUE,
  sort = FALSE
)
setorder(daily, trajectory_id, Date)

coverage_audit <- daily[, .(
  n_days = .N,
  first_date = min(Date),
  last_date = max(Date),
  expected_first_date = first(sowing_date) - (RECENT_LOOKBACK_DAYS - 1L),
  expected_last_date = first(sowing_date) + MAX_DAYS_AFTER_SOWING,
  n_calendar_gaps = sum(diff(as.integer(Date)) != 1L),
  n_missing_tmean = sum(is.na(Tmean)),
  n_missing_tmin = sum(is.na(TN)),
  n_missing_tmax = sum(is.na(TX)),
  n_missing_precipitation = sum(is.na(RR)),
  n_missing_radiation = sum(is.na(Rad)),
  n_missing_pevap = sum(is.na(potential_evaporation_mm)),
  n_missing_swvl1 = sum(is.na(swvl1)),
  n_missing_swvl2 = sum(is.na(swvl2)),
  n_missing_swvl3 = sum(is.na(swvl3)),
  n_missing_swvl4 = sum(is.na(swvl4))
), by = .(trajectory_id, Location, SeasonYear, sowing_date)]
write_table(coverage_audit, "briwecs_daily_weather_coverage.csv")
if (nrow(coverage_audit) != nrow(trajectories)) {
  stop("At least one BRIWECS trajectory has no daily weather rows.")
}
if (coverage_audit[, any(first_date != expected_first_date)]) {
  stop(
    "Daily inputs do not provide all 29 pre-sowing lookback days for every ",
    "trajectory."
  )
}
if (coverage_audit[, any(last_date != expected_last_date)]) {
  stop(
    "Daily weather does not cover the complete fixed prediction horizon for ",
    "every trajectory. Shorten PHENOLOGY_BRIWECS_MAX_DAYS_AFTER_SOWING only ",
    "with a documented, outcome-independent rule."
  )
}
if (coverage_audit[, any(n_calendar_gaps > 0L)]) {
  stop("Daily input gaps were detected; rolling predictors cannot be reproduced.")
}
if (coverage_audit[, any(
  n_missing_tmean > 0L | n_missing_tmin > 0L | n_missing_tmax > 0L
)]) {
  stop("Missing daily temperature values prevent exact feature reproduction.")
}
if (coverage_audit[, any(
  n_missing_precipitation > 0L | n_missing_radiation > 0L
)]) {
  stop(
    "Missing precipitation or radiation values prevent exact rolling-window ",
    "feature reproduction."
  )
}
if (coverage_audit[, any(
  n_missing_pevap > 0L | n_missing_swvl1 > 0L |
    n_missing_swvl2 > 0L | n_missing_swvl3 > 0L |
    n_missing_swvl4 > 0L
)]) {
  stop(
    "Missing ERA5-Land radiation, potential evaporation or soil-water values ",
    "prevent exact feature reproduction."
  )
}

daily[, days_since_sowing := as.integer(Date - sowing_date)]
daily[, CDD := pmax(COLD_THRESHOLD - TN, 0)]
daily[, HDD := pmax(TX - HEAT_THRESHOLD, 0)]
daily[, `:=`(
  soil_water_0_28_vwc = (
    SWVL1_DEPTH_MM * swvl1 + SWVL2_DEPTH_MM * swvl2
  ) / SW_DEPTH_0_28_MM,
  soil_water_0_100_mm =
    SWVL1_DEPTH_MM * swvl1 +
      SWVL2_DEPTH_MM * swvl2 +
      SWVL3_DEPTH_MM * swvl3,
  soil_water_0_100_vwc = (
    SWVL1_DEPTH_MM * swvl1 +
      SWVL2_DEPTH_MM * swvl2 +
      SWVL3_DEPTH_MM * swvl3
  ) / SW_DEPTH_0_100_MM
)]
daily[, available_soil_water_0_100_mm := pmax(
  soil_water_0_100_mm - wilting_point_0_100_mm,
  0
)]
daily[, FTSW_0_100_raw :=
  available_soil_water_0_100_mm / paw_fc33_0_100_mm]
daily[, FTSW_0_100_proxy := pmin(1, pmax(0, FTSW_0_100_raw))]
ftsw_audit <- daily[, .(
  n_days = .N,
  raw_minimum = min(FTSW_0_100_raw),
  raw_median = median(FTSW_0_100_raw),
  raw_maximum = max(FTSW_0_100_raw),
  n_raw_below_zero = sum(FTSW_0_100_raw < 0),
  n_raw_above_one = sum(FTSW_0_100_raw > 1),
  proxy_minimum = min(FTSW_0_100_proxy),
  proxy_maximum = max(FTSW_0_100_proxy)
), by = .(Location, SeasonYear, trajectory_id)]
write_table(ftsw_audit, "briwecs_ftsw_0_100_audit.csv")

# These windows include the current day and the preceding width - 1 calendar
# days, exactly as in the corrected calculate predictors scritps. All recent-weather
# predictors use the pre-specified 30-day window.
daily[, `:=`(
  sumRR_30 = rolling_sum_complete(RR, Date, 30L),
  sumRad_30 = rolling_sum_complete(Rad, Date, 30L),
  sumCDD_30 = rolling_sum_complete(CDD, Date, 30L),
  sumHDD_30 = rolling_sum_complete(HDD, Date, 30L),
  sumPEvap_30 = rolling_sum_complete(
    potential_evaporation_mm, Date, 30L
  ),
  meanSoilWater_0_100_30 = rolling_mean_complete(
    soil_water_0_100_vwc, Date, 30L
  )
), by = trajectory_id]
daily[, climaticWaterBalance_30 := sumRR_30 - sumPEvap_30]

# Remove the pre-sowing lookback only after the retrospective features have
# been formed. Process states then start at sowing and never use observed BBCH.
daily <- daily[days_since_sowing >= 0L]
setorder(daily, trajectory_id, Date)

# Cumulative weather features start at sowing. Pre-sowing rows are used only to
# construct retrospective 30-day windows and never enter cumulative sums.
daily[, `:=`(
  cumRR = cumsum_complete_to_date(RR),
  cumRad = cumsum_complete_to_date(Rad),
  cumCDD = cumsum_complete_to_date(CDD),
  cumHDD = cumsum_complete_to_date(HDD)
), by = trajectory_id]

daily[, photoperiod_civil := calc_civil_photoperiod(lat, Date)]
daily[, tempfun_beta := fcase(
  Tmean <= T_BASE | Tmean >= T_CEILING, 0,
  default =
    ((T_CEILING - Tmean) / (T_CEILING - T_OPTIMUM)) *
      ((Tmean - T_BASE) / (T_OPTIMUM - T_BASE)) ^
      ((T_OPTIMUM - T_BASE) / (T_CEILING - T_OPTIMUM))
)]
daily[, DTU_beta := (T_OPTIMUM - T_BASE) * tempfun_beta]
daily[, ppfun_raw := fifelse(
  photoperiod_civil < CPP,
  1 - PPSEN * (CPP - photoperiod_civil)^2,
  1
)]
daily[, ppfun := pmin(1, pmax(0, ppfun_raw))]
daily[, VERDAY := fcase(
  Tmean <= TBV, 0,
  Tmean < TP1V, (Tmean - TBV) / (TP1V - TBV),
  Tmean <= TP2V, 1,
  Tmean < TCV, (TCV - Tmean) / (TCV - TP2V),
  default = 0
)]
daily[, cumVER := accumulate_vernalization(VERDAY, TX), by = trajectory_id]
daily[, verfun := pmin(1, pmax(0, 1 - VSEN * (VDSAT - cumVER)))]
phase_columns <- c("ppfun_phase", "verfun_phase", "devfun_phase")
daily[, (phase_columns) := calc_soltani_phase_state(
  tempfun_beta,
  ppfun,
  verfun
), by = trajectory_id]
daily[, `:=`(
  cumdevfun_phase = cumsum_complete_to_date(devfun_phase),
  cumDTU_beta = cumsum_complete_to_date(DTU_beta)
), by = trajectory_id]

# Make the complete external feature contract visible before selecting one
# representation for prediction.
missing_daily_features <- setdiff(ALL_EXTERNAL_FEATURES, names(daily))
if (length(missing_daily_features)) {
  stop(
    "BRIWECS daily table is missing required predictor(s): ",
    paste(missing_daily_features, collapse = ", ")
  )
}
all_feature_audit <- rbindlist(lapply(
  ALL_EXTERNAL_FEATURES,
  function(feature) {
    values <- daily[[feature]]
    data.table(
      feature = feature,
      n = length(values),
      n_missing = sum(is.na(values)),
      minimum = if (all(is.na(values))) NA_real_ else min(values, na.rm = TRUE),
      median = if (all(is.na(values))) NA_real_ else median(values, na.rm = TRUE),
      maximum = if (all(is.na(values))) NA_real_ else max(values, na.rm = TRUE)
    )
  }
))
all_feature_audit[, selected_for_model := feature %in% FEATURES]
write_table(all_feature_audit, "briwecs_all_feature_audit.csv")

if (anyNA(daily[, ..ALL_EXTERNAL_FEATURES]) || any(!vapply(
  daily[, ..ALL_EXTERNAL_FEATURES],
  function(x) is.numeric(x) && all(is.finite(x)),
  logical(1)
))) {
  stop("One or more calculated BRIWECS predictors are missing/non-finite.")
}

# The fitted learner retains the TaskRegrST coordinate roles. Coordinates are
# not XGBoost features, but predict_newdata() needs them to reconstruct a valid
# spatiotemporal prediction task from the stored training-task schema.
daily[, `:=`(
  `.lon` = as.numeric(lon),
  `.lat` = as.numeric(lat)
)]
if (daily[, anyNA(.lon) | anyNA(.lat) |
    any(!is.finite(.lon)) | any(!is.finite(.lat))]) {
  stop("BRIWECS prediction coordinates contain missing/non-finite values.")
}

atomic_save_rds(
  daily,
  file.path(OUTPUT_DIR, "briwecs_daily_predictor_table.rds"),
  compress = FALSE
)

# ------------------------------------------------------------------------------
# 7. RAW EXTERNAL PREDICTIONS AND MILESTONE MATCHING
# ------------------------------------------------------------------------------

message("Predicting the complete independent BRIWECS daily trajectories.")
prediction_task_columns <- c(MODEL_FEATURE_ORDER, ".lon", ".lat")
daily_prediction <- final_learner$predict_newdata(
  newdata = daily[, ..prediction_task_columns]
)
daily[, xgboost_raw := as.numeric(daily_prediction$response)]
if (any(!is.finite(daily$xgboost_raw))) {
  stop("Final XGBoost produced missing or non-finite BRIWECS predictions.")
}
daily[, xgboost_monotone := cummax(xgboost_raw), by = trajectory_id]
if (RUN_GAM_INTERACTION_ABLATION) {
  gam_prediction_data <- daily[, ..GAM_ABLATION_FEATURES]
  daily[, gam_additive_raw := predict_gam_ablation(
    gam_additive_full,
    gam_prediction_data,
    "Full-development additive GAM"
  )]
  daily[, gam_interaction_raw := predict_gam_ablation(
    gam_interaction_full,
    gam_prediction_data,
    "Full-development interaction GAM"
  )]
  daily[, `:=`(
    gam_additive_monotone = cummax(gam_additive_raw),
    gam_interaction_monotone = cummax(gam_interaction_raw)
  ), by = trajectory_id]
}
daily[, gdd_threshold := predict_thermal_baseline(
  cumDTU_beta,
  thermal_baseline_thresholds
)]

# The target generator has already produced one observation per environment and
# milestone. Never melt the raw plots or recompute a mean after taking losses.
milestones <- copy(environment_targets)

daily_prediction_columns <- c(
  "trajectory_id", "Date", FEATURES,
  "xgboost_raw", "xgboost_monotone",
  if (RUN_GAM_INTERACTION_ABLATION) {
    c(
      "gam_additive_raw", "gam_interaction_raw",
      "gam_additive_monotone", "gam_interaction_monotone"
    )
  },
  "gdd_threshold"
)
milestone_id_columns <- setdiff(names(milestones), FEATURES)
milestone_eval <- merge(
  milestones[, ..milestone_id_columns],
  daily[, ..daily_prediction_columns],
  by.x = c("trajectory_id", "observed_date"),
  by.y = c("trajectory_id", "Date"),
  all.x = TRUE,
  sort = FALSE
)

missing_milestone_predictions <- milestone_eval[is.na(xgboost_raw)]
if (nrow(missing_milestone_predictions)) {
  write_table(
    missing_milestone_predictions[, .(
      environment_milestone_id,
      Location,
      SeasonYear,
      milestone,
      observed_date,
      trajectory_id
    )],
    "briwecs_missing_milestone_predictions.csv"
  )
  stop(
    "Some observed milestones fall outside the independent weather horizon. ",
    "See briwecs_missing_milestone_predictions.csv."
  )
}

prediction_measure_columns <- c(
  xgboost = "xgboost_raw",
  if (RUN_GAM_INTERACTION_ABLATION) {
    c(
      gam_additive = "gam_additive_raw",
      gam_interaction = "gam_interaction_raw"
    )
  },
  if (INCLUDE_MONOTONE_SENSITIVITY) {
    c(xgboost_monotone_sensitivity = "xgboost_monotone")
  },
  if (INCLUDE_GDD_BASELINE) c(gdd_threshold = "gdd_threshold")
)
validation_long <- melt(
  milestone_eval,
  id.vars = setdiff(names(milestone_eval), unname(prediction_measure_columns)),
  measure.vars = unname(prediction_measure_columns),
  variable.name = "prediction_column",
  value.name = "prediction",
  variable.factor = FALSE
)
prediction_lookup <- setNames(names(prediction_measure_columns), prediction_measure_columns)
validation_long[, learner_id := unname(prediction_lookup[prediction_column])]
validation_long[, prediction_column := NULL]
validation_long[, `:=`(
  error = prediction - BBCH_obs,
  abs_error = abs(prediction - BBCH_obs),
  squared_error = (prediction - BBCH_obs)^2,
  outside_0_99 = prediction < 0 | prediction > 99
)]

setcolorder(
  validation_long,
  c(
    "learner_id", "environment_milestone_id", "Location", "SeasonYear",
    "environment_id", "trajectory_id", "milestone", "BBCH_obs",
    "observed_date", "prediction", "error", "abs_error",
    "squared_error", "outside_0_99",
    setdiff(names(validation_long), c(
      "learner_id", "environment_milestone_id", "Location", "SeasonYear",
      "environment_id", "trajectory_id", "milestone", "BBCH_obs",
      "observed_date", "prediction", "error", "abs_error",
      "squared_error", "outside_0_99"
    ))
  )
)
write_table(validation_long, "briwecs_external_predictions_long.csv")

# ------------------------------------------------------------------------------
# 8. ENVIRONMENT-TARGET METRICS, UNCERTAINTY AND TARGET SUPPORT
# ------------------------------------------------------------------------------

environment_stage_units <- make_environment_stage_units(validation_long)
write_table(
  environment_stage_units,
  "briwecs_environment_milestone_units.csv"
)

macro_metrics <- function(data, grouping_columns = character()) {
  data[, .(
    target_definition = first(target_definition),
    target_input_md5 = first(target_input_md5),
    reference_treatment = first(reference_treatment),
    unit_weighting = "equal_environment_milestone_units",
    n_environment_milestone_units = .N,
    n_environments = uniqueN(environment_id),
    n_locations = uniqueN(Location),
    macro_mae = mean(mean_abs_error),
    macro_rmse = sqrt(mean(mean_squared_error)),
    macro_bias = mean(mean_error),
    median_unit_absolute_error = median(mean_abs_error),
    within_3 = mean(proportion_within_3),
    within_5 = mean(proportion_within_5),
    within_10 = mean(proportion_within_10),
    proportion_outside_0_99 = mean(proportion_outside_0_99)
  ), by = grouping_columns]
}

primary_metrics_overall <- macro_metrics(
  environment_stage_units,
  "learner_id"
)
primary_metrics_overall[, milestone := "Overall"]
primary_metrics_by_stage <- macro_metrics(
  environment_stage_units,
  c("learner_id", "milestone")
)
primary_metrics <- rbindlist(
  list(primary_metrics_overall, primary_metrics_by_stage),
  use.names = TRUE,
  fill = TRUE
)
setcolorder(primary_metrics, c(
  "learner_id", "milestone",
  setdiff(names(primary_metrics), c("learner_id", "milestone"))
))
write_table(primary_metrics, "briwecs_primary_macro_performance.csv")

# Paired interaction ablation. Negative alternative-minus-reference deltas
# favour adding the selected GAM interaction for prediction. The bootstrap resamples
# whole environments and therefore preserves all stages within an environment.
gam_ablation_paired_units <- data.table()
gam_ablation_summary <- data.table()

pair_gam_ablation_units <- function(long_units) {
  key_columns <- c(
    "dataset", "scheme", "unit_id", "environment_id", "milestone"
  )
  duplicate_check <- long_units[, .N, by = c(key_columns, "model")][N != 1L]
  if (nrow(duplicate_check)) {
    stop("GAM ablation units are not unique by dataset, unit and model.")
  }
  paired <- dcast(
    long_units,
    dataset + scheme + unit_id + environment_id + milestone ~ model,
    value.var = c(
      "mean_abs_error", "mean_squared_error", "mean_error",
      "mean_prediction"
    )
  )
  required <- unlist(lapply(
    c(
      "mean_abs_error", "mean_squared_error", "mean_error",
      "mean_prediction"
    ),
    function(metric) paste0(metric, c("_gam_additive", "_gam_interaction"))
  ))
  missing <- setdiff(required, names(paired))
  if (length(missing)) {
    stop(
      "The paired GAM ablation table is missing: ",
      paste(missing, collapse = ", ")
    )
  }
  if (paired[, anyNA(.SD), .SDcols = required]) {
    stop("Additive and interaction GAMs do not have identical paired units.")
  }
  paired[, `:=`(
    delta_mae_interaction_minus_additive =
      mean_abs_error_gam_interaction - mean_abs_error_gam_additive,
    delta_mse_interaction_minus_additive =
      mean_squared_error_gam_interaction - mean_squared_error_gam_additive,
    delta_bias_interaction_minus_additive =
      mean_error_gam_interaction - mean_error_gam_additive
  )]
  paired
}

summarize_gam_ablation_group <- function(data, replicates, seed) {
  calculate_metrics <- function(weights = rep(1, nrow(data))) {
    reference_mae <- weighted.mean(
      data$mean_abs_error_gam_additive,
      weights
    )
    alternative_mae <- weighted.mean(
      data$mean_abs_error_gam_interaction,
      weights
    )
    reference_mse <- weighted.mean(
      data$mean_squared_error_gam_additive,
      weights
    )
    alternative_mse <- weighted.mean(
      data$mean_squared_error_gam_interaction,
      weights
    )
    reference_bias <- weighted.mean(data$mean_error_gam_additive, weights)
    alternative_bias <- weighted.mean(
      data$mean_error_gam_interaction,
      weights
    )
    c(
      macro_mae_reference = reference_mae,
      macro_mae_alternative = alternative_mae,
      delta_macro_mae = alternative_mae - reference_mae,
      macro_rmse_reference = sqrt(reference_mse),
      macro_rmse_alternative = sqrt(alternative_mse),
      delta_macro_rmse = sqrt(alternative_mse) - sqrt(reference_mse),
      macro_bias_reference = reference_bias,
      macro_bias_alternative = alternative_bias,
      delta_macro_bias = alternative_bias - reference_bias
    )
  }

  point <- calculate_metrics()
  clusters <- unique(as.character(data$environment_id))
  delta_draws <- matrix(NA_real_, nrow = replicates, ncol = 3L)
  colnames(delta_draws) <- c(
    "delta_macro_mae", "delta_macro_rmse", "delta_macro_bias"
  )
  if (length(clusters) >= 2L) {
    cluster_index <- match(as.character(data$environment_id), clusters)
    set.seed(seed)
    for (b in seq_len(replicates)) {
      cluster_counts <- tabulate(
        sample.int(length(clusters), length(clusters), replace = TRUE),
        nbins = length(clusters)
      )
      draw <- calculate_metrics(cluster_counts[cluster_index])
      delta_draws[b, ] <- draw[colnames(delta_draws)]
    }
  }

  interval <- function(metric, probability) {
    values <- delta_draws[, metric]
    if (!any(is.finite(values))) return(NA_real_)
    quantile_or_na(values, probability)
  }
  data.table(
    reference = "gam_additive",
    alternative = "gam_interaction",
    n_units = nrow(data),
    n_environments = length(clusters),
    macro_mae_reference = unname(point["macro_mae_reference"]),
    macro_mae_alternative = unname(point["macro_mae_alternative"]),
    delta_macro_mae = unname(point["delta_macro_mae"]),
    delta_macro_mae_ci_low = interval("delta_macro_mae", 0.025),
    delta_macro_mae_ci_high = interval("delta_macro_mae", 0.975),
    macro_rmse_reference = unname(point["macro_rmse_reference"]),
    macro_rmse_alternative = unname(point["macro_rmse_alternative"]),
    delta_macro_rmse = unname(point["delta_macro_rmse"]),
    delta_macro_rmse_ci_low = interval("delta_macro_rmse", 0.025),
    delta_macro_rmse_ci_high = interval("delta_macro_rmse", 0.975),
    macro_bias_reference = unname(point["macro_bias_reference"]),
    macro_bias_alternative = unname(point["macro_bias_alternative"]),
    delta_macro_bias = unname(point["delta_macro_bias"]),
    delta_macro_bias_ci_low = interval("delta_macro_bias", 0.025),
    delta_macro_bias_ci_high = interval("delta_macro_bias", 0.975),
    proportion_units_improved = mean(
      data$mean_abs_error_gam_interaction <
        data$mean_abs_error_gam_additive
    ),
    rmse_supports_interaction = unname(point["delta_macro_rmse"]) < 0,
    rmse_ci_excludes_zero = is.finite(
      interval("delta_macro_rmse", 0.975)
    ) && interval("delta_macro_rmse", 0.975) < 0
  )
}

if (RUN_GAM_INTERACTION_ABLATION) {
  internal_ablation_long <- internal_gam_ablation_units[, .(
    dataset,
    scheme,
    unit_id = paste(split_id, environment_id, milestone, sep = "__"),
    environment_id,
    milestone,
    model,
    mean_abs_error,
    mean_squared_error,
    mean_error,
    mean_prediction
  )]
  external_ablation_long <- environment_stage_units[
    learner_id %in% c("gam_additive", "gam_interaction"),
    .(
      dataset = "External BRIWECS",
      scheme = "external_briwecs",
      unit_id = paste(environment_id, milestone, sep = "__"),
      environment_id,
      milestone,
      model = learner_id,
      mean_abs_error,
      mean_squared_error,
      mean_error,
      mean_prediction
    )
  ]
  gam_ablation_paired_units <- pair_gam_ablation_units(rbindlist(
    list(internal_ablation_long, external_ablation_long),
    use.names = TRUE
  ))
  gam_ablation_paired_units <- annotate_gam_ablation(
    gam_ablation_paired_units, GAM_ABLATION_SPEC,
    REFERENCE_TREATMENT, TARGET_DEFINITION, target_input_md5
  )
  write_table(
    gam_ablation_paired_units,
    "interaction_ablation_paired_environment_stage_units.csv"
  )

  gam_ablation_summary_input <- rbindlist(list(
    copy(gam_ablation_paired_units)[, milestone := "Overall"],
    gam_ablation_paired_units
  ))
  ablation_groups <- unique(gam_ablation_summary_input[, .(
    dataset,
    scheme,
    milestone
  )])
  setorder(ablation_groups, dataset, scheme, milestone)
  gam_ablation_summary <- rbindlist(lapply(
    seq_len(nrow(ablation_groups)),
    function(i) {
      group <- ablation_groups[i]
      subset <- gam_ablation_summary_input[
        dataset == group$dataset &
          scheme == group$scheme &
          milestone == group$milestone
      ]
      result <- summarize_gam_ablation_group(
        subset,
        BOOTSTRAP_REPLICATES,
        SEED + 26000L + i
      )
      result[, `:=`(
        dataset = group$dataset,
        scheme = group$scheme,
        milestone = group$milestone
      )]
      result
    }
  ))
  gam_ablation_summary <- annotate_gam_ablation(
    gam_ablation_summary, GAM_ABLATION_SPEC,
    REFERENCE_TREATMENT, TARGET_DEFINITION, target_input_md5
  )
  setcolorder(
    gam_ablation_summary,
    c(
      "dataset", "scheme", "milestone", "reference", "alternative",
      setdiff(
        names(gam_ablation_summary),
        c("dataset", "scheme", "milestone", "reference", "alternative")
      )
    )
  )
  write_table(
    gam_ablation_summary,
    "interaction_ablation_performance_summary.csv"
  )
}

# No replicate-weighted summary: plot counts are target-support metadata only.
# Bootstrap intervals below condition on the supplied target estimates; they
# do not propagate uncertainty from estimating the cultivar-adjusted dates.

bootstrap_metric <- function(unit_data, replicates, seed) {
  clusters <- unique(unit_data$environment_id)
  if (length(clusters) < 2L) {
    return(data.table(
      metric = c("macro_mae", "macro_rmse", "macro_bias"),
      estimate = c(
        mean(unit_data$mean_abs_error),
        sqrt(mean(unit_data$mean_squared_error)),
        mean(unit_data$mean_error)
      ),
      ci_low = NA_real_,
      ci_high = NA_real_,
      n_clusters = length(clusters)
    ))
  }
  set.seed(seed)
  values <- matrix(NA_real_, nrow = replicates, ncol = 3L)
  for (b in seq_len(replicates)) {
    sampled <- sample(clusters, length(clusters), replace = TRUE)
    sampled_data <- rbindlist(lapply(seq_along(sampled), function(i) {
      copy(unit_data[environment_id == sampled[i]])
    }))
    values[b, ] <- c(
      mean(sampled_data$mean_abs_error),
      sqrt(mean(sampled_data$mean_squared_error)),
      mean(sampled_data$mean_error)
    )
  }
  estimates <- c(
    mean(unit_data$mean_abs_error),
    sqrt(mean(unit_data$mean_squared_error)),
    mean(unit_data$mean_error)
  )
  data.table(
    metric = c("macro_mae", "macro_rmse", "macro_bias"),
    estimate = estimates,
    ci_low = apply(values, 2L, quantile, 0.025, names = FALSE),
    ci_high = apply(values, 2L, quantile, 0.975, names = FALSE),
    n_clusters = length(clusters)
  )
}

bootstrap_groups <- unique(environment_stage_units[, .(learner_id, milestone)])
bootstrap_groups <- rbind(
  data.table(
    learner_id = unique(environment_stage_units$learner_id),
    milestone = "Overall"
  ),
  bootstrap_groups
)
bootstrap_results <- rbindlist(lapply(seq_len(nrow(bootstrap_groups)), function(i) {
  learner <- bootstrap_groups$learner_id[i]
  stage <- bootstrap_groups$milestone[i]
  subset <- environment_stage_units[learner_id == learner]
  if (stage != "Overall") subset <- subset[milestone == stage]
  result <- bootstrap_metric(
    subset,
    BOOTSTRAP_REPLICATES,
    SEED + 20000L + i
  )
  result[, `:=`(learner_id = learner, milestone = stage)]
  result
}))
bootstrap_results[, `:=`(
  target_definition = TARGET_DEFINITION,
  target_input_md5 = ..target_input_md5,
  reference_treatment = REFERENCE_TREATMENT,
  target_estimation_uncertainty_propagated = FALSE
)]
setcolorder(bootstrap_results, c(
  "learner_id", "milestone", "metric", "estimate", "ci_low", "ci_high",
  "n_clusters"
))
write_table(
  bootstrap_results,
  "briwecs_cluster_bootstrap_performance.csv"
)

# Paired unit-level contrast against the conventional thermal baseline.
paired_units <- dcast(
  environment_stage_units[
    learner_id %in% c("xgboost", "gdd_threshold")
  ],
  environment_id + Location + SeasonYear + milestone ~ learner_id,
  value.var = c("mean_abs_error", "mean_squared_error", "mean_error")
)
if (nrow(paired_units) &&
    all(c("mean_abs_error_xgboost", "mean_abs_error_gdd_threshold") %in%
      names(paired_units))) {
  paired_units[, `:=`(
    delta_mae = mean_abs_error_xgboost - mean_abs_error_gdd_threshold,
    delta_mse = mean_squared_error_xgboost - mean_squared_error_gdd_threshold,
    delta_bias = mean_error_xgboost - mean_error_gdd_threshold
  )]
  paired_units <- merge(
    paired_units,
    environment_targets[, .(
      environment_id, Location, SeasonYear, milestone, reference_treatment,
      reference_sowing_date, observed_date, observed_environment_days_after_sowing,
      target_definition, target_input_md5, target_source
    )],
    by = c("environment_id", "Location", "SeasonYear", "milestone"),
    all.x = TRUE, sort = FALSE
  )
  write_table(paired_units, "briwecs_xgboost_vs_gdd_paired_units.csv")
}

target_method_support <- environment_targets[, .(
  n_environment_milestone_targets = .N,
  n_environments = uniqueN(environment_id),
  n_singular_adjustments = sum(adjustment_singular %in% TRUE),
  min_contributing_plots = min(n_plots),
  median_contributing_plots = as.numeric(median(n_plots)),
  max_contributing_plots = max(n_plots),
  min_contributing_genotypes = min(n_genotypes),
  median_contributing_genotypes = as.numeric(median(n_genotypes)),
  max_contributing_genotypes = max(n_genotypes),
  median_target_se_days = quantile_or_na(target_standard_error_days, 0.5),
  max_abs_date_rounding_days = max(abs(target_date_rounding_days))
), by = .(reference_treatment, milestone, target_source)]
write_table(target_method_support, "briwecs_target_method_support.csv")
target_method_metrics <- macro_metrics(
  environment_stage_units, c("learner_id", "milestone", "target_source")
)
write_table(target_method_metrics, "briwecs_performance_by_target_method.csv")

environment_metrics <- environment_stage_units[, .(
  target_definition = first(target_definition),
  reference_treatment = first(reference_treatment),
  n_milestones = .N,
  macro_mae = mean(mean_abs_error),
  macro_rmse = sqrt(mean(mean_squared_error)),
  macro_bias = mean(mean_error)
), by = .(learner_id, environment_id, Location, SeasonYear)]
write_table(environment_metrics, "briwecs_performance_by_environment.csv")

# ------------------------------------------------------------------------------
# 9. DATE-OFFSET ANALYSIS FROM DAILY MONOTONE TRAJECTORIES
# ------------------------------------------------------------------------------

first_crossing_date <- function(date, prediction, threshold) {
  hit <- which(is.finite(prediction) & prediction >= threshold)
  if (!length(hit)) return(as.IDate(NA))
  date[min(hit)]
}

threshold_models <- c(
  xgboost = "xgboost_monotone",
  if (RUN_GAM_INTERACTION_ABLATION) {
    c(
      gam_additive = "gam_additive_monotone",
      gam_interaction = "gam_interaction_monotone"
    )
  },
  if (INCLUDE_GDD_BASELINE) c(gdd_threshold = "gdd_threshold")
)
threshold_dates <- rbindlist(lapply(names(threshold_models), function(model) {
  column <- threshold_models[[model]]
  rbindlist(lapply(c(10, 59, 87), function(stage) {
    daily[, .(
      predicted_date = first_crossing_date(Date, get(column), stage)
    ), by = .(trajectory_id, Location, SeasonYear)][
      , `:=`(learner_id = model, BBCH_obs = stage)
    ]
  }))
}))

offset_validation <- make_environment_date_offsets(milestones, threshold_dates)
write_table(offset_validation, "briwecs_milestone_date_offsets.csv")

offset_environment_units <- copy(offset_validation)
offset_environment_units[, `:=`(
  n_target_rows = 1L,
  threshold_coverage = as.numeric(threshold_reached),
  mean_abs_offset_days = abs_offset_days,
  mean_squared_offset_days = squared_offset_days,
  mean_offset_days = offset_days
)]
write_table(
  offset_environment_units,
  "briwecs_environment_milestone_date_offset_units.csv"
)

offset_summary <- offset_environment_units[, .(
  target_definition = first(target_definition),
  target_input_md5 = first(target_input_md5),
  reference_treatment = first(reference_treatment),
  offset_definition = "predicted_DAS_minus_unrounded_target_DAS",
  n_environment_milestone_units = .N,
  n_environments = uniqueN(environment_id),
  threshold_coverage = mean(threshold_coverage),
  n_units_reached = sum(is.finite(mean_abs_offset_days)),
  median_absolute_offset_days = if (any(is.finite(mean_abs_offset_days))) {
    median(mean_abs_offset_days, na.rm = TRUE)
  } else NA_real_,
  mae_days = if (any(is.finite(mean_abs_offset_days))) {
    mean(mean_abs_offset_days, na.rm = TRUE)
  } else NA_real_,
  rmse_days = if (any(is.finite(mean_squared_offset_days))) {
    sqrt(mean(mean_squared_offset_days, na.rm = TRUE))
  } else NA_real_,
  bias_days = if (any(is.finite(mean_offset_days))) {
    mean(mean_offset_days, na.rm = TRUE)
  } else NA_real_
), by = .(learner_id, milestone, BBCH_obs)]
write_table(offset_summary, "briwecs_milestone_date_offset_summary.csv")

# ------------------------------------------------------------------------------
# 10. EXTRAPOLATION/SUPPORT AUDIT AT THE OBSERVED MILESTONES
# ------------------------------------------------------------------------------

support_rows <- milestone_eval[, c(
  "environment_milestone_id", "environment_id", "milestone",
  "trajectory_id", "observed_date", FEATURES
), with = FALSE]
validation_feature_support <- rbindlist(lapply(FEATURES, function(feature_name) {
  x <- support_rows[[feature_name]]
  train_row <- training_feature_support[
    get("feature") == feature_name
  ]
  data.table(
    feature = feature_name,
    evaluation_unit = "environment_milestone_target",
    n_validation = length(x),
    validation_min = min(x),
    validation_q01 = quantile_or_na(x, 0.01),
    validation_median = median(x),
    validation_q99 = quantile_or_na(x, 0.99),
    validation_max = max(x),
    train_min = train_row$train_min,
    train_q01 = train_row$train_q01,
    train_median = train_row$train_median,
    train_q99 = train_row$train_q99,
    train_max = train_row$train_max,
    n_outside_training_range = sum(
      x < train_row$train_min | x > train_row$train_max
    ),
    proportion_outside_training_range = mean(
      x < train_row$train_min | x > train_row$train_max
    ),
    proportion_outside_training_q01_q99 = mean(
      x < train_row$train_q01 | x > train_row$train_q99
    )
  )
}))
write_table(
  validation_feature_support,
  "briwecs_held_out_feature_support_audit.csv"
)

# ------------------------------------------------------------------------------
# 11. OPTIONAL PAIRED COMPARISON WITH THE FROZEN CUMULATIVE MODEL
# ------------------------------------------------------------------------------

cumulative_unit_path <- file.path(
  CUMULATIVE_VALIDATION_DIR,
  "briwecs_environment_milestone_units.csv"
)
cumulative_comparison_available <- FALSE
cumulative_comparison_reason <- if (FEATURE_SET_SHORT != "hybrid") {
  "Not applicable to the cumulative run."
} else if (!file.exists(cumulative_unit_path)) {
  "No completed cumulative environment-target unit table was found."
} else {
  NA_character_
}
cumulative_units <- NULL
hybrid_units <- copy(environment_stage_units[learner_id == "xgboost"])
if (FEATURE_SET_SHORT == "hybrid" && file.exists(cumulative_unit_path)) {
  cumulative_units <- fread(cumulative_unit_path)
  if (!"learner_id" %in% names(cumulative_units)) {
    cumulative_comparison_reason <- "Cumulative unit table has no learner_id."
  } else {
    cumulative_units <- cumulative_units[learner_id == "xgboost"]
    reason <- cumulative_target_compatibility(hybrid_units, cumulative_units)
    cumulative_comparison_available <- is.null(reason)
    cumulative_comparison_reason <- if (is.null(reason)) "Compatible paper-style targets." else reason
  }
}
write_table(data.table(
  comparison_available = cumulative_comparison_available,
  reason = cumulative_comparison_reason,
  cumulative_unit_file = cumulative_unit_path,
  required_target_definition = TARGET_DEFINITION,
  current_target_input_md5 = target_input_md5
), "briwecs_hybrid_vs_cumulative_comparison_status.csv")

if (cumulative_comparison_available) {
  message("Creating a paired comparison on identical paper-style targets from: ", cumulative_unit_path)
  comparison_keys <- c(
    "environment_id", "Location", "SeasonYear", "milestone"
  )
  comparison_measures <- c(
    "mean_abs_error", "mean_squared_error", "mean_error"
  )
  missing_cumulative_columns <- setdiff(
    c("learner_id", comparison_keys, comparison_measures),
    names(cumulative_units)
  )
  if (length(missing_cumulative_columns)) {
    stop(
      "The cumulative BRIWECS environment-unit table is missing: ",
      paste(missing_cumulative_columns, collapse = ", ")
    )
  }

  hybrid_compare <- hybrid_units[, c(
    comparison_keys,
    comparison_measures
  ), with = FALSE]
  cumulative_compare <- cumulative_units[, c(
    comparison_keys,
    comparison_measures
  ), with = FALSE]
  setnames(
    hybrid_compare,
    comparison_measures,
    paste0(comparison_measures, "_hybrid")
  )
  setnames(
    cumulative_compare,
    comparison_measures,
    paste0(comparison_measures, "_cumulative")
  )
  paired_feature_units <- merge(
    hybrid_compare,
    cumulative_compare,
    by = comparison_keys,
    all = TRUE,
    sort = FALSE
  )
  comparison_value_columns <- c(
    paste0(comparison_measures, "_hybrid"),
    paste0(comparison_measures, "_cumulative")
  )
  if (paired_feature_units[, anyNA(.SD), .SDcols = comparison_value_columns]) {
    stop(
      "The hybrid and cumulative BRIWECS runs do not contain identical ",
      "environment x milestone units; a paired comparison would be invalid."
    )
  }
  paired_feature_units[, `:=`(
    delta_mae_hybrid_minus_cumulative =
      mean_abs_error_hybrid - mean_abs_error_cumulative,
    delta_mse_hybrid_minus_cumulative =
      mean_squared_error_hybrid - mean_squared_error_cumulative,
    delta_bias_hybrid_minus_cumulative =
      mean_error_hybrid - mean_error_cumulative
  )]
  paired_feature_units <- merge(
    paired_feature_units,
    environment_targets[, .(
      environment_id, Location, SeasonYear, milestone, reference_treatment,
      reference_sowing_date, observed_date, observed_environment_days_after_sowing,
      target_definition, target_input_md5, target_source
    )],
    by = comparison_keys, all.x = TRUE, sort = FALSE
  )
  write_table(
    paired_feature_units,
    "briwecs_hybrid_vs_cumulative_paired_units.csv"
  )

  paired_feature_bootstrap <- function(data, replicates, seed) {
    clusters <- unique(data$environment_id)
    point <- c(
      delta_macro_mae = mean(
        data$delta_mae_hybrid_minus_cumulative
      ),
      delta_macro_rmse =
        sqrt(mean(data$mean_squared_error_hybrid)) -
        sqrt(mean(data$mean_squared_error_cumulative)),
      delta_macro_bias = mean(
        data$delta_bias_hybrid_minus_cumulative
      )
    )
    if (length(clusters) < 2L) {
      return(data.table(
        metric = names(point),
        estimate = unname(point),
        ci_low = NA_real_,
        ci_high = NA_real_,
        n_clusters = length(clusters)
      ))
    }
    set.seed(seed)
    draws <- matrix(NA_real_, nrow = replicates, ncol = 3L)
    for (b in seq_len(replicates)) {
      sampled_clusters <- sample(
        clusters,
        length(clusters),
        replace = TRUE
      )
      sampled <- rbindlist(lapply(sampled_clusters, function(cluster) {
        data[environment_id == cluster]
      }))
      draws[b, ] <- c(
        mean(sampled$delta_mae_hybrid_minus_cumulative),
        sqrt(mean(sampled$mean_squared_error_hybrid)) -
          sqrt(mean(sampled$mean_squared_error_cumulative)),
        mean(sampled$delta_bias_hybrid_minus_cumulative)
      )
    }
    data.table(
      metric = names(point),
      estimate = unname(point),
      ci_low = apply(draws, 2L, quantile, 0.025, names = FALSE),
      ci_high = apply(draws, 2L, quantile, 0.975, names = FALSE),
      n_clusters = length(clusters)
    )
  }

  comparison_groups <- c("Overall", sort(unique(paired_feature_units$milestone)))
  paired_feature_summary <- rbindlist(lapply(
    seq_along(comparison_groups),
    function(i) {
      stage <- comparison_groups[i]
      subset <- if (stage == "Overall") {
        paired_feature_units
      } else {
        paired_feature_units[milestone == stage]
      }
      result <- paired_feature_bootstrap(
        subset,
        BOOTSTRAP_REPLICATES,
        SEED + 50000L + i
      )
      result[, milestone := stage]
      result
    }
  ))
  paired_feature_summary[, `:=`(
    target_definition = TARGET_DEFINITION,
    target_input_md5 = ..target_input_md5,
    reference_treatment = REFERENCE_TREATMENT
  )]
  setcolorder(
    paired_feature_summary,
    c("milestone", "metric", "estimate", "ci_low", "ci_high", "n_clusters")
  )
  write_table(
    paired_feature_summary,
    "briwecs_hybrid_vs_cumulative_paired_summary.csv"
  )
} else {
  message(
    "Skipping the optional hybrid-versus-cumulative comparison: ",
    cumulative_comparison_reason
  )
}

# ------------------------------------------------------------------------------
# 12. PUBLICATION-READY FIGURES
# ------------------------------------------------------------------------------

if (RUN_GAM_INTERACTION_ABLATION) {
  overall_ablation <- gam_ablation_summary[milestone == "Overall"]
  ablation_plot_data <- rbindlist(list(
    overall_ablation[, .(
      dataset,
      metric = "RMSE",
      estimate = delta_macro_rmse,
      ci_low = delta_macro_rmse_ci_low,
      ci_high = delta_macro_rmse_ci_high
    )],
    overall_ablation[, .(
      dataset,
      metric = "MAE",
      estimate = delta_macro_mae,
      ci_low = delta_macro_mae_ci_low,
      ci_high = delta_macro_mae_ci_high
    )]
  ))
  ablation_plot_data[, dataset := factor(
    dataset,
    levels = c("Internal development CV", "External BRIWECS")
  )]
  ablation_plot <- ggplot(
    ablation_plot_data,
    aes(x = dataset, y = estimate)
  ) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
    geom_errorbar(
      aes(ymin = ci_low, ymax = ci_high),
      width = 0.10,
      linewidth = 0.55,
      color = "#3A3A3A"
    ) +
    geom_point(size = 2.4, color = "#0072B2") +
    facet_wrap(vars(metric), scales = "free_y") +
    labs(
      x = NULL,
      y = "Interaction - additive error (BBCH stages)",
      title = paste("Predictive GAM comparison:", GAM_ABLATION_SPEC$label),
      subtitle = paste0(
        paste(GAM_ABLATION_FEATURES, collapse = " x "),
        "; external reference treatment: ", REFERENCE_TREATMENT
      ),
      caption = paste0(
        "Negative values favour the interaction; intervals use environment-cluster bootstrap.\n",
        "Internal CV uses DE+AT observations. This comparison does not establish causal effects."
      )
    ) +
    theme_pheno() +
    theme(axis.text.x = element_text(angle = 15, hjust = 1))
  save_plot(
    ablation_plot,
    GAM_ABLATION_SPEC$figure_stem,
    width = 7.4,
    height = 4.8
  )
}

learner_labels <- c(
  xgboost = "XGBoost",
  gam_additive = "GAM: additive",
  gam_interaction = paste0("GAM: ", tolower(GAM_ABLATION_SPEC$label)),
  xgboost_monotone_sensitivity = "XGBoost monotone sensitivity",
  gdd_threshold = "Thermal-time baseline"
)
learner_colors <- c(
  xgboost = "#2878B5",
  gam_additive = "#E69F00",
  gam_interaction = "#009E73",
  xgboost_monotone_sensitivity = "#74ADD1",
  gdd_threshold = "#777777"
)
milestone_levels <- c("BBCH10", "BBCH59", "BBCH87")

plot_units <- copy(environment_stage_units)
plot_units[, milestone := factor(milestone, levels = milestone_levels)]
plot_units[, learner_display := factor(
  learner_id,
  levels = names(learner_labels),
  labels = unname(learner_labels)
)]

prediction_plot <- ggplot(
  plot_units,
  aes(x = milestone, y = mean_prediction, fill = learner_id)
) +
  geom_boxplot(
    width = 0.68,
    outlier.alpha = 0.35,
    position = position_dodge(width = 0.78)
  ) +
  geom_point(
    data = data.table(
      milestone = factor(milestone_levels, levels = milestone_levels),
      truth = c(10, 59, 87)
    ),
    aes(x = milestone, y = truth),
    inherit.aes = FALSE,
    shape = 4,
    stroke = 1.2,
    size = 3.2,
    color = "black"
  ) +
  scale_fill_manual(
    values = learner_colors,
    labels = learner_labels,
    breaks = names(learner_labels)
  ) +
  labs(
    x = NULL,
    y = "Predicted BBCH",
    fill = NULL,
    title = "BRIWECS environment-level milestone predictions",
    subtitle = paste0(
      "One prediction per environment at its target date; reference treatment: ",
      REFERENCE_TREATMENT
    ),
    caption = paste(
      "Crosses mark the nominal BBCH stage. Targets use paper-style cultivar",
      "BLUE averages or the recorded equal-cultivar-mean fallback."
    )
  ) +
  theme_pheno()
save_plot(
  prediction_plot,
  "figure_briwecs_predictions_by_milestone",
  width = 8.8,
  height = 5.4
)

error_plot <- ggplot(
  plot_units,
  aes(x = milestone, y = mean_abs_error, fill = learner_id)
) +
  geom_boxplot(
    width = 0.68,
    outlier.alpha = 0.35,
    position = position_dodge(width = 0.78)
  ) +
  scale_fill_manual(
    values = learner_colors,
    labels = learner_labels,
    breaks = names(learner_labels)
  ) +
  labs(
    x = NULL,
    y = "Absolute environment-target error (BBCH units)",
    fill = NULL,
    title = "External prediction error by milestone",
    subtitle = paste("One loss per environment target; reference treatment:", REFERENCE_TREATMENT)
  ) +
  theme_pheno()
save_plot(
  error_plot,
  "figure_briwecs_absolute_error_by_milestone",
  width = 8.8,
  height = 5.4
)

offset_plot_data <- copy(offset_environment_units[
  is.finite(mean_offset_days)
])
offset_plot_data[, milestone := factor(milestone, levels = milestone_levels)]
offset_plot <- ggplot(
  offset_plot_data,
  aes(x = milestone, y = mean_offset_days, fill = learner_id)
) +
  geom_hline(yintercept = 0, linetype = 2, color = "#A33A2B") +
  geom_boxplot(
    width = 0.68,
    outlier.alpha = 0.3,
    position = position_dodge(width = 0.78)
  ) +
  scale_fill_manual(
    values = learner_colors,
    labels = learner_labels,
    breaks = names(learner_labels)
  ) +
  labs(
    x = NULL,
    y = "Predicted - environment-target timing (days)",
    fill = NULL,
    title = "Timing of predicted BRIWECS milestones",
    subtitle = paste0(
      "First crossing minus unrounded target DAS; reference treatment: ",
      REFERENCE_TREATMENT
    ),
    caption = "Non-crossing trajectories are excluded from timing errors; crossing coverage is reported separately."
  ) +
  theme_pheno()
save_plot(
  offset_plot,
  "figure_briwecs_milestone_date_offsets",
  width = 8.8,
  height = 5.4
)

support_plot_data <- copy(validation_feature_support)
feature_label_lookup <- c(
  cumdevfun_phase = "Phase-aware development index",
  cumDTU_beta = "Cumulative beta thermal time",
  sowing_DOY = "Sowing day of year",
  photoperiod_civil = "Civil-twilight photoperiod",
  paw_fc33_0_100_mm = "Plant-available water capacity (0-100 cm)",
  FTSW_0_100_proxy = "Fraction of transpirable soil water (0-100 cm)",
  sumRR_30 = "Precipitation, preceding 30 days",
  sumRad_30 = "Radiation, preceding 30 days",
  sumCDD_30 = "Cold degree-days, preceding 30 days",
  sumHDD_30 = "Heat degree-days, preceding 30 days",
  sumPEvap_30 = "Potential evaporation, preceding 30 days",
  climaticWaterBalance_30 = "Climatic water balance, preceding 30 days",
  meanSoilWater_0_100_30 = "Mean soil water (0-100 cm), preceding 30 days",
  cumRR = "Cumulative precipitation",
  cumRad = "Cumulative radiation",
  cumCDD = "Cumulative cold degree-days",
  cumHDD = "Cumulative heat degree-days"
)
support_plot_data[, feature_display := factor(
  feature,
  levels = rev(FEATURES),
  labels = unname(feature_label_lookup[rev(FEATURES)])
)]
support_plot <- ggplot(
  support_plot_data,
  aes(
    x = proportion_outside_training_q01_q99,
    y = feature_display
  )
) +
  geom_col(fill = "#2878B5", width = 0.72) +
  scale_x_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    x = "Environment targets outside training 1st-99th percentile",
    y = NULL,
    title = "BRIWECS predictor-support audit",
    subtitle = "Predictors evaluated at the paper-style environment target dates"
  ) +
  theme_pheno() +
  theme(legend.position = "none")
save_plot(
  support_plot,
  "figure_briwecs_predictor_support",
  width = 8.8,
  height = 5.2
)

heatmap_data <- environment_stage_units[learner_id == "xgboost"]
heatmap_data[, environment_display := paste(Location, SeasonYear)]
heatmap_data[, milestone := factor(milestone, levels = milestone_levels)]
environment_order <- heatmap_data[
  , .(mean_absolute_error = mean(mean_abs_error)),
  by = environment_display
][order(mean_absolute_error)]
heatmap_data[, environment_display := factor(
  environment_display,
  levels = environment_order$environment_display
)]
heatmap_plot <- ggplot(
  heatmap_data,
  aes(x = milestone, y = environment_display, fill = mean_error)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#3569A8",
    mid = "white",
    high = "#C84A3A",
    midpoint = 0
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Error\n(BBCH)",
    title = "XGBoost error across independent BRIWECS environments",
    subtitle = "Positive values indicate predictions above the nominal BBCH at the environment target date"
  ) +
  theme_pheno() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 7)
  )
save_plot(
  heatmap_plot,
  "figure_briwecs_environment_error_heatmap",
  width = 8.8,
  height = max(5.4, 1.6 + 0.22 * length(environment_order$environment_display))
)

# ------------------------------------------------------------------------------
# 13. PROVENANCE AND COMPLETION CHECKS
# ------------------------------------------------------------------------------

analysis_configuration <- list(
  validation_run_id = VALIDATION_RUN_ID,
  benchmark_run_id = BENCHMARK_RUN_ID,
  benchmark_dir = normalizePath(BENCHMARK_DIR),
  development_input_file = normalizePath(training_input_file),
  development_input_md5 = training_input_md5,
  external_environment_targets_file = normalizePath(BRIWECS_ENVIRONMENT_TARGETS_FILE),
  external_environment_targets_md5 = target_input_md5,
  target_definition = TARGET_DEFINITION,
  reference_treatment = REFERENCE_TREATMENT,
  target_sources = unique(environment_targets$target_source),
  target_construction_status = paste(
    "environment-level target definition adopted after examining BRIWECS plot structure;",
    "no external outcome used to fit or tune prediction models"
  ),
  target_uncertainty_propagated_in_performance_intervals = FALSE,
  raw_plot_outcomes_loaded = FALSE,
  plot_counts_used_as_weights = FALSE,
  reference_sowing_date_rule = "read from target CSV; must be constant within Location x SeasonYear",
  bbch_evaluation_date_rule = "supplied calendar date equals reference_sowing_date + round(target DAS)",
  date_offset_rule = "first-crossing DAS minus unrounded observed_environment_days_after_sowing",
  deprecated_outputs = c("replicate-weighted metrics", "treatment-factor metrics/plots", "genotype metrics"),
  external_location_file = normalizePath(BRIWECS_LOCATION_FILE),
  external_location_md5 = unname(tools::md5sum(BRIWECS_LOCATION_FILE)),
  external_weather_files = normalizePath(weather_files),
  external_weather_md5 = unname(tools::md5sum(weather_files)),
  external_era5_directory = normalizePath(BRIWECS_ERA5_DIR),
  external_era5_required_variables = REQUIRED_ERA5_VARIABLES,
  external_era5_files = era5_file_audit$file,
  external_era5_md5 = era5_file_audit$md5,
  potential_evaporation_rule = paste(
    "potential_evaporation_mm if supplied, otherwise -1000 * native pev;",
    "negative converted values retained without clipping to match training"
  ),
  external_soil_file = normalizePath(BRIWECS_SOIL_FILE),
  external_soil_md5 = unname(tools::md5sum(BRIWECS_SOIL_FILE)),
  model_role = "full-development-cohort deployment model",
  external_outcomes_used_in_tuning = FALSE,
  inferential_status = if (FEATURE_SET_SHORT == "hybrid") {
    paste(
      "post hoc external sensitivity analysis because the cumulative-model",
      "BRIWECS results were inspected first"
    )
  } else {
    "independent external validation of the cumulative representation"
  },
  feature_set_request = FEATURE_SET_REQUEST,
  feature_set = FEATURE_SET,
  feature_set_definitions = FEATURE_SETS,
  all_external_features_calculated = ALL_EXTERNAL_FEATURES,
  configured_features = FEATURES,
  model_feature_order = MODEL_FEATURE_ORDER,
  final_tuning_scheme = FINAL_TUNING_SCHEME,
  final_tuning_split_ids = tuning_specs$split_id,
  final_tuning_evaluations = TUNING_EVALUATIONS,
  tuning_measure = paste0(
    XGBOOST_TUNING_MEASURE_ID,
    " within development-only internal folds"
  ),
  tuned_parameter_values = tuned_parameter_values,
  prediction_primary = "raw XGBoost response at one rounded environment milestone target date; loss calculated afterward",
  predictions_clipped = FALSE,
  monotone_sensitivity = INCLUDE_MONOTONE_SENSITIVITY,
  monotone_rule = "cumulative maximum over predictions available to date",
  monotone_rule_used_for_primary_stage_metrics = FALSE,
  milestone_date_rule = "first date monotone prediction reaches threshold",
  gam_interaction_ablation = list(
    enabled = RUN_GAM_INTERACTION_ABLATION,
    feature_set = FEATURE_SET,
    features = GAM_ABLATION_FEATURES,
    interaction_pair = GAM_ABLATION_SPEC$pair_id,
    reference = paste(sprintf("s(%s)", GAM_ABLATION_FEATURES), collapse = " + "),
    alternative = paste0(
      paste(sprintf("s(%s)", GAM_ABLATION_FEATURES), collapse = " + "),
      " + ti(", paste(GAM_ABLATION_FEATURES, collapse = ", "), ")"
    ),
    model_scope = "two-predictor GAM comparison; not a full-XGBoost ablation",
    cache_filename = GAM_ABLATION_SPEC$cache_filename,
    refit_allowed = ALLOW_GAM_ABLATION_REFIT,
    external_reference_treatment = REFERENCE_TREATMENT,
    internal_reference_treatment = NA_character_,
    main_effect_basis = "cubic regression spline",
    main_effect_k = GAM_MAIN_K,
    interaction_basis = "tensor-product cubic regression spline",
    interaction_marginal_k = GAM_INTERACTION_K,
    fit_method = "mgcv::bam with fREML, discrete=TRUE and select=TRUE",
    internal_assessment = paste(
      "held-out predictions from the same development-only grouped split",
      "specifications used for final tuning"
    ),
    external_assessment = "models frozen on DE+AT before BRIWECS prediction",
    primary_contrast = "alternative minus reference macro RMSE",
    uncertainty = "environment-cluster bootstrap",
    smooth_interpretation_allowed = FALSE,
    reason = paste(
      "concurvity is audited; performance is interpreted, not smooth",
      "surfaces, coefficients or term p-values"
    )
  ),
  rolling_window_rule = paste0(
    "current day plus preceding width-1 complete calendar days; ",
    RECENT_LOOKBACK_DAYS - 1L,
    " pre-sowing days loaded before process states begin"
  ),
  prediction_horizon_rule = paste0(
    "sowing date through sowing + ", MAX_DAYS_AFTER_SOWING,
    " days, independent of observed milestones"
  ),
  primary_weighting = paste(
    "equal weight per available Location x SeasonYear x milestone target;",
    "Overall retains unit weighting, not equal environment weighting when stage coverage differs"
  ),
  bootstrap_cluster = "Location x SeasonYear environment",
  bootstrap_replicates = BOOTSTRAP_REPLICATES,
  thermal_response = list(
    response_function = "SSM-Wheat beta temperature response",
    base = T_BASE,
    optimum = T_OPTIMUM,
    ceiling = T_CEILING
  ),
  cold_degree_day_threshold = COLD_THRESHOLD,
  heat_degree_day_threshold = HEAT_THRESHOLD,
  process_feature_definition = list(
    feature = "cumdevfun_phase",
    civil_twilight_degrees = SUN_ANGLE_BELOW_HORIZON,
    critical_photoperiod_hours = CPP,
    photoperiod_sensitivity = PPSEN,
    vernalization_saturation_days = VDSAT,
    vernalization_sensitivity = VSEN,
    phase_boundary_biological_days = BD_PHOTO_END,
    observed_bbch_used = FALSE
  ),
  soil_paw_definition = paste0(
    "70*(theta_fc33_l1-theta_wp_l1) + ",
    "210*(theta_fc33_l2-theta_wp_l2) + ",
    "720*(theta_fc33_l3-theta_wp_l3)"
  ),
  ftsw_0_100_definition = paste0(
    "clip(((70*swvl1 + 210*swvl2 + 720*swvl3) - ",
    "wilting_point_0_100_mm) / paw_fc33_0_100_mm, 0, 1)"
  ),
  cumulative_comparison_directory = CUMULATIVE_VALIDATION_DIR,
  cumulative_comparison_available = cumulative_comparison_available,
  cumulative_comparison_status = cumulative_comparison_reason,
  package_versions = as.list(vapply(
    required_packages,
    package_version_or_na,
    character(1)
  )),
  model_cache_version = MODEL_CACHE_VERSION,
  model_source_directory = MODEL_SOURCE_DIR,
  model_refit_allowed = ALLOW_MODEL_REFIT,
  gam_refit_allowed = ALLOW_GAM_ABLATION_REFIT,
  xgboost_cache_loaded_from = model_loaded_from,
  gam_cache_loaded_from = gam_loaded_from,
  n_development_rows = nrow(training),
  n_briwecs_environment_milestone_targets = nrow(environment_targets),
  n_briwecs_environments = uniqueN(milestones[, .(Location, SeasonYear)]),
  n_briwecs_trajectories = nrow(trajectories),
  n_briwecs_daily_rows = nrow(daily),
  n_external_milestone_targets = nrow(milestones)
)
atomic_save_rds(
  analysis_configuration,
  file.path(OUTPUT_DIR, "briwecs_validation_configuration.rds")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "session_info.txt")
)

required_outputs <- c(
  "briwecs_environment_targets_used.csv",
  "briwecs_environment_target_input_audit.csv",
  "briwecs_environment_trajectories.csv",
  "briwecs_environment_milestone_units.csv",
  "briwecs_target_method_support.csv",
  "briwecs_performance_by_target_method.csv",
  "briwecs_hybrid_vs_cumulative_comparison_status.csv",
  "briwecs_validation_configuration.rds",
  "briwecs_primary_macro_performance.csv",
  "briwecs_cluster_bootstrap_performance.csv",
  "briwecs_external_predictions_long.csv",
  "briwecs_era5_file_audit.csv",
  "briwecs_era5_join_audit.csv",
  "briwecs_era5_pev_value_audit.csv",
  "briwecs_ftsw_0_100_audit.csv",
  "briwecs_all_feature_audit.csv",
  "briwecs_milestone_date_offset_summary.csv",
  "briwecs_environment_milestone_date_offset_units.csv",
  "briwecs_held_out_feature_support_audit.csv",
  file.path("model", "final_xgboost_model.rds"),
  file.path("plots", "figure_briwecs_absolute_error_by_milestone.pdf")
)
if (cumulative_comparison_available) {
  required_outputs <- c(
    required_outputs,
    "briwecs_hybrid_vs_cumulative_paired_units.csv",
    "briwecs_hybrid_vs_cumulative_paired_summary.csv"
  )
}
if (RUN_GAM_INTERACTION_ABLATION) {
  required_outputs <- c(
    required_outputs,
    "gam_interaction_ablation_internal_environment_stage_units.csv",
    "interaction_ablation_paired_environment_stage_units.csv",
    "interaction_ablation_performance_summary.csv",
    file.path(
      "model",
      GAM_ABLATION_SPEC$cache_filename
    ),
    file.path(
      "model",
      "gam_interaction_ablation_concurvity.csv"
    ),
    file.path(
      "plots",
      paste0(GAM_ABLATION_SPEC$figure_stem, ".pdf")
    )
  )
}
missing_outputs <- required_outputs[
  !file.exists(file.path(OUTPUT_DIR, required_outputs))
]
if (length(missing_outputs)) {
  stop(
    "External validation finished incompletely. Missing: ",
    paste(missing_outputs, collapse = ", ")
  )
}

message("Environment-target BRIWECS validation finished.")
message("Targets: ", nrow(environment_targets), " across ", nrow(trajectories),
        " environments; reference treatment: ", REFERENCE_TREATMENT)
message("Results written to: ", OUTPUT_DIR)
message(
  "Primary performance: ",
  file.path(OUTPUT_DIR, "briwecs_primary_macro_performance.csv")
)
if (cumulative_comparison_available) {
  message(
    "Important: no BRIWECS outcome entered fitting or tuning"
  )
}

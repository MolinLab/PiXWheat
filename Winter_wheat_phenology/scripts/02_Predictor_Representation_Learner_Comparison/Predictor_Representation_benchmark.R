# ==============================================================================
# Predictor representation and developmental-clock benchmark
# ==============================================================================
# Compare five developmental clocks in a GAM and five predictor representations
# using fixed GAM, random forest, MARS and XGBoost learners. Spatial and
# contiguous temporal folds share a complete-case cohort within this script.
# The primary endpoint is environment-stage macro-MAE.
# Run from the repository root or set PHENOLOGY_PROJECT_DIR to its path.

# ------------------------------------------------------------------------------
# 1. USER CONFIGURATION
# ------------------------------------------------------------------------------
PROJECT_DIR <- Sys.getenv("PHENOLOGY_PROJECT_DIR", unset = getwd())
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = TRUE)
INPUT_FILE <- Sys.getenv(
  "PHENOLOGY_POINT4_INPUT",
  unset = file.path(PROJECT_DIR, "Data", "pheno_df_CV.csv")
)

RUN_ID <- Sys.getenv(
  "PHENOLOGY_POINT4_RUN_ID",
  unset = "point4_hypothesis_benchmark_paper"
)
OUTPUT_DIR <- file.path(PROJECT_DIR, "benchmark_results", RUN_ID)
CHECKPOINT_DIR <- file.path(OUTPUT_DIR, "checkpoints")

SEED <- 123L
TARGET <- "BBCH"
MIN_BBCH <- 0
MAX_BBCH <- 90

CV_SCHEMES <- c("spatial", "temporal_block")
SPATIAL_FOLDS <- 5L
SPATIAL_BLOCK_SIZE_M <- 119946.2
SPATIAL_SELECTION <- "random"
TEMPORAL_FOLDS <- 5L

GAM_ENGINE <- "bam"       # "bam" or "gam"
GAM_FAMILY <- "gaussian"      # "scat" or "gaussian"
GAM_DISCRETE <- TRUE

RF_NUM_TREES <- 500L
XGB_NROUNDS <- 800L

BOOTSTRAP_REPLICATES <- 1000L
WORKERS_REQUESTED <- 10L
RESUME_FROM_CHECKPOINTS <- TRUE
CHECKPOINT_VERSION <- 1L


# ------------------------------------------------------------------------------
# 2. PACKAGES AND REPRODUCIBILITY
# ------------------------------------------------------------------------------

required_packages <- c(
  "data.table", "R6", "mlr3", "mlr3learners", "mlr3extralearners",
  "future", "future.apply", "parallelly", "paradox", "mgcv", "earth", "ranger",
  "xgboost", "mlr3spatiotempcv", "sf", "blockCV"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install these packages before running the representation benchmark: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(R6)
  library(mlr3)
  library(mlr3learners)
  library(mlr3extralearners)
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

# Outer-fold jobs are parallel; each fitted learner stays single-threaded.
data.table::setDTthreads(1L)


# ------------------------------------------------------------------------------
# 3. PRE-SPECIFIED FEATURE SETS
# ------------------------------------------------------------------------------
# All features are available on the observation date. Full-season summaries and
# observed-BBCH-derived phase indicators are excluded.

feature_sets <- list(
  # Mechanism sets each contain one developmental clock.
  elapsed_time = c("days_since_sowing"),
  gdd_linear = c("cumGDD"),
  thermal_beta = c("cumDTU_beta"),
  process_phase = c("cumdevfun_phase"),
  process_phase_water = c("cumdevfun_phase_water_0_100"),
  
  # Predictive representation families.
  cumulative_weather = c(
    "cumDTU_beta", "sowing_DOY", "photoperiod_civil", "cumRR", "cumRad",
    "cumCDD", "cumHDD", "paw_fc33_0_100_mm"
  ),

  calendar_recent_weather = c(
    "days_since_sowing", "sowing_DOY", "TMP", "photoperiod_civil",
    "sumRR_30", "sumRad_30", "sumCDD_30", "sumHDD_30",
    "sumPEvap_30", "climaticWaterBalance_30",
    "meanSoilWater_0_100_30", 
    "paw_fc33_0_100_mm"
  ),

  process_components = c(
    "cumDTU_beta", "sowing_DOY", "meanppfun_phase_to_date",
    "meanvernfun_phase_to_date", "FTSW_0_100_proxy",
    "paw_fc33_0_100_mm"
  ),

  process_context = c(
    "cumdevfun_phase", "FTSW_0_100_proxy","sowing_DOY", "paw_fc33_0_100_mm"
  ),
 
  hybrid_process_recent = c(
    "cumdevfun_phase", "sowing_DOY", "paw_fc33_0_100_mm", "FTSW_0_100_proxy",
    "sumRR_30", "sumRad_30", "sumCDD_30", "sumHDD_30",
    "sumPEvap_30", "climaticWaterBalance_30",
    "meanSoilWater_0_100_30"
  )
)

feature_set_family <- c(
  elapsed_time = "mechanism",
  gdd_linear = "mechanism",
  thermal_beta = "mechanism",
  process_phase = "mechanism",
  process_phase_water = "mechanism",
  cumulative_weather = "representation",
  calendar_recent_weather = "representation",
  process_components = "representation",
  process_context = "representation",
  hybrid_process_recent = "representation"
)

gam_formulas <- list(
  elapsed_time =
    BBCH ~ s(days_since_sowing, k = 50, bs = "ts"),

  gdd_linear =
    BBCH ~ s(cumGDD, k = 50, bs = "ts"),

  thermal_beta =
    BBCH ~ s(cumDTU_beta, k = 50, bs = "ts"),

  process_phase =
    BBCH ~ s(cumdevfun_phase, k = 50, bs = "ts"),

  process_phase_water =
    BBCH ~ s(cumdevfun_phase_water_0_100, k = 50, bs = "ts"),

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

  calendar_recent_weather =
    BBCH ~
      s(days_since_sowing, k = 50, bs = "ts") +
      s(sowing_DOY, k = 12, bs = "ts") +
      s(TMP, k = 12, bs = "ts") +
      s(photoperiod_civil, k = 12, bs = "ts") +
      s(sumRR_30, k = 12, bs = "ts") +
      s(sumRad_30, k = 12, bs = "ts") +
      s(sumCDD_30, k = 12, bs = "ts") +
      s(sumHDD_30, k = 12, bs = "ts") +
      s(sumPEvap_30, k = 12, bs = "ts") +
      s(climaticWaterBalance_30, k = 12, bs = "ts") +
      s(meanSoilWater_0_100_30, k = 12, bs = "ts") +
      s(paw_fc33_0_100_mm, k = 12, bs = "ts") +
      ti(
        days_since_sowing, paw_fc33_0_100_mm,
        k = c(12, 5), bs = c("ts", "ts")
      )+
      ti(
      days_since_sowing, sowing_DOY,
      k = c(12, 5), bs = c("ts", "ts")
     ),

  process_components =
    BBCH ~
      s(cumDTU_beta, k = 50, bs = "ts") +
      s(sowing_DOY, k = 12, bs = "ts") +
      s(meanppfun_phase_to_date, k = 12, bs = "ts") +
      s(meanvernfun_phase_to_date, k = 12, bs = "ts") +
      s(FTSW_0_100_proxy, k = 12, bs = "ts") +
      s(paw_fc33_0_100_mm, k = 12, bs = "ts") +
      ti(
        cumDTU_beta, meanppfun_phase_to_date,
        k = c(12, 5), bs = c("ts", "ts")
      ) +
      ti(
        cumDTU_beta, meanvernfun_phase_to_date,
        k = c(12, 5), bs = c("ts", "ts")
      ) +
      ti(
        cumDTU_beta, FTSW_0_100_proxy,
        k = c(12, 5), bs = c("ts", "ts")
      ),

  process_context =
    BBCH ~
    s(cumdevfun_phase, k = 50, bs = "ts") +
    s(FTSW_0_100_proxy, k=12, bs="ts")+
    s(sowing_DOY, k = 12, bs = "ts") +
    s(paw_fc33_0_100_mm, k = 12, bs = "ts") +
    ti(
      cumdevfun_phase,
      sowing_DOY,
      k = c(12, 5),
      bs = c("ts", "ts")
    ),

  hybrid_process_recent =
    BBCH ~
      s(cumdevfun_phase, k = 50, bs = "ts") +
      s(sowing_DOY, k = 12, bs = "ts") +
      s(paw_fc33_0_100_mm, k = 12, bs = "ts") +
      s(sumRR_30, k = 12, bs = "ts") +
      s(sumRad_30, k = 12, bs = "ts") +
      s(sumCDD_30, k = 12, bs = "ts") +
      s(sumHDD_30, k = 12, bs = "ts") +
      s(sumPEvap_30, k = 12, bs = "ts") +
      s(climaticWaterBalance_30, k = 12, bs = "ts") +
      s(meanSoilWater_0_100_30, k = 12, bs = "ts")+
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

stopifnot(
  identical(names(feature_sets), names(feature_set_family)),
  identical(names(feature_sets), names(gam_formulas))
)


# ------------------------------------------------------------------------------
# 4. PRE-SPECIFIED HYPOTHESES
# ------------------------------------------------------------------------------
# A negative delta means the alternative has lower error than its reference.

hypotheses <- data.table(
  hypothesis_id = paste0("H", 1:8),
  analysis_family = c(
    rep("mechanism", 4),
    rep("representation", 4)
  ),
  reference = c(
    "elapsed_time",
    "gdd_linear",
    "thermal_beta",
    "process_phase",
    "process_components",
    "process_context",
    "cumulative_weather",
    "cumulative_weather"
  ),
  alternative = c(
    "gdd_linear",
    "thermal_beta",
    "process_phase",
    "process_phase_water",
    "process_context",
    "hybrid_process_recent",
    "calendar_recent_weather",
    "hybrid_process_recent"
  ),
  hypothesis = c(
    "Thermal accumulation improves over elapsed calendar time.",
    "The beta thermal response improves over conventional linear GDD.",
    "Photoperiod and vernalization improve over temperature alone.",
    "The water modifier improves phase-aware development.",
    "The integrated process index improves over separately supplied components.",
    "Recent weather adds transferable information to the process index.",
    "Recent-window weather transfers better than cumulative weather.",
    "The process-weather hybrid improves over cumulative weather."
  )
)


# ------------------------------------------------------------------------------
# 5. UTILITIES
# ------------------------------------------------------------------------------

safe_id <- function(x) gsub("[^A-Za-z0-9_-]+", "_", as.character(x))

resolve_column <- function(data, candidates, label) {
  hit <- candidates[candidates %in% names(data)]
  if (!length(hit)) {
    stop(
      "Could not resolve ", label, ". Expected one of: ",
      paste(candidates, collapse = ", ")
    )
  }
  hit[1L]
}

coalesce_character_columns <- function(data, candidates) {
  available <- candidates[candidates %in% names(data)]
  if (!length(available)) return(rep(NA_character_, nrow(data)))

  out <- rep(NA_character_, nrow(data))
  for (column in available) {
    value <- trimws(as.character(data[[column]]))
    value[!nzchar(value)] <- NA_character_
    replace <- is.na(out) & !is.na(value)
    out[replace] <- value[replace]
  }
  out
}

make_coordinate_site_id <- function(longitude, latitude, digits = 5L) {
  paste0(
    formatC(latitude, format = "f", digits = digits),
    "__",
    formatC(longitude, format = "f", digits = digits)
  )
}

make_contiguous_year_folds <- function(year, k) {
  year <- as.integer(year)
  available_years <- sort(unique(year))
  if (length(available_years) < k) {
    stop(
      "Temporal blocking requested ", k, " folds but only ",
      length(available_years), " season years are available."
    )
  }

  # Complete years stay together, and neighbouring years form contiguous blocks.
  year_folds <- pmin(
    k,
    ceiling(seq_along(available_years) * k / length(available_years))
  )
  lookup <- setNames(year_folds, available_years)
  unname(lookup[as.character(year)])
}

prediction_to_table <- function(prediction) {
  row_ids <- prediction$row_ids
  truth <- prediction$truth
  response <- prediction$response

  if (is.null(row_ids) || is.null(truth) || is.null(response)) {
    stop(
      "Prediction lacks row_ids, truth or response. Class: ",
      paste(class(prediction), collapse = "/")
    )
  }
  if (is.data.frame(response)) response <- as.matrix(response)
  if (is.matrix(response) && ncol(response) != 1L) {
    stop("Expected one response column; received ", ncol(response), ".")
  }

  response <- as.numeric(response)
  truth <- as.numeric(truth)

  if (length(row_ids) != length(truth) || length(truth) != length(response)) {
    stop("Prediction fields have inconsistent lengths.")
  }
  if (anyNA(row_ids) || any(!is.finite(truth)) || any(!is.finite(response))) {
    stop("Prediction contains missing IDs or non-finite values.")
  }

  data.table(row_ids = row_ids, truth = truth, response = response)
}

write_table <- function(x, filename) {
  fwrite(x, file.path(OUTPUT_DIR, filename), na = "NA")
}

first_value <- function(x) x[1L]


# ------------------------------------------------------------------------------
# 6. LOAD DATA AND DEFINE ONE COMMON COMPARISON COHORT
# ------------------------------------------------------------------------------

if (!file.exists(INPUT_FILE)) {
  stop(
    "Representation benchmark input not found: ", INPUT_FILE, "\n",
    "Set PHENOLOGY_POINT4_INPUT to the combined German-Austrian modelling table."
  )
}

message("Reading: ", INPUT_FILE)
raw_data <- fread(
  INPUT_FILE,
  na.strings = c("", "NA", "NaN", "Inf", "-Inf"),
  showProgress = TRUE
)
raw_data[, .source_row_id := .I]

longitude_column <- resolve_column(
  raw_data, c("lon", "lon_orig", "longitude"), "longitude"
)
latitude_column <- resolve_column(
  raw_data, c("lat", "lat_orig", "latitude"), "latitude"
)
year_column <- resolve_column(
  raw_data, c("season_year", "SeasonYear", "year"), "season year"
)
date_column <- resolve_column(
  raw_data, c("Date", "phen_date", "season_date"), "observation date"
)
sowing_column <- resolve_column(
  raw_data, c("sowing_DOY", "sow_doy"), "sowing day of year"
)

# Harmonise the two source names for sowing day of year.
raw_data[, sowing_DOY := as.numeric(get(sowing_column))]

required_features <- unique(unlist(feature_sets, use.names = FALSE))
missing_model_columns <- setdiff(
  c(TARGET, required_features),
  names(raw_data)
)
if (length(missing_model_columns)) {
  stop(
    "Required benchmark columns are absent: ",
    paste(missing_model_columns, collapse = ", ")
  )
}

raw_data[, .longitude := as.numeric(get(longitude_column))]
raw_data[, .latitude := as.numeric(get(latitude_column))]
raw_data[, .season_year := as.integer(get(year_column))]
raw_data[, .date := as.IDate(get(date_column))]
raw_data[, .site_id := make_coordinate_site_id(.longitude, .latitude)]

raw_environment_id <- coalesce_character_columns(
  raw_data,
  c("env_id", "traj_id", "environment_id")
)
fallback_environment_id <- paste0(
  raw_data$.site_id, "__", raw_data$.season_year
)
raw_environment_id[is.na(raw_environment_id)] <-
  fallback_environment_id[is.na(raw_environment_id)]

# Prefixing with site and year prevents accidental ID collisions between sources.
raw_data[, .env_id := paste(.site_id, .season_year, raw_environment_id, sep = "__")]

raw_data <- raw_data[
  is.finite(get(TARGET)) &
    get(TARGET) >= MIN_BBCH &
    get(TARGET) <= MAX_BBCH
]

non_numeric_features <- required_features[
  !vapply(raw_data[, ..required_features], is.numeric, logical(1))
]
if (length(non_numeric_features)) {
  stop(
    "All benchmark predictors must be numeric. Non-numeric predictors: ",
    paste(non_numeric_features, collapse = ", ")
  )
}

required_common_columns <- unique(c(
  TARGET, required_features, ".longitude", ".latitude", ".env_id", ".site_id",
  ".season_year", ".date"
))

missingness_audit <- rbindlist(lapply(required_common_columns, function(column) {
  value <- raw_data[[column]]
  data.table(
    column = column,
    n_missing = sum(is.na(value)),
    proportion_missing = mean(is.na(value))
  )
}))
write_table(missingness_audit, "missingness_required_columns.csv")

complete_mask <- complete.cases(raw_data[, ..required_common_columns])
finite_columns <- unique(c(
  TARGET, required_features, ".longitude", ".latitude", ".season_year"
))
for (column in finite_columns) {
  complete_mask <- complete_mask & is.finite(raw_data[[column]])
}

cohort <- copy(raw_data[complete_mask])
cohort[, .obs_id := .I]

cohort_flow <- data.table(
  step = c(
    "input_after_BBCH_range",
    "complete_common_cohort",
    "removed_for_missing_or_nonfinite"
  ),
  n_rows = c(
    nrow(raw_data),
    nrow(cohort),
    nrow(raw_data) - nrow(cohort)
  ),
  proportion_of_filtered_input = c(
    1,
    nrow(cohort) / nrow(raw_data),
    (nrow(raw_data) - nrow(cohort)) / nrow(raw_data)
  )
)
write_table(cohort_flow, "cohort_flow.csv")

if (nrow(cohort) < 100L) {
  stop("The common benchmark cohort contains fewer than 100 observations.")
}


# ------------------------------------------------------------------------------
# 7. QUALITY CONTROL BEFORE FITTING
# ------------------------------------------------------------------------------

if (any(cohort$sowing_DOY < 1 | cohort$sowing_DOY > 366)) {
  stop("sowing_DOY contains values outside 1-366.")
}
if (any(cohort$days_since_sowing < 0)) {
  stop("days_since_sowing contains negative values.")
}
if (any(cohort$photoperiod_civil <= 0 | cohort$photoperiod_civil > 24)) {
  stop("photoperiod_civil contains values outside (0, 24] hours.")
}
if (any(cohort$FTSW_0_100_proxy < -1e-6 | cohort$FTSW_0_100_proxy > 1 + 1e-6)) {
  stop("FTSW_0_100_proxy is outside [0, 1]. Fix the upstream units.")
}

bounded_process_columns <- c(
  "meanppfun_phase_to_date",
  "meanvernfun_phase_to_date"
)
for (column in bounded_process_columns) {
  if (any(cohort[[column]] < -1e-6 | cohort[[column]] > 1 + 1e-6)) {
    stop(column, " is outside [0, 1].")
  }
}

mapping_qc <- cohort[, .(
  n_sites = uniqueN(.site_id),
  n_years = uniqueN(.season_year)
), by = .env_id]

invalid_environment_ids <- mapping_qc[n_sites != 1L | n_years != 1L]
if (nrow(invalid_environment_ids)) {
  write_table(invalid_environment_ids, "invalid_environment_ids.csv")
  stop(
    "At least one env_id maps to multiple physical sites or season years. ",
    "See invalid_environment_ids.csv."
  )
}

monotonic_columns <- c(
  "days_since_sowing", "cumGDD", "cumDTU_beta",
  "cumdevfun_phase", "cumdevfun_phase_water_0_100"
)
monotonic_qc <- copy(cohort)[order(.env_id, .date), {
  decreases <- vapply(
    monotonic_columns,
    function(column) sum(diff(get(column)) < -1e-8),
    integer(1)
  )
  as.list(setNames(decreases, paste0("n_decreases_", monotonic_columns)))
}, by = .env_id]

decrease_columns <- setdiff(names(monotonic_qc), ".env_id")
invalid_monotonic_rows <- rowSums(
  as.matrix(monotonic_qc[, ..decrease_columns])
) > 0L
invalid_monotonic <- monotonic_qc[invalid_monotonic_rows]
if (nrow(invalid_monotonic)) {
  write_table(invalid_monotonic, "non_monotonic_input_trajectories.csv")
  stop("Non-monotonic development-clock trajectories were detected.")
}

duplicate_observations <- cohort[
  ,
  .N,
  by = .(.env_id, .date, observed_BBCH = get(TARGET))
][N > 1L]
write_table(duplicate_observations, "duplicate_environment_date_stage_rows.csv")
if (nrow(duplicate_observations)) {
  stop(
    "Duplicate environment-date-stage observations would change model weights. ",
    "Resolve them upstream; see duplicate_environment_date_stage_rows.csv."
  )
}

feature_qc <- rbindlist(lapply(required_features, function(column) {
  x <- cohort[[column]]
  data.table(
    feature = column,
    n_unique = uniqueN(x),
    minimum = min(x),
    median = median(x),
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

feature_set_table <- rbindlist(lapply(names(feature_sets), function(feature_set) {
  data.table(
    feature_set = feature_set,
    analysis_family = unname(feature_set_family[[feature_set]]),
    feature = feature_sets[[feature_set]],
    n_features = length(feature_sets[[feature_set]])
  )
}))
write_table(feature_set_table, "feature_set_definitions.csv")
write_table(hypotheses, "pre_specified_hypotheses.csv")


# ------------------------------------------------------------------------------
# 8. MLR3 TASKS
# ------------------------------------------------------------------------------

make_task <- function(feature_set, features) {
  backend_columns <- unique(c(
    ".obs_id", ".env_id", ".site_id", TARGET,
    ".longitude", ".latitude", features
  ))
  backend <- cohort[, ..backend_columns]

  task <- TaskRegrST$new(
    id = feature_set,
    backend = backend,
    target = TARGET,
    coordinate_names = c(".longitude", ".latitude"),
    crs = "EPSG:4326"
  )

  # This role is not a predictor. It protects later resampling extensions from
  # ever splitting observations from the same physical site. 
  task$set_col_roles(".site_id", roles = "group")
  task$select(features)
  task
}

tasks <- Map(make_task, names(feature_sets), feature_sets)
names(tasks) <- names(feature_sets)

reference_row_ids <- tasks[[1L]]$row_ids
if (!identical(as.integer(reference_row_ids), cohort$.obs_id)) {
  stop("Task row IDs do not match the common-cohort observation IDs.")
}
if (!all(vapply(
  tasks,
  function(task) identical(task$row_ids, reference_row_ids),
  logical(1)
))) {
  stop("Feature-set tasks do not have identical row IDs.")
}


# ------------------------------------------------------------------------------
# 9. REUSABLE OUTER-FOLD MANIFEST
# ------------------------------------------------------------------------------

fold_vectors <- list()

if ("spatial" %in% CV_SCHEMES) {
  # Assign one spatial fold per unique coordinate-defined physical site.
  site_table <- cohort[, .(
    .longitude = first_value(.longitude),
    .latitude = first_value(.latitude),
    spatial_target = mean(get(TARGET))
  ), by = .site_id]
  site_table[, .site_row_id := .I]

  site_backend <- site_table[, .(
    spatial_target,
    .longitude,
    .latitude
  )]

  site_task <- TaskRegrST$new(
    id = "physical_sites_for_spatial_folds",
    backend = site_backend,
    target = "spatial_target",
    coordinate_names = c(".longitude", ".latitude"),
    crs = "EPSG:4326"
  )

  set.seed(SEED + 20L)
  spatial_resampling <- rsmp(
    "spcv_block",
    range = as.integer(SPATIAL_BLOCK_SIZE_M),
    folds = as.integer(SPATIAL_FOLDS),
    selection = SPATIAL_SELECTION,
    hexagon = FALSE
  )
  spatial_resampling$instantiate(site_task)

  spatial_instance <- as.data.table(spatial_resampling$instance)
  row_id_column <- intersect(
    c("row_id", "row_ids"),
    names(spatial_instance)
  )[1L]
  if (!length(row_id_column) || is.na(row_id_column)) {
    stop("Could not identify row IDs in the spatial resampling instance.")
  }

  site_folds <- spatial_instance$fold[
    match(site_task$row_ids, spatial_instance[[row_id_column]])
  ]
  if (anyNA(site_folds)) {
    stop("Spatial CV did not assign every physical site to a fold.")
  }

  site_fold_lookup <- setNames(site_folds, site_table$.site_id)
  fold_vectors$spatial <- unname(site_fold_lookup[cohort$.site_id])
}

if ("temporal_block" %in% CV_SCHEMES) {
  fold_vectors$temporal_block <- make_contiguous_year_folds(
    cohort$.season_year,
    TEMPORAL_FOLDS
  )
}

if (!identical(sort(names(fold_vectors)), sort(CV_SCHEMES))) {
  stop("One or more requested CV schemes were not constructed.")
}

for (scheme in names(fold_vectors)) {
  env_check <- data.table(
    env_id = cohort$.env_id,
    fold = as.character(fold_vectors[[scheme]])
  )[, .(n_folds = uniqueN(fold)), by = env_id]
  if (any(env_check$n_folds != 1L)) {
    stop("env_id was split across folds in scheme: ", scheme)
  }
}

spatial_site_check <- data.table(
  site_id = cohort$.site_id,
  fold = as.character(fold_vectors$spatial)
)[, .(n_folds = uniqueN(fold)), by = site_id]
if (any(spatial_site_check$n_folds != 1L)) {
  stop("A physical site was split across spatial folds.")
}

temporal_year_check <- data.table(
  year = cohort$.season_year,
  fold = as.character(fold_vectors$temporal_block)
)[, .(n_folds = uniqueN(fold)), by = year]
if (any(temporal_year_check$n_folds != 1L)) {
  stop("A season year was split across temporal blocks.")
}

fold_manifest <- cohort[, .(
  obs_id = .obs_id,
  source_row_id = .source_row_id,
  env_id = .env_id,
  site_id = .site_id,
  season_year = .season_year,
  Date = .date,
  BBCH = get(TARGET),
  lon = .longitude,
  lat = .latitude
)]
for (scheme in names(fold_vectors)) {
  fold_manifest[
    ,
    (paste0("fold_", scheme)) := as.character(fold_vectors[[scheme]])
  ]
}
write_table(fold_manifest, "outer_fold_manifest.csv")

fold_sizes <- rbindlist(lapply(names(fold_vectors), function(scheme) {
  ans <- data.table(
    fold = as.character(fold_vectors[[scheme]])
  )[, .N, by = fold]
  ans[, scheme := scheme]
  ans
}))
setcolorder(fold_sizes, c("scheme", "fold", "N"))
write_table(fold_sizes, "outer_fold_sizes.csv")

temporal_fold_years <- unique(
  fold_manifest[, .(scheme = "temporal_block", fold_temporal_block, season_year)]
)
setnames(temporal_fold_years, "fold_temporal_block", "fold")
setorder(temporal_fold_years, fold, season_year)
write_table(temporal_fold_years, "temporal_block_years.csv")


# ------------------------------------------------------------------------------
# 10. TRAINING-ONLY CALIBRATED GDD-THRESHOLD BASELINE
# ------------------------------------------------------------------------------
# create learner that calculates the median for each stage and uses isotonic regression to force non-decreasing GDD with increasing BBCH
LearnerRegrGDDThreshold <- R6Class(
  "LearnerRegrGDDThreshold",
  inherit = LearnerRegr,
  public = list(
    gdd_column = NULL,

    initialize = function(gdd_column = "cumGDD", id = "gdd_threshold") {
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
      list(response = response)
    }
  )
)


# ------------------------------------------------------------------------------
# 11. CUSTOM MGCV LEARNER
# ------------------------------------------------------------------------------

LearnerRegrGamMgcv <- R6Class(
  "LearnerRegrGamMgcv",
  inherit = LearnerRegr,
  public = list(
    initialize = function(formula, family, engine = "bam", id = "gam") {
      parameter_set <- paradox::ps(
        formula = paradox::p_uty(),
        family = paradox::p_uty(),
        engine = paradox::p_fct(c("gam", "bam")),
        discrete = paradox::p_lgl()
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
        family = family,
        engine = engine,
        discrete = GAM_DISCRETE
      )
    }
  ),
  private = list(
    .train = function(task) {
      data <- task$data()
      values <- self$param_set$values

      if (identical(values$engine, "bam")) {
        mgcv::bam(
          formula = values$formula,
          data = data,
          family = values$family,
          method = "fREML",
          select = TRUE,
          discrete = values$discrete,
          nthreads = 1L
        )
      } else {
        mgcv::gam(
          formula = values$formula,
          data = data,
          family = values$family,
          method = "REML",
          select = TRUE
        )
      }
    },

    .predict = function(task) {
      newdata <- task$data(cols = task$feature_names)
      list(response = as.numeric(
        predict(self$model, newdata = newdata, type = "response")
      ))
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
# 12. FIXED SCREENING LEARNERS
# ------------------------------------------------------------------------------

make_learners <- function(feature_set, task) {
  p <- length(task$feature_names)

  learner_gam <- LearnerRegrGamMgcv$new(
    formula = gam_formulas[[feature_set]],
    family = gam_family,
    engine = GAM_ENGINE,
    id = "gam"
  )

  if (identical(feature_set_family[[feature_set]], "mechanism")) {
    learners <- list(gam = learner_gam)
    if (identical(feature_set, "gdd_linear")) {
      learners <- c(
        list(gdd_threshold = LearnerRegrGDDThreshold$new("cumGDD")),
        learners
      )
    }
    return(learners)
  }

  learner_rf <- lrn(
    "regr.ranger",
    num.trees = RF_NUM_TREES,
    mtry = max(1L, floor(sqrt(p))),
    min.node.size = 5L,
    sample.fraction = 0.8,
    importance = "none",
    num.threads = 1L,
    seed = SEED,
    predict_type = "response"
  )

  learner_mars <- lrn(
    "regr.earth",
    degree = min(2L, p),
    penalty = 3,
    nprune = max(8L, min(40L, 3L * p + 5L)),
    predict_type = "response"
  )

  learner_xgb <- lrn(
    "regr.xgboost",
    objective = "reg:squarederror",
    booster = "gbtree",
    tree_method = "hist",
    nrounds = XGB_NROUNDS,
    eta = 0.03,
    max_depth = 6L,
    min_child_weight = 5,
    subsample = 0.8,
    colsample_bytree = 0.8,
    lambda = 1,
    alpha = 0,
    nthread = 1L,
    seed = SEED,
    predict_type = "response"
  )

  list(
    gam = learner_gam,
    random_forest = learner_rf,
    mars = learner_mars,
    xgboost = learner_xgb
  )
}

learners_by_task <- Map(make_learners, names(tasks), tasks)
names(learners_by_task) <- names(tasks)


# ------------------------------------------------------------------------------
# 13. EXPLICIT OUTER-CV JOBS
# ------------------------------------------------------------------------------

jobs <- rbindlist(lapply(names(tasks), function(feature_set) {
  rbindlist(lapply(names(learners_by_task[[feature_set]]), function(learner_id) {
    rbindlist(lapply(names(fold_vectors), function(scheme) {
      fold_ids <- sort(unique(as.character(fold_vectors[[scheme]])))
      data.table(
        feature_set = feature_set,
        analysis_family = unname(feature_set_family[[feature_set]]),
        learner_id = learner_id,
        scheme = scheme,
        fold_id = fold_ids,
        iteration = seq_along(fold_ids)
      )
    }))
  }))
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
    cached_result <- tryCatch(readRDS(checkpoint_file), error = function(e) NULL)
    if (
      !is.null(cached_result) &&
        identical(cached_result$checkpoint_version, CHECKPOINT_VERSION)
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
    stop("Empty train or test set in outer job ", job$job_id, ".")
  }

  if (length(intersect(
    unique(cohort$.env_id[train_ids]),
    unique(cohort$.env_id[test_ids])
  ))) {
    stop("Environment leakage in outer job ", job$job_id, ".")
  }

  if (
    identical(job$scheme, "spatial") &&
      length(intersect(
        unique(cohort$.site_id[train_ids]),
        unique(cohort$.site_id[test_ids])
      ))
  ) {
    stop("Physical-site leakage in spatial outer job ", job$job_id, ".")
  }

  if (
    identical(job$scheme, "temporal_block") &&
      length(intersect(
        unique(cohort$.season_year[train_ids]),
        unique(cohort$.season_year[test_ids])
      ))
  ) {
    stop("Season-year leakage in temporal outer job ", job$job_id, ".")
  }

  warning_messages <- character()
  start_time <- proc.time()[["elapsed"]]

  result <- tryCatch(
    withCallingHandlers({
      learner$train(task, row_ids = train_ids)
      prediction <- learner$predict(task, row_ids = test_ids)
      prediction_table <- prediction_to_table(prediction)
      prediction_table[, feature_set := job$feature_set]
      prediction_table[, analysis_family := job$analysis_family]
      prediction_table[, learner_id := job$learner_id]
      prediction_table[, scheme := job$scheme]
      prediction_table[, fold_id := as.character(job$fold_id)]
      prediction_table[, iteration := as.integer(job$iteration)]
      error <- prediction_table$response - prediction_table$truth

      list(
        status = "ok",
        job = copy(job),
        predictions = prediction_table,
        fold_metrics = data.table(
          feature_set = job$feature_set,
          analysis_family = job$analysis_family,
          learner_id = job$learner_id,
          scheme = job$scheme,
          fold_id = as.character(job$fold_id),
          iteration = as.integer(job$iteration),
          n_test = nrow(prediction_table),
          rmse = sqrt(mean(error^2)),
          mae = mean(abs(error)),
          bias = mean(error),
          within_5_bbch = mean(abs(error) <= 5),
          out_of_bounds = mean(
            prediction_table$response < MIN_BBCH |
              prediction_table$response > MAX_BBCH
          )
        ),
        warnings = warning_messages,
        error = NA_character_
      )
    }, warning = function(w) {
      warning_messages <<- unique(c(warning_messages, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }),
    error = function(e) {
      list(
        status = "error",
        job = copy(job),
        predictions = NULL,
        fold_metrics = data.table(
          feature_set = job$feature_set,
          analysis_family = job$analysis_family,
          learner_id = job$learner_id,
          scheme = job$scheme,
          fold_id = as.character(job$fold_id),
          iteration = as.integer(job$iteration),
          n_test = length(test_ids),
          rmse = NA_real_,
          mae = NA_real_,
          bias = NA_real_,
          within_5_bbch = NA_real_,
          out_of_bounds = NA_real_
        ),
        warnings = warning_messages,
        error = conditionMessage(e)
      )
    }
  )

  result$checkpoint_version <- CHECKPOINT_VERSION
  result$elapsed_seconds <- proc.time()[["elapsed"]] - start_time

  temporary_checkpoint <- paste0(checkpoint_file, ".tmp_", Sys.getpid())
  saveRDS(result, temporary_checkpoint, compress = FALSE)
  if (!file.rename(temporary_checkpoint, checkpoint_file)) {
    stop("Could not atomically write checkpoint: ", checkpoint_file)
  }

  result
}

message(
  "Starting ", nrow(jobs), " outer-fold jobs with ", workers,
  " worker(s)."
)
options(future.globals.maxSize = 2 * 1024^3)

future::plan(future::multisession, workers = workers)
outer_results <- future.apply::future_lapply(
  seq_len(nrow(jobs)),
  function(i) run_outer_job(jobs[i]),
  future.seed = TRUE,
  future.scheduling = 2
)
future::plan(future::sequential)

saveRDS(
  outer_results,
  file.path(OUTPUT_DIR, "outer_results_compact.rds"),
  compress = FALSE
)


# ------------------------------------------------------------------------------
# 14. COMBINE PREDICTIONS AND JOB RESULTS
# ------------------------------------------------------------------------------

job_status <- rbindlist(lapply(outer_results, function(result) {
  data.table(
    job_id = result$job$job_id,
    feature_set = result$job$feature_set,
    analysis_family = result$job$analysis_family,
    learner_id = result$job$learner_id,
    scheme = result$job$scheme,
    fold_id = as.character(result$job$fold_id),
    iteration = result$job$iteration,
    status = result$status,
    elapsed_seconds = result$elapsed_seconds,
    warnings = paste(result$warnings, collapse = " | "),
    error = result$error
  )
}), fill = TRUE)
write_table(job_status, "job_status.csv")

fold_metrics <- rbindlist(
  lapply(outer_results, function(x) x[["fold_metrics"]]),
  use.names = TRUE,
  fill = TRUE
)
write_table(fold_metrics, "fold_metrics.csv")

failed_jobs <- job_status[status != "ok"]
if (nrow(failed_jobs)) {
  write_table(failed_jobs, "failed_jobs.csv")
  stop(
    nrow(failed_jobs),
    " outer-fold jobs failed. Inspect failed_jobs.csv before interpretation."
  )
}

oof_predictions <- rbindlist(
  lapply(outer_results, function(x) x[["predictions"]]),
  use.names = TRUE,
  fill = TRUE
)

metadata_columns <- unique(c(
  ".obs_id", ".source_row_id", ".env_id", ".site_id", ".season_year",
  ".date",
  intersect(
    c(
      "region", "Location", "location", "Stationsname", "Source_name",
      "Experiment_name", "country"
    ),
    names(cohort)
  )
))
prediction_metadata <- cohort[, ..metadata_columns]

oof_predictions <- merge(
  oof_predictions,
  prediction_metadata,
  by.x = "row_ids",
  by.y = ".obs_id",
  all.x = TRUE,
  sort = FALSE
)
setnames(
  oof_predictions,
  c(".source_row_id", ".env_id", ".site_id", ".season_year", ".date"),
  c("source_row_id", "env_id", "site_id", "season_year", "Date")
)
oof_predictions[, error := response - truth]
oof_predictions[, absolute_error := abs(error)]
oof_predictions[, squared_error := error^2]
write_table(oof_predictions, "out_of_fold_predictions.csv")


# ------------------------------------------------------------------------------
# 15. PRIMARY AND SECONDARY PERFORMANCE SUMMARIES
# ------------------------------------------------------------------------------

overall_metrics <- oof_predictions[, .(
  n = .N,
  rmse = sqrt(mean(squared_error)),
  mae = mean(absolute_error),
  bias = mean(error),
  within_5_bbch = mean(absolute_error <= 5),
  out_of_bounds = mean(response < MIN_BBCH | response > MAX_BBCH)
), by = .(
  analysis_family, feature_set, learner_id, scheme
)]

stage_metrics <- oof_predictions[, .(
  n = .N,
  rmse = sqrt(mean(squared_error)),
  mae = mean(absolute_error),
  bias = mean(error),
  within_5_bbch = mean(absolute_error <= 5)
), by = .(
  analysis_family, feature_set, learner_id, scheme, BBCH_stage = truth
)]

environment_stage_metrics <- oof_predictions[, .(
  n = .N,
  mae = mean(absolute_error),
  mse = mean(squared_error),
  bias = mean(error)
), by = .(
  analysis_family, feature_set, learner_id, scheme,
  env_id, site_id, season_year, BBCH_stage = truth
)]
environment_stage_metrics[, rmse := sqrt(mse)]


macro_environment_stage <- environment_stage_metrics[, .(
  n_environment_stage_units = .N,
  macro_mae = mean(mae),
  macro_rmse = sqrt(mean(mse)),
  macro_absolute_bias = mean(abs(bias))
), by = .(
  analysis_family, feature_set, learner_id, scheme
)]

stage_extremes <- stage_metrics[, .(
  n_stages = .N,
  worst_stage_rmse = max(rmse),
  worst_stage_mae = max(mae),
  worst_rmse_stage = BBCH_stage[which.max(rmse)],
  worst_mae_stage = BBCH_stage[which.max(mae)]
), by = .(
  analysis_family, feature_set, learner_id, scheme
)]

predictor_set_summary <- merge(
  overall_metrics,
  macro_environment_stage,
  by = c("analysis_family", "feature_set", "learner_id", "scheme"),
  all = TRUE
)
predictor_set_summary <- merge(
  predictor_set_summary,
  stage_extremes,
  by = c("analysis_family", "feature_set", "learner_id", "scheme"),
  all = TRUE
)

predictor_set_summary[, rank_macro_mae := frank(
  macro_mae,
  ties.method = "average"
), by = .(analysis_family, learner_id, scheme)]
predictor_set_summary[, rank_macro_rmse := frank(
  macro_rmse,
  ties.method = "average"
), by = .(analysis_family, learner_id, scheme)]

setorder(
  predictor_set_summary,
  analysis_family,
  scheme,
  learner_id,
  rank_macro_mae
)

representation_rank_summary <- predictor_set_summary[
  analysis_family == "representation",
  .(
    mean_rank_macro_mae = mean(rank_macro_mae),
    worst_rank_macro_mae = max(rank_macro_mae),
    mean_macro_mae = mean(macro_mae),
    mean_macro_rmse = mean(macro_rmse),
    n_learner_scheme_combinations = .N
  ),
  by = feature_set
][order(mean_rank_macro_mae, mean_macro_mae)]

write_table(overall_metrics, "overall_row_weighted_metrics.csv")
write_table(stage_metrics, "stage_specific_metrics.csv")
write_table(environment_stage_metrics, "environment_stage_metrics.csv")
write_table(predictor_set_summary, "predictor_set_summary.csv")
write_table(
  representation_rank_summary,
  "representation_rank_summary.csv"
)


# ------------------------------------------------------------------------------
# 16. PRE-SPECIFIED PAIRED HYPOTHESIS CONTRASTS
# ------------------------------------------------------------------------------

make_hypothesis_contrast <- function(hypothesis_row) {
  reference_predictions <- oof_predictions[
    feature_set == hypothesis_row$reference,
    .(
      row_ids, learner_id, scheme, fold_id, truth,
      reference_response = response,
      reference_absolute_error = absolute_error,
      reference_squared_error = squared_error,
      env_id, site_id, season_year
    )
  ]

  alternative_predictions <- oof_predictions[
    feature_set == hypothesis_row$alternative,
    .(
      row_ids, learner_id, scheme, fold_id, truth,
      alternative_response = response,
      alternative_absolute_error = absolute_error,
      alternative_squared_error = squared_error
    )
  ]

  paired <- merge(
    reference_predictions,
    alternative_predictions,
    by = c("row_ids", "learner_id", "scheme", "fold_id", "truth"),
    all = FALSE
  )

  paired[, hypothesis_id := hypothesis_row$hypothesis_id]
  paired[, analysis_family := hypothesis_row$analysis_family]
  paired[, hypothesis := hypothesis_row$hypothesis]
  paired[, reference := hypothesis_row$reference]
  paired[, alternative := hypothesis_row$alternative]
  paired[
    ,
    delta_absolute_error :=
      alternative_absolute_error - reference_absolute_error
  ]
  paired[
    ,
    delta_squared_error :=
      alternative_squared_error - reference_squared_error
  ]
  paired
}

paired_contrasts <- rbindlist(
  lapply(seq_len(nrow(hypotheses)), function(i) {
    make_hypothesis_contrast(hypotheses[i])
  }),
  use.names = TRUE,
  fill = TRUE
)

if (!nrow(paired_contrasts)) {
  stop("No paired predictions were available for the hypotheses.")
}

contrast_units <- paired_contrasts[, .(
  n_observations = .N,
  reference_mae = mean(reference_absolute_error),
  alternative_mae = mean(alternative_absolute_error),
  reference_mse = mean(reference_squared_error),
  alternative_mse = mean(alternative_squared_error)
), by = .(
  hypothesis_id, analysis_family, hypothesis, reference, alternative,
  learner_id, scheme, fold_id, env_id, site_id, season_year,
  BBCH_stage = truth
)]

hypothesis_summary <- contrast_units[, .(
  n_environment_stage_units = .N,
  reference_macro_mae = mean(reference_mae),
  alternative_macro_mae = mean(alternative_mae),
  delta_macro_mae = mean(alternative_mae) - mean(reference_mae),
  reference_macro_rmse = sqrt(mean(reference_mse)),
  alternative_macro_rmse = sqrt(mean(alternative_mse)),
  delta_macro_rmse =
    sqrt(mean(alternative_mse)) - sqrt(mean(reference_mse)),
  proportion_units_alternative_better =
    mean(alternative_mae < reference_mae),
  proportion_units_tied =
    mean(abs(alternative_mae - reference_mae) < 1e-12)
), by = .(
  hypothesis_id, analysis_family, hypothesis, reference, alternative,
  learner_id, scheme
)]

hypothesis_stage_summary <- contrast_units[, .(
  n_environment_stage_units = .N,
  reference_macro_mae = mean(reference_mae),
  alternative_macro_mae = mean(alternative_mae),
  delta_macro_mae = mean(alternative_mae) - mean(reference_mae),
  reference_macro_rmse = sqrt(mean(reference_mse)),
  alternative_macro_rmse = sqrt(mean(alternative_mse)),
  delta_macro_rmse =
    sqrt(mean(alternative_mse)) - sqrt(mean(reference_mse))
), by = .(
  hypothesis_id, analysis_family, hypothesis, reference, alternative,
  learner_id, scheme, BBCH_stage
)]

hypothesis_fold_summary <- contrast_units[, .(
  n_environment_stage_units = .N,
  delta_macro_mae = mean(alternative_mae) - mean(reference_mae),
  delta_macro_rmse =
    sqrt(mean(alternative_mse)) - sqrt(mean(reference_mse))
), by = .(
  hypothesis_id, analysis_family, hypothesis, reference, alternative,
  learner_id, scheme, fold_id
)]

bootstrap_one_contrast <- function(data, bootstrap_replicates, seed) {
  if (uniqueN(data$scheme) != 1L) {
    stop("bootstrap_one_contrast requires one CV scheme.")
  }

  scheme_value <- data$scheme[1L]
  data[, bootstrap_cluster := if (
    identical(scheme_value, "spatial")
  ) site_id else as.character(season_year)]

  cluster_statistics <- data[, .(
    n_units = .N,
    sum_reference_mae = sum(reference_mae),
    sum_alternative_mae = sum(alternative_mae),
    sum_reference_mse = sum(reference_mse),
    sum_alternative_mse = sum(alternative_mse)
  ), by = bootstrap_cluster]

  n_clusters <- nrow(cluster_statistics)
  if (n_clusters < 8L || bootstrap_replicates < 1L) {
    return(data.table(
      n_bootstrap_clusters = n_clusters,
      delta_macro_mae_ci_low = NA_real_,
      delta_macro_mae_ci_high = NA_real_,
      delta_macro_rmse_ci_low = NA_real_,
      delta_macro_rmse_ci_high = NA_real_,
      bootstrap_probability_alternative_better_mae = NA_real_
    ))
  }

  set.seed(seed)
  boot_mae <- numeric(bootstrap_replicates)
  boot_rmse <- numeric(bootstrap_replicates)

  for (b in seq_len(bootstrap_replicates)) {
    sampled <- sample.int(n_clusters, n_clusters, replace = TRUE)
    denominator <- sum(cluster_statistics$n_units[sampled])

    reference_mae <- sum(
      cluster_statistics$sum_reference_mae[sampled]
    ) / denominator
    alternative_mae <- sum(
      cluster_statistics$sum_alternative_mae[sampled]
    ) / denominator

    reference_rmse <- sqrt(
      sum(cluster_statistics$sum_reference_mse[sampled]) / denominator
    )
    alternative_rmse <- sqrt(
      sum(cluster_statistics$sum_alternative_mse[sampled]) / denominator
    )

    boot_mae[b] <- alternative_mae - reference_mae
    boot_rmse[b] <- alternative_rmse - reference_rmse
  }

  data.table(
    n_bootstrap_clusters = n_clusters,
    delta_macro_mae_ci_low = unname(quantile(boot_mae, 0.025)),
    delta_macro_mae_ci_high = unname(quantile(boot_mae, 0.975)),
    delta_macro_rmse_ci_low = unname(quantile(boot_rmse, 0.025)),
    delta_macro_rmse_ci_high = unname(quantile(boot_rmse, 0.975)),
    bootstrap_probability_alternative_better_mae = mean(boot_mae < 0)
  )
}

bootstrap_groups <- unique(
  contrast_units[, .(hypothesis_id, learner_id, scheme)]
)
bootstrap_results <- rbindlist(lapply(seq_len(nrow(bootstrap_groups)), function(i) {
  group <- bootstrap_groups[i]
  subset_data <- contrast_units[
    hypothesis_id == group$hypothesis_id &
      learner_id == group$learner_id &
      scheme == group$scheme
  ]

  result <- bootstrap_one_contrast(
    copy(subset_data),
    BOOTSTRAP_REPLICATES,
    SEED + 1000L + i
  )
  result[, hypothesis_id := group$hypothesis_id]
  result[, learner_id := group$learner_id]
  result[, scheme := group$scheme]
  result
}))

hypothesis_summary <- merge(
  hypothesis_summary,
  bootstrap_results,
  by = c("hypothesis_id", "learner_id", "scheme"),
  all.x = TRUE
)
hypothesis_summary[, conclusion := fifelse(
  !is.na(delta_macro_mae_ci_high) & delta_macro_mae_ci_high < 0,
  "supports_alternative",
  fifelse(
    !is.na(delta_macro_mae_ci_low) & delta_macro_mae_ci_low > 0,
    "supports_reference",
    "inconclusive"
  )
)]

setorder(hypothesis_summary, hypothesis_id, scheme, learner_id)

write_table(paired_contrasts, "paired_hypothesis_predictions.csv")
write_table(contrast_units, "paired_hypothesis_environment_stage_units.csv")
write_table(hypothesis_summary, "hypothesis_summary.csv")
write_table(hypothesis_stage_summary, "hypothesis_stage_summary.csv")
write_table(hypothesis_fold_summary, "hypothesis_fold_summary.csv")


# ------------------------------------------------------------------------------
# 17. DEVELOPMENT-TRAJECTORY DIAGNOSTICS
# ------------------------------------------------------------------------------

trajectory_diagnostics <- oof_predictions[order(
  analysis_family, feature_set, learner_id, scheme, env_id, Date
), .(
  n_observations = .N,
  n_prediction_decreases = if (.N > 1L) {
    sum(diff(response) < -1e-8)
  } else {
    0L
  },
  largest_prediction_decrease = if (.N > 1L) {
    min(c(0, diff(response)))
  } else {
    0
  }
), by = .(
  analysis_family, feature_set, learner_id, scheme, env_id
)]

trajectory_summary <- trajectory_diagnostics[, .(
  n_environments = .N,
  environments_with_decrease = sum(n_prediction_decreases > 0L),
  proportion_environments_with_decrease = mean(n_prediction_decreases > 0L),
  total_prediction_decreases = sum(n_prediction_decreases),
  largest_prediction_decrease = min(largest_prediction_decrease)
), by = .(
  analysis_family, feature_set, learner_id, scheme
)]

write_table(
  trajectory_diagnostics,
  "trajectory_diagnostics_by_environment.csv"
)
write_table(trajectory_summary, "trajectory_diagnostics_summary.csv")


# ------------------------------------------------------------------------------
# 18. CONFIGURATION AND PROVENANCE
# ------------------------------------------------------------------------------

configuration <- list(
  run_id = RUN_ID,
  input_file = normalizePath(INPUT_FILE),
  input_md5 = unname(tools::md5sum(INPUT_FILE)),
  seed = SEED,
  target = TARGET,
  BBCH_range = c(MIN_BBCH, MAX_BBCH),
  feature_sets = feature_sets,
  feature_set_family = feature_set_family,
  hypotheses = hypotheses,
  gam_formulas = lapply(
    gam_formulas,
    function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
  ),
  cv_schemes = CV_SCHEMES,
  spatial_folds = SPATIAL_FOLDS,
  spatial_block_size_m = SPATIAL_BLOCK_SIZE_M,
  spatial_selection = SPATIAL_SELECTION,
  temporal_folds = TEMPORAL_FOLDS,
  rf_num_trees = RF_NUM_TREES,
  xgb_nrounds = XGB_NROUNDS,
  gam_engine = GAM_ENGINE,
  gam_family = GAM_FAMILY,
  gam_discrete = GAM_DISCRETE,
  bootstrap_replicates = BOOTSTRAP_REPLICATES,
  checkpoint_version = CHECKPOINT_VERSION,
  workers = workers,
  n_input_rows_after_BBCH_range = nrow(raw_data),
  n_common_cohort_rows = nrow(cohort),
  n_environments = uniqueN(cohort$.env_id),
  n_sites = uniqueN(cohort$.site_id),
  n_years = uniqueN(cohort$.season_year)
)

saveRDS(configuration, file.path(OUTPUT_DIR, "benchmark_configuration.rds"))
writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "session_info.txt")
)

message("Predictor-representation benchmark finished successfully.")
message("Results written to: ", OUTPUT_DIR)
message(
  "Primary hypothesis table: ",
  file.path(OUTPUT_DIR, "hypothesis_summary.csv")
)
message(
  "Representation ranking: ",
  file.path(OUTPUT_DIR, "representation_rank_summary.csv")
)


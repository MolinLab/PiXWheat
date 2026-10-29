# ==============================================================================
# Out-of-fold TreeSHAP for the selected winter-wheat phenology model
# ==============================================================================
# Refit the saved XGBoost model for each requested outer fold, verify its
# predictions against the generalisation benchmark, and explain held-out
# observations only. Results are reported separately by transfer scheme.
# TreeSHAP describes fitted-model attribution, not causal effects.
# Run from the repository root or set PHENOLOGY_PROJECT_DIR to its path.

# ------------------------------------------------------------------------------
# 1. USER CONFIGURATION
# ------------------------------------------------------------------------------
PROJECT_DIR <- Sys.getenv("PHENOLOGY_PROJECT_DIR", unset = getwd())
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = TRUE)

BENCHMARK_RUN_ID <- Sys.getenv(
  "PHENOLOGY_GENERALISATION_RUN_ID",
  unset = "multischeme_xgboost_paper_cum_2"
)

BENCHMARK_DIR <- Sys.getenv(
  "PHENOLOGY_GENERALISATION_BENCHMARK_DIR",
  unset = file.path(
    PROJECT_DIR,
    "benchmark_results",
    BENCHMARK_RUN_ID
  )
)

SHAP_RUN_ID <- Sys.getenv(
  "PHENOLOGY_SHAP_RUN_ID",
  unset = "xgboost_generalisation_oof_treeshap_paper_cum_2"
)

OUTPUT_DIR <- file.path(
  PROJECT_DIR,
  "shap_results",
  SHAP_RUN_ID
)
CHECKPOINT_DIR <- file.path(OUTPUT_DIR, "checkpoints")
PLOT_DIR <- file.path(OUTPUT_DIR, "plots")

# The selected representation is read from benchmark_configuration.rds by
# default. Override only after the representation decision has been frozen.
FEATURE_SET_OVERRIDE <- Sys.getenv(
  "PHENOLOGY_SHAP_FEATURE_SET",
  unset = ""
)

# The primary paper analysis addresses two distinct transfer questions. Use the
# environment variable to add sensitivity schemes without pooling schemes.
# Example: PHENOLOGY_SHAP_SCHEMES="spatial_block,expanding_origin,environmental_zone"
SCHEMES_OVERRIDE <- Sys.getenv(
  "PHENOLOGY_SHAP_SCHEMES",
  unset = ""
)
DEFAULT_SHAP_SCHEMES <- c("spatial_block", "expanding_origin")

SEED <- 123L
TARGET <- "BBCH"
LEARNER_ID <- "xgboost"

read_env_numeric <- function(name, default) {
  value <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) return(as.numeric(default))
  parsed <- suppressWarnings(as.numeric(value))
  if (!is.finite(parsed)) stop(name, " must be a finite number.")
  parsed
}

read_env_integer <- function(name, default, minimum = 1L) {
  value <- read_env_numeric(name, default)
  parsed <- as.integer(value)
  if (!identical(value, as.numeric(parsed)) || parsed < minimum) {
    stop(name, " must be an integer >= ", minimum, ".")
  }
  parsed
}

# Exact refitted predictions must reproduce the stored benchmark predictions.
PREDICTION_MATCH_TOLERANCE <- read_env_numeric(
  "PHENOLOGY_SHAP_PREDICTION_TOLERANCE",
  1e-5
)
SHAP_ADDITIVITY_TOLERANCE <- read_env_numeric(
  "PHENOLOGY_SHAP_ADDITIVITY_TOLERANCE",
  2e-3
)
SHAP_ADDITIVITY_RELATIVE_TOLERANCE <- read_env_numeric(
  "PHENOLOGY_SHAP_ADDITIVITY_RELATIVE_TOLERANCE",
  1e-5
)
INTERACTION_ADDITIVITY_TOLERANCE <- read_env_numeric(
  "PHENOLOGY_SHAP_INTERACTION_ADDITIVITY_TOLERANCE",
  2e-3
)
INTERACTION_ADDITIVITY_RELATIVE_TOLERANCE <- read_env_numeric(
  "PHENOLOGY_SHAP_INTERACTION_ADDITIVITY_RELATIVE_TOLERANCE",
  1e-5
)

# XGBoost documents equality between predictions and summed TreeSHAP values as
# holding up to floating-point precision. Across many trees, a single fixed
# cutoff is brittle. Each row therefore uses an absolute-plus-relative bound:
# 2e-3 + 1e-5 * max(1, abs(prediction)). On the 0--99 BBCH scale this remains
# below 0.003 BBCH and is sensitive to any scientifically substantive failure.

# Exact interaction SHAP is substantially more expensive than ordinary SHAP.
# It is therefore evaluated on one representative per environment x exact BBCH
# stage, balanced across folds and stages, and capped within each scheme.
RUN_INTERACTIONS <- tolower(trimws(Sys.getenv(
  "PHENOLOGY_SHAP_RUN_INTERACTIONS",
  unset = "true"
))) %in% c("1", "true", "yes")
INTERACTION_MAX_ROWS_PER_SCHEME <- read_env_integer(
  "PHENOLOGY_SHAP_INTERACTION_MAX_ROWS_PER_SCHEME",
  10000L
)

# Plotting uses an environment-stage-balanced sample.
SUMMARY_PLOT_MAX_ROWS_PER_SCHEME <- read_env_integer(
  "PHENOLOGY_SHAP_SUMMARY_MAX_ROWS_PER_SCHEME",
  4000L
)
DEPENDENCE_TOP_N <- 6L
INTERACTION_TOP_N <- read_env_integer(
  "PHENOLOGY_SHAP_INTERACTION_TOP_N",
  10L
)
FOCAL_BBCH_STAGES <- c(75)

BOOTSTRAP_REPLICATES <- read_env_integer(
  "PHENOLOGY_SHAP_BOOTSTRAP_REPLICATES",
  2000L,
  minimum = 100L
)
CORRELATION_WARNING_THRESHOLD <- 0.80

# The saved learner parameters retain the benchmark's one-thread XGBoost
# setting; outer-fold jobs can run in parallel. Two workers is conservative for
# server memory.
WORKERS_REQUESTED <- read_env_integer(
  "PHENOLOGY_SHAP_WORKERS",
  3L
)

RESUME_FROM_CHECKPOINTS <- TRUE
CHECKPOINT_VERSION <- 1L

# `checkpoints_only` is a recovery mode for an interrupted plotting run. It
# never refits a model or recomputes TreeSHAP: every expected fold checkpoint
# must already exist and contain the requested feature set. It also accepts
# older checkpoint-version metadata when their row-level SHAP contents remain
# structurally compatible, which is appropriate for post-processing-only fixes.
CHECKPOINT_MODE <- tolower(trimws(Sys.getenv(
  "PHENOLOGY_SHAP_CHECKPOINT_MODE",
  unset = "strict"
)))
if (!CHECKPOINT_MODE %in% c("strict", "checkpoints_only")) {
  stop(
    "PHENOLOGY_SHAP_CHECKPOINT_MODE must be either 'strict' or ",
    "'checkpoints_only'."
  )
}

# ------------------------------------------------------------------------------
# 2. PACKAGES AND REPRODUCIBILITY
# ------------------------------------------------------------------------------

required_packages <- c(
  "data.table",
  "R6",
  "mlr3",
  "mlr3learners",
  "mlr3spatiotempcv",
  "xgboost",
  "ggplot2",
  "future",
  "future.apply",
  "parallelly",
  "mgcv"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install these packages before running SHAP: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(R6)
  library(mlr3)
  library(mlr3learners)
  library(mlr3spatiotempcv)
  library(xgboost)
  library(ggplot2)
  library(future)
  library(future.apply)
  library(parallelly)
  library(mgcv)
})

RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)
data.table::setDTthreads(1L)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

workers <- max(
  1L,
  min(WORKERS_REQUESTED, parallelly::availableCores(omit = 1L))
)

# ------------------------------------------------------------------------------
# 3. GENERAL UTILITIES
# ------------------------------------------------------------------------------

safe_id <- function(x) {
  gsub("[^A-Za-z0-9_-]+", "_", as.character(x))
}

parse_comma_vector <- function(x) {
  if (!nzchar(trimws(x))) return(character())
  answer <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  unique(answer[nzchar(answer)])
}

assert_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " does not exist: ", path)
  }
  invisible(path)
}

write_table <- function(x, filename) {
  fwrite(x, file.path(OUTPUT_DIR, filename), na = "NA")
}

atomic_save_rds <- function(object, path, compress = FALSE) {
  temporary <- paste0(path, ".tmp_", Sys.getpid())
  saveRDS(object, temporary, compress = compress)
  if (!file.rename(temporary, path)) {
    stop("Could not atomically write: ", path)
  }
  invisible(path)
}

# Recommended mlr3 learners are registered through namespace load hooks. This
# explicit check is required inside clean multisession workers as well as in the
# main R process.
ensure_xgboost_learner_registered <- function() {
  loadNamespace("mlr3learners")
  if (!"regr.xgboost" %in% mlr3::mlr_learners$keys()) {
    stop(
      "Learner 'regr.xgboost' is not registered. Installed versions: ",
      "mlr3=", as.character(utils::packageVersion("mlr3")),
      ", mlr3learners=",
      as.character(utils::packageVersion("mlr3learners")),
      ", xgboost=", as.character(utils::packageVersion("xgboost")),
      ". Restore the benchmark renv before running SHAP."
    )
  }
  invisible(TRUE)
}

ensure_xgboost_learner_registered()

phase_from_bbch <- function(x) {
  fcase(
    x <= 10, "establishment",
    x <= 29, "tillering",
    x <= 39, "stem_elongation",
    x <= 69, "heading_flowering",
    x <= 89, "grain_development_ripening",
    default = "maturity_post_maturity"
  )
}

scheme_bootstrap_cluster <- function(
  scheme,
  site_id,
  season_year,
  env_zone,
  fold_id,
  env_id
) {
  fcase(
    scheme == "spatial_block",
    paste0("site:", site_id),
    scheme %in% c("leave_one_year_out", "expanding_origin"),
    paste0("year:", season_year),
    scheme == "environmental_zone",
    paste0("zone:", env_zone),
    scheme == "joint_space_time",
    paste0("joint-fold:", fold_id),
    default = paste0("environment:", env_id)
  )
}

cor_or_na <- function(x, y, method = "pearson") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L) return(NA_real_)
  if (sd(x[keep]) == 0 || sd(y[keep]) == 0) return(NA_real_)
  unname(cor(x[keep], y[keep], method = method))
}

cluster_bootstrap_mean <- function(
  unit_table,
  value_column,
  cluster_column,
  replicates,
  seed
) {
  cluster_summary <- unit_table[, .(
    value_sum = sum(get(value_column)),
    n_units = .N
  ), by = cluster_column]

  observed <- mean(unit_table[[value_column]])

  if (nrow(cluster_summary) < 2L) {
    return(data.table(
      estimate = observed,
      ci_low = NA_real_,
      ci_high = NA_real_,
      n_clusters = nrow(cluster_summary),
      n_units = nrow(unit_table)
    ))
  }

  set.seed(seed)
  bootstrap_values <- replicate(replicates, {
    sampled <- sample.int(
      nrow(cluster_summary),
      nrow(cluster_summary),
      replace = TRUE
    )
    sum(cluster_summary$value_sum[sampled]) /
      sum(cluster_summary$n_units[sampled])
  })

  data.table(
    estimate = observed,
    ci_low = as.numeric(quantile(bootstrap_values, 0.025)),
    ci_high = as.numeric(quantile(bootstrap_values, 0.975)),
    n_clusters = nrow(cluster_summary),
    n_units = nrow(unit_table)
  )
}

# Return one representative row per environment x exact BBCH unit, then cap
# with round-robin selection over fold x stage strata. This preserves broad
# fold/stage coverage when a fixed plotting or interaction budget is required.
balanced_environment_stage_sample <- function(
  metadata,
  maximum_rows,
  seed
) {
  candidates <- copy(metadata)
  candidates[, local_row := .I]

  set.seed(seed)
  representatives <- candidates[
    ,
    .(local_row = sample(local_row, 1L)),
    by = .(env_id, BBCH_stage)
  ]
  representatives <- candidates[representatives$local_row]

  if (nrow(representatives) <= maximum_rows) {
    representatives[, local_row := NULL]
    return(representatives[])
  }

  representatives[, random_order := runif(.N)]
  representatives[
    ,
    stratum_round := frank(random_order, ties.method = "first"),
    by = .(fold_id, BBCH_stage)
  ]
  setorder(
    representatives,
    stratum_round,
    random_order,
    fold_id,
    BBCH_stage
  )
  representatives <- representatives[seq_len(maximum_rows)]
  representatives[
    ,
    c("local_row", "random_order", "stratum_round") := NULL
  ]
  representatives[]
}

balanced_plot_sample <- function(shap_table, maximum_rows, seed) {
  rbindlist(lapply(seq_along(unique(shap_table$scheme)), function(i) {
    scheme_id <- unique(shap_table$scheme)[i]
    balanced_environment_stage_sample(
      metadata = shap_table[scheme == scheme_id],
      maximum_rows = maximum_rows,
      seed = seed + i
    )
  }), use.names = TRUE)
}

# ------------------------------------------------------------------------------
# 4. BIOLOGICAL FEATURE GROUPS AND PRE-SPECIFIED INTERACTIONS
# ------------------------------------------------------------------------------

feature_group_lookup <- c(
  cumDTU_beta = "developmental_clock",
  cumdevfun_phase = "developmental_clock",
  cumdevfun_phase_water_0_28 = "developmental_clock",
  cumdevfun_phase_water_0_100 = "developmental_clock",
  sowing_DOY = "sowing_context",
  photoperiod_civil = "photoperiod",
  Daylength = "photoperiod",
  cumRR = "water_supply",
  sumRR_30 = "water_supply",
  cumRad = "radiation",
  sumRad_30 = "radiation",
  sumRad_14 = "radiation",
  cumCDD = "cold_exposure",
  sumCDD_30 = "cold_exposure",
  sumCDD_14 = "cold_exposure",
  cumHDD = "heat_exposure",
  sumHDD_30 = "heat_exposure",
  sumHDD_14 = "heat_exposure",
  paw_fc33_0_28_mm = "soil_water_capacity",
  paw_fc33_0_100_mm = "soil_water_capacity",
  FTSW_0_28_proxy = "soil_water_status",
  FTSW_0_100_proxy = "soil_water_status",
  meanSoilWater_0_28_14 = "soil_water_status",
  meanSoilWater_0_100_30 = "soil_water_status",
  soil_water_0_28_vwc = "soil_water_status",
  soil_water_0_100_vwc = "soil_water_status",
  sumPEvap_30 = "atmospheric_water_demand",
  sumPEvap_14 = "atmospheric_water_demand",
  soilWaterChange_14 = "soil_water_change",
  climaticWaterBalance_14 = "climatic_water_balance",
  climaticWaterBalance_30 = "climatic_water_balance"
)

predeclared_interactions <- data.table(
  feature_1 = c(
    "cumGDD_from_sow_Soltani",
    "cumGDD_from_sow_Soltani",
    "cumGDD_from_sow_Soltani",
    "cumGDD_from_sow_Soltani",
    "cumGDD_from_sow_Soltani",
    "cumRR",
    "cumdevfun_phase",
    "cumdevfun_phase",
    "cumdevfun_phase",
    "cumdevfun_phase",
    "cumdevfun_phase",
    "FTSW_0_28_proxy"
  ),
  feature_2 = c(
    "sowing_DOY",
    "Daylength",
    "cumRR",
    "cumHDD",
    "paw_fc33_0_28_mm",
    "paw_fc33_0_28_mm",
    "sowing_DOY",
    "FTSW_0_28_proxy",
    "sumHDD_14",
    "sumCDD_14",
    "climaticWaterBalance_14",
    "paw_fc33_0_28_mm"
  ),
  hypothesis = c(
    "thermal progression x sowing timing",
    "thermal progression x photoperiod",
    "thermal progression x accumulated water supply",
    "thermal progression x accumulated heat exposure",
    "thermal progression x soil water capacity",
    "water supply x soil water capacity",
    "process clock x sowing timing",
    "process clock x soil water status",
    "process clock x recent heat exposure",
    "process clock x recent cold exposure",
    "process clock x recent climatic water balance",
    "soil water status x soil water capacity"
  )
)

canonical_pair <- function(feature_1, feature_2) {
  paste(
    pmin(feature_1, feature_2),
    pmax(feature_1, feature_2),
    sep = "__X__"
  )
}

predeclared_interactions[
  ,
  pair_id := canonical_pair(feature_1, feature_2)
]

# ------------------------------------------------------------------------------
# 5. LOAD AND VALIDATE THE BENCHMARK
# ------------------------------------------------------------------------------

configuration_path <- file.path(
  BENCHMARK_DIR,
  "benchmark_configuration.rds"
)
fold_manifest_path <- file.path(
  BENCHMARK_DIR,
  "outer_fold_manifest.csv"
)
split_specs_path <- file.path(
  BENCHMARK_DIR,
  "outer_split_specs.csv"
)
oof_predictions_path <- file.path(
  BENCHMARK_DIR,
  "oof_predictions.csv"
)
model_specifications_path <- file.path(
  BENCHMARK_DIR,
  "fold_model_specifications.rds"
)

assert_file(configuration_path, "Benchmark configuration")
assert_file(fold_manifest_path, "Outer-fold manifest")
assert_file(split_specs_path, "Outer-split specifications")
assert_file(oof_predictions_path, "OOF predictions")
assert_file(model_specifications_path, "Fold model specifications")

benchmark_configuration <- readRDS(configuration_path)

if (
  is.null(benchmark_configuration$selected_learner) ||
    !identical(
      tolower(as.character(benchmark_configuration$selected_learner)),
      LEARNER_ID
    )
) {
  stop(
    "The generalisation benchmark did not select XGBoost. Saved learner: ",
    as.character(benchmark_configuration$selected_learner)
  )
}

feature_set <- if (nzchar(FEATURE_SET_OVERRIDE)) {
  FEATURE_SET_OVERRIDE
} else {
  as.character(benchmark_configuration$feature_set)
}

if (!identical(feature_set, as.character(benchmark_configuration$feature_set))) {
  stop(
    "The SHAP feature-set override differs from the frozen generalisation ",
    "feature set. SHAP reconstruction requires: ",
    as.character(benchmark_configuration$feature_set)
  )
}
features <- as.character(benchmark_configuration$features)
if (!length(features)) {
  stop("The benchmark configuration contains no frozen predictor list.")
}
predeclared_interactions <- predeclared_interactions[
  feature_1 %in% features & feature_2 %in% features
]

missing_group_mappings <- setdiff(
  features,
  names(feature_group_lookup)
)
if (length(missing_group_mappings)) {
  stop(
    "Add biological group mappings for: ",
    paste(missing_group_mappings, collapse = ", ")
  )
}
feature_groups <- data.table(
  feature = features,
  biological_group = unname(feature_group_lookup[features])
)
write_table(feature_groups, "feature_group_definitions.csv")

requested_schemes <- parse_comma_vector(SCHEMES_OVERRIDE)
if (!length(requested_schemes)) {
  requested_schemes <- DEFAULT_SHAP_SCHEMES
}

fold_manifest <- fread(fold_manifest_path)
outer_split_specs <- fread(split_specs_path)
benchmark_oof <- fread(oof_predictions_path)
fold_model_specifications <- readRDS(model_specifications_path)

available_schemes <- intersect(
  as.character(benchmark_configuration$cv_schemes),
  unique(as.character(outer_split_specs$scheme))
)
missing_schemes <- setdiff(requested_schemes, available_schemes)
if (length(missing_schemes)) {
  stop(
    "Requested SHAP scheme(s) are absent from the benchmark: ",
    paste(missing_schemes, collapse = ", ")
  )
}
schemes <- requested_schemes
selected_feature_set_value <- feature_set

outer_split_specs <- outer_split_specs[scheme %in% schemes]
outer_split_specs[, scheme_order__ := match(scheme, schemes)]
setorder(
  outer_split_specs,
  scheme_order__,
  split_id
)
outer_split_specs[, scheme_order__ := NULL]
if (!nrow(outer_split_specs)) stop("No requested outer splits remain.")

# An archived run may record the original server path. Use the repository copy
# if that path is unavailable; the MD5 check below still requires identical data.
recorded_input_file <- benchmark_configuration$input_file
if (is.null(recorded_input_file) || length(recorded_input_file) != 1L) {
  stop("The benchmark configuration does not specify one input file.")
}
repo_input_file <- file.path(PROJECT_DIR, "Data", "pheno_df_CV.csv")
default_input_file <- if (file.exists(recorded_input_file)) {
  recorded_input_file
} else {
  repo_input_file
}
input_file <- Sys.getenv("PHENOLOGY_INPUT_FILE", unset = default_input_file)
assert_file(input_file, "Benchmark input data")

input_md5 <- unname(tools::md5sum(input_file))
if (
  !is.null(benchmark_configuration$input_md5) &&
    !identical(input_md5, benchmark_configuration$input_md5)
) {
  stop(
    "Input MD5 differs from the generalisation benchmark. ",
    "SHAP must use the exact benchmark cohort."
  )
}

reconstruction_packages <- c("mlr3", "mlr3learners", "xgboost")
benchmark_package_versions <- vapply(
  reconstruction_packages,
  function(package) {
    value <- benchmark_configuration$package_versions[[package]]
    if (is.null(value)) NA_character_ else as.character(value)
  },
  character(1)
)
current_package_versions <- vapply(
  reconstruction_packages,
  function(package) as.character(utils::packageVersion(package)),
  character(1)
)
package_version_audit <- data.table(
  package = reconstruction_packages,
  benchmark_version = unname(benchmark_package_versions),
  current_version = unname(current_package_versions)
)
package_version_audit[
  ,
  exact_match := is.na(benchmark_version) |
    benchmark_version == current_version
]
write_table(package_version_audit, "shap_package_version_audit.csv")

known_version_mismatches <- package_version_audit[
  !is.na(benchmark_version) & !exact_match
]
if (nrow(known_version_mismatches)) {
  mismatch_text <- paste0(
    known_version_mismatches$package,
    " (benchmark ",
    known_version_mismatches$benchmark_version,
    ", current ",
    known_version_mismatches$current_version,
    ")",
    collapse = "; "
  )
  stop(
    "Exact fold reconstruction requires the benchmark package versions. ",
    "Version mismatch: ",
    mismatch_text,
    ". Restore the benchmark renv before running SHAP."
  )
}

benchmark_xgboost_version <- unname(
  benchmark_package_versions[["xgboost"]]
)
current_xgboost_version <- unname(current_package_versions[["xgboost"]])

raw_data <- fread(
  input_file,
  na.strings = c("", "NA", "NaN", "Inf", "-Inf"),
  showProgress = TRUE
)
raw_data[, source_row_id_internal := .I]

if (anyDuplicated(fold_manifest$source_row_id)) {
  stop("source_row_id is duplicated in the fold manifest.")
}
if (anyDuplicated(fold_manifest$obs_id)) {
  stop("obs_id is duplicated in the fold manifest.")
}
if (
  any(
    fold_manifest$source_row_id < 1L |
      fold_manifest$source_row_id > nrow(raw_data)
  )
) {
  stop("Fold manifest contains invalid source_row_id values.")
}

source_match <- match(
  fold_manifest$source_row_id,
  raw_data$source_row_id_internal
)
if (anyNA(source_match)) {
  stop("Could not reconstruct all manifest rows from source_row_id.")
}
cohort <- raw_data[source_match]
cohort[, obs_id := fold_manifest$obs_id]
cohort[, source_row_id := fold_manifest$source_row_id]
cohort[, env_id := as.character(fold_manifest$env_id)]
cohort[, site_id := as.character(fold_manifest$site_id)]
cohort[, season_year := as.integer(fold_manifest$season_year)]
cohort[, env_zone := as.character(fold_manifest$env_zone)]
cohort[, observation_date := as.IDate(fold_manifest$Date)]
cohort[, longitude := as.numeric(fold_manifest$lon)]
cohort[, latitude := as.numeric(fold_manifest$lat)]
cohort[, spatial_block := as.character(fold_manifest$spatial_block)]
cohort[
  ,
  joint_temporal_block := as.character(
    fold_manifest$joint_temporal_block
  )
]

missing_columns <- setdiff(c(TARGET, features), names(cohort))
if (length(missing_columns)) {
  stop(
    "Required data columns are absent: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (!isTRUE(all.equal(
  as.numeric(cohort[[TARGET]]),
  as.numeric(fold_manifest$BBCH),
  tolerance = 0
))) {
  stop("Target values do not match the saved fold manifest.")
}

if (anyNA(cohort[, c(TARGET, features), with = FALSE])) {
  stop("The reconstructed benchmark cohort contains missing model values.")
}
if (
  any(!vapply(
    cohort[, ..features],
    is.numeric,
    logical(1)
  ))
) {
  stop("All XGBoost predictors must be numeric.")
}
if (
  any(!vapply(
    cohort[, c(TARGET, features), with = FALSE],
    function(x) all(is.finite(x)),
    logical(1)
  ))
) {
  stop("The reconstructed benchmark cohort contains non-finite values.")
}

# Reconstruct the exact mlr3 task used by the generalisation benchmark. The
# learner trains on task$feature_names, whose order need not equal the order in
# benchmark_configuration$features. That distinction changes stochastic column
# sampling and must therefore be preserved for an exact booster refit.
cohort[, `.obs_id` := as.integer(obs_id)]
cohort[, `.env_id` := as.character(env_id)]
cohort[, `.site_id` := as.character(site_id)]
cohort[, `.lon` := as.numeric(longitude)]
cohort[, `.lat` := as.numeric(latitude)]

benchmark_backend_columns <- unique(c(
  ".obs_id",
  ".env_id",
  ".site_id",
  TARGET,
  ".lon",
  ".lat",
  features
))
benchmark_backend <- cohort[, ..benchmark_backend_columns]
benchmark_task <- TaskRegrST$new(
  id = feature_set,
  backend = benchmark_backend,
  target = TARGET,
  coordinate_names = c(".lon", ".lat")
)
benchmark_task$select(features)

if (!identical(
  as.integer(benchmark_task$row_ids),
  as.integer(cohort$obs_id)
)) {
  stop(
    "Reconstructed mlr3 task row IDs differ from the benchmark manifest."
  )
}
model_features <- as.character(benchmark_task$feature_names)
if (!setequal(model_features, features)) {
  stop("Reconstructed mlr3 task does not contain the frozen feature set.")
}
configured_feature_index <- match(features, model_features)
if (anyNA(configured_feature_index)) {
  stop("Could not map model features to the configured feature order.")
}
model_feature_order <- data.table(
  model_position = seq_along(model_features),
  feature = model_features,
  configured_position = match(model_features, features),
  same_position_as_configuration = model_features == features
)
write_table(model_feature_order, "shap_model_feature_order.csv")
if (!identical(model_features, features)) {
  message(
    "The mlr3 model-feature order differs from the configured feature order. ",
    "The mlr3 order will be used for refitting and TreeSHAP, then outputs ",
    "will be restored to the configured order."
  )
}

benchmark_oof <- benchmark_oof[
  get("feature_set") == selected_feature_set_value &
    get("learner_id") == LEARNER_ID &
    get("scheme") %in% schemes
]
if (!nrow(benchmark_oof)) {
  stop("No matching XGBoost OOF predictions were found.")
}
expected_oof_splits <- sort(unique(as.integer(outer_split_specs$split_id)))
observed_oof_splits <- sort(unique(as.integer(benchmark_oof$split_id)))
if (!identical(expected_oof_splits, observed_oof_splits)) {
  stop("OOF predictions do not cover exactly the requested outer splits.")
}
if (benchmark_oof[, anyDuplicated(row_ids), by = scheme][, any(V1 > 0L)]) {
  stop("An observation has duplicate XGBoost OOF predictions within a scheme.")
}

model_spec_key <- function(specification) {
  if (is.null(specification$job)) return(NA_character_)
  paste(
    as.character(specification$job$model_id),
    as.character(specification$job$scheme),
    as.integer(specification$job$split_id),
    sep = "__"
  )
}

saved_model_keys <- vapply(
  fold_model_specifications,
  model_spec_key,
  character(1)
)
required_model_keys <- paste(
  LEARNER_ID,
  outer_split_specs$scheme,
  outer_split_specs$split_id,
  sep = "__"
)
if (anyDuplicated(saved_model_keys[!is.na(saved_model_keys)])) {
  stop("Fold model specifications contain duplicated model/scheme/split keys.")
}
missing_model_keys <- setdiff(required_model_keys, saved_model_keys)
if (length(missing_model_keys)) {
  stop(
    "Missing saved XGBoost specification(s): ",
    paste(missing_model_keys, collapse = ", ")
  )
}

scheme_claims <- benchmark_configuration$scheme_claims
if (is.null(scheme_claims)) scheme_claims <- list()
scheme_estimands <- data.table(
  scheme = schemes,
  transfer_claim = vapply(schemes, function(x) {
    value <- scheme_claims[[x]]
    if (is.null(value)) NA_character_ else as.character(value)
  }, character(1)),
  bootstrap_cluster = c(
    spatial_block = "physical site",
    leave_one_year_out = "season year",
    expanding_origin = "season year",
    environmental_zone = "environmental zone",
    joint_space_time = "joint outer space-time fold"
  )[schemes]
)
write_table(scheme_estimands, "shap_scheme_estimands.csv")

interaction_sample_manifest <- data.table()
if (RUN_INTERACTIONS) {
  interaction_candidates <- benchmark_oof[, .(
    scheme = as.character(scheme),
    split_id = as.integer(split_id),
    fold_id = as.character(fold_id),
    obs_id = as.integer(row_ids),
    source_row_id = as.integer(source_row_id),
    env_id = as.character(env_id),
    site_id = as.character(site_id),
    season_year = as.integer(season_year),
    env_zone = as.character(env_zone),
    BBCH_stage = as.numeric(truth)
  )]
  interaction_sample_manifest <- rbindlist(lapply(
    seq_along(schemes),
    function(i) {
      scheme_id <- schemes[i]
      balanced_environment_stage_sample(
        metadata = interaction_candidates[scheme == scheme_id],
        maximum_rows = INTERACTION_MAX_ROWS_PER_SCHEME,
        seed = SEED + 10000L + i
      )
    }
  ), use.names = TRUE)
  setorder(interaction_sample_manifest, scheme, split_id, obs_id)
  write_table(
    interaction_sample_manifest,
    "shap_interaction_sample_manifest.csv"
  )
}

# ------------------------------------------------------------------------------
# 6. XGBOOST PARAMETER AND SHAP UTILITIES
# ------------------------------------------------------------------------------

make_refit_xgboost_learner <- function(parameter_values) {
  if (is.null(parameter_values) || !length(parameter_values)) {
    stop("Saved fold specification contains no XGBoost parameters.")
  }
  values <- as.list(parameter_values)
  if (is.null(names(values)) || any(!nzchar(names(values)))) {
    stop("Saved XGBoost parameters must be a fully named list.")
  }
  if (anyDuplicated(names(values))) {
    stop("Saved XGBoost parameter names are duplicated.")
  }

  ensure_xgboost_learner_registered()
  learner <- mlr3::lrn(
    "regr.xgboost",
    predict_type = "response"
  )
  unknown_parameters <- setdiff(
    names(values),
    learner$param_set$ids()
  )
  if (length(unknown_parameters)) {
    stop(
      "The current mlr3learners wrapper does not recognise saved parameter(s): ",
      paste(unknown_parameters, collapse = ", "),
      ". Restore the benchmark renv."
    )
  }

  learner$param_set$values <- values
  learner$predict_type <- "response"
  learner
}

parameter_signature <- function(parameter_values) {
  values <- as.list(parameter_values)
  values <- values[order(names(values))]
  paste(capture.output(dput(values)), collapse = "\n")
}

normalise_interaction_array <- function(
  interaction_result,
  n_rows,
  n_features
) {
  expected <- c(
    as.integer(n_rows),
    as.integer(n_features + 1L),
    as.integer(n_features + 1L)
  )
  observed <- dim(interaction_result)
  if (identical(observed, expected)) return(interaction_result)

  # Recent XGBoost versions can retain a singleton output-group dimension
  # under strict output shaping. Squeeze only that documented singleton.
  if (
    length(observed) == 4L &&
      observed[2L] == 1L &&
      identical(observed[c(1L, 3L, 4L)], expected)
  ) {
    return(array(
      interaction_result[, 1L, , , drop = FALSE],
      dim = expected
    ))
  }

  stop(
    "Unexpected predinteraction shape: ",
    paste(observed, collapse = " x "),
    "; expected ",
    paste(expected, collapse = " x "),
    " (optionally with a singleton output-group dimension)"
  )
}

normalise_contribution_matrix <- function(
  contribution_result,
  n_rows,
  n_features
) {
  expected <- c(as.integer(n_rows), as.integer(n_features + 1L))
  observed <- dim(contribution_result)
  if (identical(observed, expected)) return(contribution_result)

  if (
    length(observed) == 3L &&
      observed[2L] == 1L &&
      identical(observed[c(1L, 3L)], expected)
  ) {
    return(matrix(
      contribution_result[, 1L, , drop = FALSE],
      nrow = n_rows,
      ncol = n_features + 1L
    ))
  }

  stop(
    "Unexpected predcontrib shape: ",
    paste(observed, collapse = " x "),
    "; expected ",
    paste(expected, collapse = " x "),
    " (optionally with a singleton output-group dimension)"
  )
}

validate_contribution_shape <- function(
  contribution_matrix,
  n_rows,
  n_features
) {
  if (
    !is.matrix(contribution_matrix) ||
      !identical(
        dim(contribution_matrix),
        c(as.integer(n_rows), as.integer(n_features + 1L))
      )
  ) {
    stop(
      "Unexpected predcontrib shape: ",
      paste(dim(contribution_matrix), collapse = " x "),
      "; expected ",
      n_rows,
      " x ",
      n_features + 1L
    )
  }
}

interaction_array_to_long <- function(
  interaction_array,
  feature_matrix,
  metadata,
  features
) {
  n <- nrow(metadata)
  p <- length(features)

  if (!identical(dim(feature_matrix), c(as.integer(n), as.integer(p)))) {
    stop("Interaction feature matrix does not match metadata and features.")
  }

  if (
    length(dim(interaction_array)) != 3L ||
      !identical(
        dim(interaction_array),
        c(as.integer(n), as.integer(p + 1L), as.integer(p + 1L))
      )
  ) {
    stop(
      "Unexpected predinteraction shape: ",
      paste(dim(interaction_array), collapse = " x "),
      "; expected ",
      paste(c(n, p + 1L, p + 1L), collapse = " x ")
    )
  }

  pair_indices <- which(
    upper.tri(matrix(FALSE, p, p), diag = FALSE),
    arr.ind = TRUE
  )

  rbindlist(lapply(seq_len(nrow(pair_indices)), function(i) {
    feature_i <- pair_indices[i, 1L]
    feature_j <- pair_indices[i, 2L]

    # The interaction array is symmetric. Multiplying one off-diagonal cell by
    # two yields the contribution of the unordered pair to the prediction.
    data.table(
      obs_id = metadata$obs_id,
      source_row_id = metadata$source_row_id,
      env_id = metadata$env_id,
      site_id = metadata$site_id,
      season_year = metadata$season_year,
      env_zone = metadata$env_zone,
      observation_date = metadata$observation_date,
      BBCH_stage = metadata$BBCH_stage,
      BBCH_phase = metadata$BBCH_phase,
      feature_1 = features[feature_i],
      feature_2 = features[feature_j],
      pair_id = canonical_pair(
        features[feature_i],
        features[feature_j]
      ),
      feature_1_value = feature_matrix[, feature_i],
      feature_2_value = feature_matrix[, feature_j],
      interaction_value =
        2 * interaction_array[, feature_i, feature_j]
    )
  }))
}

feature_support_audit <- function(
  train_matrix,
  test_matrix,
  scheme,
  split_id,
  fold_id,
  features
) {
  rbindlist(lapply(seq_along(features), function(i) {
    train_x <- train_matrix[, i]
    test_x <- test_matrix[, i]
    train_q01 <- as.numeric(quantile(train_x, 0.01))
    train_q99 <- as.numeric(quantile(train_x, 0.99))

    data.table(
      scheme = scheme,
      split_id = split_id,
      fold_id = fold_id,
      feature = features[i],
      train_min = min(train_x),
      train_q01 = train_q01,
      train_median = median(train_x),
      train_q99 = train_q99,
      train_max = max(train_x),
      test_min = min(test_x),
      test_median = median(test_x),
      test_max = max(test_x),
      n_test = length(test_x),
      n_outside_train_range = sum(
        test_x < min(train_x) | test_x > max(train_x)
      ),
      proportion_outside_train_range = mean(
        test_x < min(train_x) | test_x > max(train_x)
      ),
      proportion_outside_train_q01_q99 = mean(
        test_x < train_q01 | test_x > train_q99
      )
    )
  }))
}

# ------------------------------------------------------------------------------
# 7. REFIT EACH OUTER MODEL AND EXPLAIN ONLY ITS TEST FOLD
# ------------------------------------------------------------------------------

jobs <- copy(outer_split_specs)
jobs[, scheme_order__ := match(scheme, schemes)]
setorder(jobs, scheme_order__, split_id)
jobs[, scheme_order__ := NULL]
jobs[, job_id := .I]
setcolorder(jobs, c("job_id", setdiff(names(jobs), "job_id")))
write_table(jobs, "shap_jobs.csv")

resolve_outer_rows <- function(spec) {
  scheme <- as.character(spec$scheme)

  if (identical(scheme, "spatial_block")) {
    test_mask <- cohort$spatial_block == as.character(
      spec$holdout_spatial_fold
    )
    train_mask <- !test_mask
  } else if (identical(scheme, "leave_one_year_out")) {
    test_mask <- cohort$season_year == as.integer(spec$holdout_year)
    train_mask <- !test_mask
  } else if (identical(scheme, "expanding_origin")) {
    test_mask <- cohort$season_year == as.integer(spec$holdout_year)
    train_mask <- cohort$season_year < as.integer(spec$holdout_year)
  } else if (identical(scheme, "environmental_zone")) {
    test_mask <- cohort$env_zone == as.character(spec$holdout_zone)
    train_mask <- !test_mask
  } else if (identical(scheme, "joint_space_time")) {
    heldout_space <- as.character(spec$holdout_spatial_fold)
    heldout_time <- as.character(spec$holdout_temporal_block)
    test_mask <- (
      cohort$spatial_block == heldout_space &
        cohort$joint_temporal_block == heldout_time
    )
    train_mask <- (
      cohort$spatial_block != heldout_space &
        cohort$joint_temporal_block != heldout_time
    )
  } else {
    stop("Unknown outer scheme: ", scheme)
  }

  train_mask[is.na(train_mask)] <- FALSE
  test_mask[is.na(test_mask)] <- FALSE
  list(
    train_rows = which(train_mask),
    test_rows = which(test_mask),
    excluded_rows = which(!train_mask & !test_mask)
  )
}

get_saved_fold_specification <- function(job) {
  key <- paste(
    LEARNER_ID,
    as.character(job$scheme),
    as.integer(job$split_id),
    sep = "__"
  )
  index <- match(key, saved_model_keys)
  if (is.na(index)) stop("No saved fold model specification for ", key)
  specification <- fold_model_specifications[[index]]
  if (!identical(as.character(specification$status), "ok")) {
    stop(
      "Saved fold model did not finish successfully for ",
      key,
      ": ",
      as.character(specification$error)
    )
  }
  if (is.null(specification$tuned_parameter_values)) {
    stop("Saved fold model contains no tuned parameters for ", key)
  }
  specification
}

run_shap_job <- function(job) {
  job_scheme <- as.character(job$scheme)
  job_split_id <- as.integer(job$split_id)
  job_fold_id <- as.character(job$fold_id)
  saved_specification <- get_saved_fold_specification(job)
  tuned_parameter_values <- saved_specification$tuned_parameter_values
  benchmark_job_id <- as.integer(saved_specification$job$job_id)
  if (length(benchmark_job_id) != 1L || is.na(benchmark_job_id)) {
    stop("Saved fold specification contains no valid benchmark job ID.")
  }
  tuned_parameter_signature <- parameter_signature(tuned_parameter_values)
  split_signature <- paste(capture.output(dput(as.list(job))), collapse = "\n")
  selected_interaction_ids <- if (RUN_INTERACTIONS) {
    interaction_sample_manifest[
      scheme == job_scheme & split_id == job_split_id,
      obs_id
    ]
  } else {
    integer()
  }
  interaction_signature <- paste(
    RUN_INTERACTIONS,
    INTERACTION_MAX_ROWS_PER_SCHEME,
    paste(selected_interaction_ids, collapse = ","),
    sep = "__"
  )

  row_sets <- resolve_outer_rows(job)
  train_rows <- row_sets$train_rows
  test_rows <- row_sets$test_rows
  excluded_rows <- row_sets$excluded_rows

  checkpoint_file <- file.path(
    CHECKPOINT_DIR,
    paste0(
      safe_id(feature_set),
      "__",
      safe_id(job_scheme),
      "__split_",
      sprintf("%03d", job_split_id),
      "__fold_",
      safe_id(job_fold_id),
      ".rds"
    )
  )

  if (RESUME_FROM_CHECKPOINTS && file.exists(checkpoint_file)) {
    cached <- tryCatch(
      readRDS(checkpoint_file),
      error = function(error) NULL
    )
    if (
      !is.null(cached) &&
        identical(cached$checkpoint_version, CHECKPOINT_VERSION) &&
        identical(cached$input_md5, input_md5) &&
        identical(cached$feature_set, feature_set) &&
        identical(cached$scheme, job_scheme) &&
        identical(cached$split_id, job_split_id) &&
        identical(cached$fold_id, job_fold_id) &&
        identical(cached$benchmark_job_id, benchmark_job_id) &&
        identical(cached$split_signature, split_signature) &&
        identical(cached$interaction_signature, interaction_signature) &&
        identical(
          cached$tuned_parameter_signature,
          tuned_parameter_signature
        ) &&
        identical(cached$xgboost_version, as.character(
          utils::packageVersion("xgboost")
        ))
    ) {
      return(cached)
    }
  }

  if (!length(train_rows) || !length(test_rows)) {
    stop(
      "Empty training or test set for ",
      job_scheme,
      " fold ",
      job_fold_id
    )
  }
  if (length(intersect(
    unique(cohort$env_id[train_rows]),
    unique(cohort$env_id[test_rows])
  ))) {
    stop(
      "Environment leakage in ",
      job_scheme,
      " fold ",
      job_fold_id
    )
  }
  if (
    job_scheme %in% c("spatial_block", "joint_space_time") &&
      length(intersect(
        unique(cohort$site_id[train_rows]),
        unique(cohort$site_id[test_rows])
      ))
  ) {
    stop(
      "Physical-site leakage in ",
      job_scheme,
      " fold ",
      job_fold_id
    )
  }
  if (
    job_scheme %in% c(
      "leave_one_year_out",
      "expanding_origin",
      "joint_space_time"
    ) &&
      length(intersect(
        unique(cohort$season_year[train_rows]),
        unique(cohort$season_year[test_rows])
      ))
  ) {
    stop("Season-year leakage in ", job_scheme, " fold ", job_fold_id)
  }
  if (
    identical(job_scheme, "environmental_zone") &&
      length(intersect(
        unique(cohort$env_zone[train_rows]),
        unique(cohort$env_zone[test_rows])
      ))
  ) {
    stop("Environmental-zone leakage in fold ", job_fold_id)
  }
  if (
    identical(job_scheme, "expanding_origin") &&
      max(cohort$season_year[train_rows]) >=
        min(cohort$season_year[test_rows])
  ) {
    stop("Expanding-origin training is not strictly earlier than its test year.")
  }

  train_ids <- as.integer(cohort$obs_id[train_rows])
  test_ids <- as.integer(cohort$obs_id[test_rows])

  fold_train_task <- benchmark_task$clone(deep = TRUE)
  fold_train_task$filter(train_ids)
  if (!identical(
    as.character(fold_train_task$feature_names),
    model_features
  )) {
    stop("Fold task changed the benchmark model-feature order.")
  }

  # Refit through the same mlr3learners wrapper as the benchmark. Replacing
  # this with a direct xgb.train() call is not equivalent when task feature
  # order or wrapper-managed training arguments differ.
  set.seed(SEED + benchmark_job_id)
  fold_learner <- make_refit_xgboost_learner(tuned_parameter_values)
  fold_learner$train(fold_train_task)
  booster <- fold_learner$model
  if (!inherits(booster, "xgb.Booster")) {
    stop("The refitted mlr3 learner did not expose an xgb.Booster model.")
  }

  wrapper_prediction <- fold_learner$predict(
    benchmark_task,
    row_ids = test_ids
  )
  if (!identical(
    as.integer(wrapper_prediction$row_ids),
    test_ids
  )) {
    stop("Refitted mlr3 prediction rows do not match the outer test IDs.")
  }
  prediction <- as.numeric(wrapper_prediction$response)

  # Use the task-defined model order for both support auditing and TreeSHAP.
  # The final output is reordered to the configured scientific feature order.
  train_matrix <- as.matrix(
    fold_train_task$data(cols = model_features)
  )
  test_matrix <- as.matrix(
    benchmark_task$data(rows = test_ids, cols = model_features)
  )
  storage.mode(train_matrix) <- "double"
  storage.mode(test_matrix) <- "double"

  dtest <- xgb.DMatrix(
    data = test_matrix,
    label = cohort[[TARGET]][test_rows],
    missing = NA_real_
  )
  booster_prediction <- as.numeric(predict(
    booster,
    dtest,
    strict_shape = FALSE
  ))
  maximum_wrapper_booster_difference <- max(abs(
    prediction - booster_prediction
  ))
  if (
    !is.finite(maximum_wrapper_booster_difference) ||
      maximum_wrapper_booster_difference > PREDICTION_MATCH_TOLERANCE
  ) {
    stop(
      "The mlr3 wrapper and its extracted booster disagree in fold ",
      job_fold_id,
      ". Maximum absolute difference = ",
      signif(maximum_wrapper_booster_difference, 8)
    )
  }

  benchmark_reference <- benchmark_oof[
    get("scheme") == job_scheme &
      get("split_id") == job_split_id
  ]
  benchmark_reference <- benchmark_reference[
    match(cohort$obs_id[test_rows], row_ids)
  ]
  if (
    nrow(benchmark_reference) != length(test_rows) ||
      anyNA(benchmark_reference$row_ids) ||
      !identical(
        as.integer(benchmark_reference$row_ids),
        as.integer(cohort$obs_id[test_rows])
      )
  ) {
    stop(
      "Could not align refitted predictions to benchmark OOF rows for ",
      job_scheme,
      " fold ",
      job_fold_id
    )
  }
  if (!isTRUE(all.equal(
    as.numeric(benchmark_reference$truth),
    as.numeric(cohort[[TARGET]][test_rows]),
    tolerance = 0
  ))) {
    stop("Benchmark truth values differ in fold ", job_fold_id)
  }

  maximum_prediction_difference <- max(
    abs(prediction - benchmark_reference$response)
  )
  if (
    !is.finite(maximum_prediction_difference) ||
      maximum_prediction_difference > PREDICTION_MATCH_TOLERANCE
  ) {
    mismatch_diagnostic <- data.table(
      obs_id = test_ids,
      benchmark_prediction = as.numeric(benchmark_reference$response),
      refitted_prediction = prediction
    )
    mismatch_diagnostic[
      ,
      difference := refitted_prediction - benchmark_prediction
    ]
    mismatch_diagnostic[, absolute_difference := abs(difference)]
    setorder(mismatch_diagnostic, -absolute_difference)
    mismatch_path <- file.path(
      OUTPUT_DIR,
      paste0(
        "prediction_reconstruction_mismatch__",
        safe_id(job_scheme),
        "__",
        safe_id(job_fold_id),
        ".csv"
      )
    )
    fwrite(mismatch_diagnostic, mismatch_path, na = "NA")

    stop(
      "Refitted XGBoost predictions do not reproduce the benchmark for ",
      job_scheme,
      " fold ",
      job_fold_id,
      ". Maximum absolute difference = ",
      signif(maximum_prediction_difference, 8),
      ". Reconstruction used the benchmark mlr3 wrapper, task row IDs, ",
      "task feature order, saved parameter list, and package versions. ",
      "Diagnostic: ",
      mismatch_path
    )
  }

  contribution <- predict(
    booster,
    dtest,
    predcontrib = TRUE,
    approxcontrib = FALSE,
    strict_shape = FALSE
  )
  contribution <- normalise_contribution_matrix(
    contribution_result = contribution,
    n_rows = length(test_rows),
    n_features = length(model_features)
  )
  validate_contribution_shape(
    contribution,
    length(test_rows),
    length(model_features)
  )
  contribution_names <- colnames(contribution)
  if (
    !is.null(contribution_names) &&
      !identical(
        as.character(contribution_names[seq_along(model_features)]),
        model_features
      )
  ) {
    stop("TreeSHAP contribution columns do not match the frozen feature order.")
  }

  shap_matrix <- contribution[
    ,
    seq_along(model_features),
    drop = FALSE
  ]
  colnames(shap_matrix) <- model_features
  base_value <- contribution[, length(model_features) + 1L]

  shap_additivity_residual <- abs(
    rowSums(shap_matrix) + base_value - booster_prediction
  )
  shap_additivity_allowed <- SHAP_ADDITIVITY_TOLERANCE +
    SHAP_ADDITIVITY_RELATIVE_TOLERANCE * pmax(
      1,
      abs(booster_prediction)
    )
  maximum_shap_additivity_error <- max(shap_additivity_residual)
  maximum_shap_scaled_error <- max(
    shap_additivity_residual / shap_additivity_allowed
  )
  if (
    !is.finite(maximum_shap_scaled_error) ||
      maximum_shap_scaled_error > 1
  ) {
    failing_row <- which.max(
      shap_additivity_residual / shap_additivity_allowed
    )
    stop(
      "TreeSHAP additivity failed for ",
      job_scheme,
      " fold ",
      job_fold_id,
      ". Maximum absolute error = ",
      signif(maximum_shap_additivity_error, 8),
      "; largest scaled error = ",
      signif(maximum_shap_scaled_error, 8),
      "; allowed error at worst row = ",
      signif(shap_additivity_allowed[failing_row], 8),
      " (absolute tolerance + relative tolerance x prediction magnitude)"
    )
  }

  metadata <- cohort[test_rows, .(
    obs_id,
    source_row_id,
    env_id,
    site_id,
    season_year,
    env_zone,
    observation_date,
    longitude,
    latitude,
    spatial_block,
    joint_temporal_block
  )]
  metadata[, BBCH_stage := as.numeric(cohort[[TARGET]][test_rows])]
  metadata[, BBCH_phase := phase_from_bbch(BBCH_stage)]
  metadata[, response := prediction]
  metadata[, error := response - BBCH_stage]
  metadata[, base_value := base_value]
  metadata[, scheme := job_scheme]
  metadata[, split_id := job_split_id]
  metadata[, fold_id := job_fold_id]
  metadata[
    ,
    bootstrap_cluster := scheme_bootstrap_cluster(
      scheme,
      site_id,
      season_year,
      env_zone,
      fold_id,
      env_id
    )
  ]

  x_table <- as.data.table(test_matrix[, configured_feature_index, drop = FALSE])
  setnames(x_table, paste0("x__", features))
  shap_table <- as.data.table(
    shap_matrix[, configured_feature_index, drop = FALSE]
  )
  setnames(shap_table, paste0("shap__", features))
  row_shap <- cbind(metadata, x_table, shap_table)

  support <- feature_support_audit(
    train_matrix = train_matrix,
    test_matrix = test_matrix,
    scheme = job_scheme,
    split_id = job_split_id,
    fold_id = job_fold_id,
    features = model_features
  )

  interaction_long <- NULL
  interaction_additivity_error <- NA_real_
  n_interaction_rows <- 0L

  if (RUN_INTERACTIONS) {
    interaction_rows_local <- match(
      selected_interaction_ids,
      metadata$obs_id
    )
    if (anyNA(interaction_rows_local)) {
      stop("Interaction sample could not be aligned to fold ", job_fold_id)
    }

    if (length(interaction_rows_local)) {
      interaction_matrix <- test_matrix[
        interaction_rows_local,
        ,
        drop = FALSE
      ]
      dinteraction <- xgb.DMatrix(
        data = interaction_matrix,
        missing = NA_real_
      )
      interaction_array <- predict(
        booster,
        dinteraction,
        predinteraction = TRUE,
        approxcontrib = FALSE,
        strict_shape = FALSE
      )
      interaction_array <- normalise_interaction_array(
        interaction_result = interaction_array,
        n_rows = length(interaction_rows_local),
        n_features = length(model_features)
      )

      interaction_prediction <- booster_prediction[interaction_rows_local]
      interaction_additivity_residual <- abs(
        apply(interaction_array, 1L, sum) - interaction_prediction
      )
      interaction_additivity_allowed <-
        INTERACTION_ADDITIVITY_TOLERANCE +
        INTERACTION_ADDITIVITY_RELATIVE_TOLERANCE * pmax(
          1,
          abs(interaction_prediction)
        )
      interaction_additivity_error <- max(
        interaction_additivity_residual
      )
      maximum_interaction_scaled_error <- max(
        interaction_additivity_residual / interaction_additivity_allowed
      )
      if (
        !is.finite(maximum_interaction_scaled_error) ||
          maximum_interaction_scaled_error > 1
      ) {
        failing_interaction_row <- which.max(
          interaction_additivity_residual / interaction_additivity_allowed
        )
        stop(
          "TreeSHAP interaction additivity failed for ",
          job_scheme,
          " fold ",
          job_fold_id,
          ". Maximum absolute error = ",
          signif(interaction_additivity_error, 8),
          "; largest scaled error = ",
          signif(maximum_interaction_scaled_error, 8),
          "; allowed error at worst row = ",
          signif(
            interaction_additivity_allowed[failing_interaction_row],
            8
          ),
          " (absolute tolerance + relative tolerance x prediction magnitude)"
        )
      }

      interaction_metadata <- metadata[interaction_rows_local]
      interaction_long <- interaction_array_to_long(
        interaction_array = interaction_array,
        feature_matrix = interaction_matrix,
        metadata = interaction_metadata,
        features = model_features
      )
      interaction_long[, scheme := job_scheme]
      interaction_long[, split_id := job_split_id]
      interaction_long[, fold_id := job_fold_id]
      n_interaction_rows <- length(interaction_rows_local)
    }
  }

  result <- list(
    checkpoint_version = CHECKPOINT_VERSION,
    input_md5 = input_md5,
    feature_set = feature_set,
    scheme = job_scheme,
    split_id = job_split_id,
    fold_id = job_fold_id,
    split_signature = split_signature,
    interaction_signature = interaction_signature,
    tuned_parameter_signature = tuned_parameter_signature,
    tuned_parameter_values = tuned_parameter_values,
    reconstruction_engine = "mlr3learners::regr.xgboost",
    model_features = model_features,
    xgboost_version = as.character(utils::packageVersion("xgboost")),
    n_train = length(train_rows),
    n_test = length(test_rows),
    n_excluded = length(excluded_rows),
    n_interaction_rows = n_interaction_rows,
    maximum_prediction_difference = maximum_prediction_difference,
    maximum_wrapper_booster_difference =
      maximum_wrapper_booster_difference,
    maximum_shap_additivity_error =
      maximum_shap_additivity_error,
    interaction_additivity_error =
      interaction_additivity_error,
    row_shap = row_shap,
    interaction_long = interaction_long,
    support = support
  )

  atomic_save_rds(result, checkpoint_file, compress = FALSE)
  result
}

# Standardise older checkpoint schemas before aggregation. Early SHAP script
# versions used `row_id` and `truth`; later versions use `obs_id` and
# `BBCH_stage`. The scientific contents are the same, but post-processing must
# work with a single schema.
rename_checkpoint_column <- function(data, target, aliases) {
  if (target %in% names(data)) return(invisible(NULL))
  source <- intersect(aliases, names(data))
  if (length(source)) setnames(data, source[[1L]], target)
  invisible(NULL)
}

normalise_checkpoint_table <- function(
  data,
  job,
  label,
  require_predictions = FALSE
) {
  if (is.null(data)) return(NULL)

  table <- as.data.table(copy(data))
  rename_checkpoint_column(table, "obs_id", c("row_id", "row_ids"))
  rename_checkpoint_column(table, "BBCH_stage", c("truth", "BBCH"))
  rename_checkpoint_column(table, "response", c("prediction", "pred"))

  if (!"obs_id" %in% names(table)) {
    stop("Checkpoint ", label, " has no obs_id or legacy row_id column.")
  }
  table[, obs_id := as.integer(obs_id)]

  cohort_match <- match(table$obs_id, cohort$obs_id)
  if (anyNA(cohort_match)) {
    stop(
      "Checkpoint ", label,
      " contains observation IDs absent from the current benchmark cohort."
    )
  }

  cohort_metadata <- c(
    "source_row_id",
    "env_id",
    "site_id",
    "season_year",
    "env_zone",
    "observation_date",
    "longitude",
    "latitude",
    "spatial_block",
    "joint_temporal_block"
  )
  for (column in cohort_metadata) {
    if (!column %in% names(table)) {
      table[, (column) := cohort[[column]][cohort_match]]
    }
  }
  if (!"BBCH_stage" %in% names(table)) {
    table[, BBCH_stage := as.numeric(cohort[[TARGET]][cohort_match])]
  }
  table[, BBCH_stage := as.numeric(BBCH_stage)]
  if (!"BBCH_phase" %in% names(table)) {
    table[, BBCH_phase := phase_from_bbch(BBCH_stage)]
  }

  expected_scheme <- as.character(job$scheme[[1L]])
  expected_split_id <- as.integer(job$split_id[[1L]])
  expected_fold_id <- as.character(job$fold_id[[1L]])
  if (!"scheme" %in% names(table)) table[, scheme := expected_scheme]
  if (!"split_id" %in% names(table)) table[, split_id := expected_split_id]
  if (!"fold_id" %in% names(table)) table[, fold_id := expected_fold_id]
  if (
    any(as.character(table$scheme) != expected_scheme) ||
      any(as.integer(table$split_id) != expected_split_id) ||
      any(as.character(table$fold_id) != expected_fold_id)
  ) {
    stop("Checkpoint ", label, " does not match its requested outer fold.")
  }

  if (!"bootstrap_cluster" %in% names(table)) {
    table[, bootstrap_cluster := scheme_bootstrap_cluster(
      scheme,
      site_id,
      season_year,
      env_zone,
      fold_id,
      env_id
    )]
  }

  if (require_predictions) {
    required_shap_columns <- paste0("shap__", features)
    missing_shap_columns <- setdiff(required_shap_columns, names(table))
    if (length(missing_shap_columns)) {
      stop(
        "Checkpoint row-level SHAP table is missing: ",
        paste(missing_shap_columns, collapse = ", ")
      )
    }
    if (!"response" %in% names(table)) {
      stop("Checkpoint row-level SHAP table has no prediction column.")
    }
    table[, response := as.numeric(response)]
    if (!"error" %in% names(table)) {
      table[, error := response - BBCH_stage]
    }
    if (!"base_value" %in% names(table)) {
      stop("Checkpoint row-level SHAP table has no TreeSHAP base_value column.")
    }
  }

  table[]
}

load_checkpoint_only_job <- function(job) {
  job_scheme <- as.character(job$scheme[[1L]])
  job_split_id <- as.integer(job$split_id[[1L]])
  job_fold_id <- as.character(job$fold_id[[1L]])
  checkpoint_file <- file.path(
    CHECKPOINT_DIR,
    paste0(
      safe_id(feature_set),
      "__",
      safe_id(job_scheme),
      "__split_",
      sprintf("%03d", job_split_id),
      "__fold_",
      safe_id(job_fold_id),
      ".rds"
    )
  )
  if (!file.exists(checkpoint_file)) {
    stop(
      "Checkpoint-only mode will not refit missing work. Missing checkpoint: ",
      checkpoint_file
    )
  }

  cached <- tryCatch(readRDS(checkpoint_file), error = function(error) error)
  if (inherits(cached, "error")) {
    stop("Could not read checkpoint ", checkpoint_file, ": ",
         conditionMessage(cached))
  }
  if (
    !identical(as.character(cached$feature_set), feature_set) ||
      !identical(as.character(cached$scheme), job_scheme) ||
      !identical(as.integer(cached$split_id), job_split_id) ||
      !identical(as.character(cached$fold_id), job_fold_id)
  ) {
    stop("Checkpoint metadata does not match the requested fold: ",
         checkpoint_file)
  }
  if (is.null(cached$row_shap)) {
    stop("Checkpoint contains no row-level SHAP values: ", checkpoint_file)
  }

  cached$row_shap <- normalise_checkpoint_table(
    cached$row_shap,
    job = job,
    label = checkpoint_file,
    require_predictions = TRUE
  )
  cached$interaction_long <- normalise_checkpoint_table(
    cached$interaction_long,
    job = job,
    label = checkpoint_file,
    require_predictions = FALSE
  )
  cached$checkpoint_reused <- TRUE
  cached
}

message("Feature set: ", feature_set)
message("Schemes: ", paste(schemes, collapse = ", "))

if (CHECKPOINT_MODE == "checkpoints_only") {
  message(
    "Checkpoint-only mode: loading ", nrow(jobs),
    " saved fold objects; no XGBoost refitting or TreeSHAP calculation will run."
  )
  shap_results <- lapply(seq_len(nrow(jobs)), function(i) {
    load_checkpoint_only_job(jobs[i])
  })
} else {
  message(
    "Starting ", nrow(jobs), " outer-fold SHAP jobs with ",
    workers, " worker(s)."
  )
  options(future.globals.maxSize = 8 * 1024^3)
  future::plan(future::multisession, workers = workers)
  shap_results <- future.apply::future_lapply(
    seq_len(nrow(jobs)),
    function(i) run_shap_job(jobs[i]),
    future.packages = c(
      "data.table",
      "R6",
      "mlr3",
      "mlr3learners",
      "mlr3spatiotempcv",
      "xgboost"
    ),
    future.seed = TRUE,
    future.scheduling = 1
  )
  future::plan(future::sequential)
}

# ------------------------------------------------------------------------------
# 8. COMBINE ROW-LEVEL SHAP AND TECHNICAL DIAGNOSTICS
# ------------------------------------------------------------------------------

row_shap <- rbindlist(
  lapply(shap_results, `[[`, "row_shap"),
  use.names = TRUE
)
setorder(row_shap, scheme, fold_id, obs_id)

if (
  row_shap[
    ,
    anyDuplicated(obs_id),
    by = scheme
  ][, any(V1 > 0L)]
) {
  stop("An observation has duplicate OOF SHAP values within a scheme.")
}

expected_rows <- nrow(benchmark_oof)
if (nrow(row_shap) != expected_rows) {
  stop(
    "OOF SHAP row count is incomplete. Expected ",
    expected_rows,
    ", found ",
    nrow(row_shap)
  )
}
for (scheme_id in schemes) {
  expected_ids <- sort(benchmark_oof[scheme == scheme_id, row_ids])
  observed_ids <- sort(row_shap[scheme == scheme_id, obs_id])
  if (!identical(as.integer(observed_ids), as.integer(expected_ids))) {
    stop("OOF SHAP observation coverage differs for scheme: ", scheme_id)
  }
}

atomic_save_rds(
  row_shap,
  file.path(OUTPUT_DIR, "oof_shap_values_wide.rds"),
  compress = FALSE
)

technical_diagnostics <- rbindlist(lapply(shap_results, function(x) {
  data.table(
    feature_set = x$feature_set,
    scheme = x$scheme,
    split_id = x$split_id,
    fold_id = x$fold_id,
    benchmark_job_id = x$benchmark_job_id,
    n_train = x$n_train,
    n_test = x$n_test,
    n_excluded = x$n_excluded,
    n_interaction_rows = x$n_interaction_rows,
    reconstruction_engine = x$reconstruction_engine,
    maximum_prediction_difference =
      x$maximum_prediction_difference,
    maximum_wrapper_booster_difference =
      x$maximum_wrapper_booster_difference,
    maximum_shap_additivity_error =
      x$maximum_shap_additivity_error,
    interaction_additivity_error =
      x$interaction_additivity_error,
    xgboost_version = x$xgboost_version
  )
}))
write_table(technical_diagnostics, "shap_technical_diagnostics.csv")

support_audit <- rbindlist(
  lapply(shap_results, `[[`, "support"),
  use.names = TRUE
)
write_table(support_audit, "held_out_feature_support_audit.csv")

base_value_summary <- row_shap[, .(
  n = .N,
  mean_base_value = mean(base_value),
  sd_base_value = sd(base_value),
  min_base_value = min(base_value),
  max_base_value = max(base_value)
), by = .(scheme, split_id, fold_id)]
write_table(base_value_summary, "fold_base_value_summary.csv")

# ------------------------------------------------------------------------------
# 9. FEATURE-LEVEL MACRO SHAP
# ------------------------------------------------------------------------------

shap_columns <- paste0("shap__", features)

shap_long <- melt(
  row_shap,
  id.vars = c(
    "obs_id",
    "source_row_id",
    "env_id",
    "site_id",
    "season_year",
    "env_zone",
    "observation_date",
    "BBCH_stage",
    "BBCH_phase",
    "response",
    "error",
    "bootstrap_cluster",
    "scheme",
    "split_id",
    "fold_id"
  ),
  measure.vars = shap_columns,
  variable.name = "feature",
  value.name = "shap_value",
  variable.factor = FALSE
)
shap_long[
  ,
  feature := sub("^shap__", "", feature)
]

# Primary statistical unit: one environment x exact observed BBCH stage.
feature_unit <- shap_long[, .(
  mean_abs_shap = mean(abs(shap_value)),
  mean_shap = mean(shap_value),
  sd_shap = sd(shap_value),
  n_observations = .N
), by = .(
  scheme,
  split_id,
  fold_id,
  env_id,
  site_id,
  season_year,
  env_zone,
  bootstrap_cluster,
  BBCH_stage,
  BBCH_phase,
  feature
)]

global_importance_macro <- feature_unit[, .(
  n_environment_stage_units = .N,
  mean_abs_shap = mean(mean_abs_shap),
  mean_shap = mean(mean_shap),
  sd_unit_abs_shap = sd(mean_abs_shap)
), by = .(scheme, feature)]
global_importance_macro[
  ,
  importance_rank := frank(
    -mean_abs_shap,
    ties.method = "average"
  ),
  by = scheme
]
global_importance_macro[
  ,
  importance_share := mean_abs_shap / sum(mean_abs_shap),
  by = scheme
]
setorder(global_importance_macro, scheme, importance_rank)

bootstrap_feature_importance <- rbindlist(lapply(
  seq_len(nrow(global_importance_macro)),
  function(i) {
    group_row <- global_importance_macro[i]
    unit_subset <- feature_unit[
      scheme == group_row$scheme &
        feature == group_row$feature
    ]
    result <- cluster_bootstrap_mean(
      unit_table = unit_subset,
      value_column = "mean_abs_shap",
      cluster_column = "bootstrap_cluster",
      replicates = BOOTSTRAP_REPLICATES,
      seed = SEED + 20000L + i
    )
    result[, scheme := group_row$scheme]
    result[, feature := group_row$feature]
    result
  }
))

global_importance_macro <- merge(
  global_importance_macro,
  bootstrap_feature_importance[
    ,
    .(
      scheme,
      feature,
      ci_low,
      ci_high,
      n_bootstrap_clusters = n_clusters
    )
  ],
  by = c("scheme", "feature"),
  all.x = TRUE
)
setorder(global_importance_macro, scheme, importance_rank)

write_table(
  global_importance_macro,
  "global_shap_importance_macro_environment_stage.csv"
)
write_table(
  feature_unit,
  "shap_environment_stage_units.csv"
)

# Supplementary row-weighted importance.
global_importance_row_weighted <- shap_long[, .(
  n_rows = .N,
  mean_abs_shap = mean(abs(shap_value)),
  mean_shap = mean(shap_value)
), by = .(scheme, feature)]
global_importance_row_weighted[
  ,
  importance_rank := frank(
    -mean_abs_shap,
    ties.method = "average"
  ),
  by = scheme
]
setorder(
  global_importance_row_weighted,
  scheme,
  importance_rank
)
write_table(
  global_importance_row_weighted,
  "global_shap_importance_row_weighted_supplement.csv"
)

# Stage- and phase-specific attribution.
stage_importance <- feature_unit[, .(
  n_environment_stage_units = .N,
  mean_abs_shap = mean(mean_abs_shap),
  mean_shap = mean(mean_shap),
  sd_unit_abs_shap = sd(mean_abs_shap)
), by = .(scheme, BBCH_stage, feature)]
stage_importance[
  ,
  stage_rank := frank(
    -mean_abs_shap,
    ties.method = "average"
  ),
  by = .(scheme, BBCH_stage)
]
stage_importance[
  ,
  importance_share := mean_abs_shap / sum(mean_abs_shap),
  by = .(scheme, BBCH_stage)
]
setorder(stage_importance, scheme, BBCH_stage, stage_rank)
write_table(stage_importance, "shap_importance_by_exact_bbch.csv")

phase_importance <- feature_unit[, .(
  n_environment_stage_units = .N,
  mean_abs_shap = mean(mean_abs_shap),
  mean_shap = mean(mean_shap),
  sd_unit_abs_shap = sd(mean_abs_shap)
), by = .(scheme, BBCH_phase, feature)]
phase_importance[
  ,
  phase_rank := frank(
    -mean_abs_shap,
    ties.method = "average"
  ),
  by = .(scheme, BBCH_phase)
]
phase_importance[
  ,
  importance_share := mean_abs_shap / sum(mean_abs_shap),
  by = .(scheme, BBCH_phase)
]
setorder(phase_importance, scheme, BBCH_phase, phase_rank)
write_table(phase_importance, "shap_importance_by_bbch_phase.csv")

# The generalisation benchmark identified BBCH 75 as the weakest stage. This
# table describes which features push predictions earlier/later at that stage;
# it is descriptive attribution and is not a causal stress-effect estimate.
focal_stage_unit <- feature_unit[BBCH_stage %in% FOCAL_BBCH_STAGES]
focal_stage_summary <- focal_stage_unit[, .(
  n_environment_stage_units = .N,
  mean_abs_shap = mean(mean_abs_shap),
  mean_shap = mean(mean_shap),
  median_shap = median(mean_shap),
  q25_shap = as.numeric(quantile(mean_shap, 0.25)),
  q75_shap = as.numeric(quantile(mean_shap, 0.75))
), by = .(scheme, BBCH_stage, feature)]
focal_stage_summary[
  ,
  importance_rank := frank(-mean_abs_shap, ties.method = "average"),
  by = .(scheme, BBCH_stage)
]
setorder(focal_stage_summary, scheme, BBCH_stage, importance_rank)
write_table(focal_stage_summary, "shap_focal_bbch_stage_summary.csv")

focal_error_unit <- shap_long[
  BBCH_stage %in% FOCAL_BBCH_STAGES,
  .(
    mean_shap = mean(shap_value),
    mean_prediction_error = mean(error)
  ),
  by = .(
    scheme,
    split_id,
    fold_id,
    env_id,
    BBCH_stage,
    feature
  )
]
focal_error_association <- focal_error_unit[, .(
  n_environment_stage_units = .N,
  spearman_shap_error = cor_or_na(
    mean_shap,
    mean_prediction_error,
    method = "spearman"
  ),
  pearson_shap_error = cor_or_na(
    mean_shap,
    mean_prediction_error,
    method = "pearson"
  )
), by = .(scheme, BBCH_stage, feature)]
write_table(
  focal_error_association,
  "shap_focal_bbch_error_association_descriptive.csv"
)

# Fold stability.
fold_feature_importance <- feature_unit[, .(
  n_environment_stage_units = .N,
  mean_abs_shap = mean(mean_abs_shap),
  mean_shap = mean(mean_shap)
), by = .(scheme, split_id, fold_id, feature)]
fold_feature_importance[
  ,
  fold_rank := frank(
    -mean_abs_shap,
    ties.method = "average"
  ),
  by = .(scheme, split_id, fold_id)
]

feature_fold_stability <- fold_feature_importance[, .(
  n_folds = .N,
  mean_fold_importance = mean(mean_abs_shap),
  sd_fold_importance = if (.N > 1L) sd(mean_abs_shap) else NA_real_,
  min_fold_importance = min(mean_abs_shap),
  max_fold_importance = max(mean_abs_shap),
  mean_fold_rank = mean(fold_rank),
  sd_fold_rank = if (.N > 1L) sd(fold_rank) else NA_real_,
  best_fold_rank = min(fold_rank),
  worst_fold_rank = max(fold_rank)
), by = .(scheme, feature)]
setorder(
  feature_fold_stability,
  scheme,
  mean_fold_rank
)
write_table(fold_feature_importance, "shap_importance_by_fold.csv")
write_table(feature_fold_stability, "shap_fold_stability.csv")

# ------------------------------------------------------------------------------
# 10. BIOLOGICAL GROUP ATTRIBUTION
# ------------------------------------------------------------------------------

group_id_columns <- c(
  "obs_id",
  "env_id",
  "site_id",
  "season_year",
  "env_zone",
  "bootstrap_cluster",
  "BBCH_stage",
  "BBCH_phase",
  "scheme",
  "split_id",
  "fold_id"
)
missing_group_id_columns <- setdiff(group_id_columns, names(row_shap))
if (length(missing_group_id_columns)) {
  stop(
    "Cannot construct biological-group SHAP attribution. Missing row-level ",
    "identifier columns: ",
    paste(missing_group_id_columns, collapse = ", ")
  )
}

group_shap_wide <- copy(row_shap[, ..group_id_columns])

for (group_name in unique(feature_groups$biological_group)) {
  group_features <- feature_groups[
    biological_group == group_name,
    feature
  ]
  group_feature_columns <- paste0("shap__", group_features)
  missing_group_feature_columns <- setdiff(
    group_feature_columns,
    names(row_shap)
  )
  if (length(missing_group_feature_columns)) {
    stop(
      "Cannot construct biological group '", group_name,
      "'. Missing SHAP columns: ",
      paste(missing_group_feature_columns, collapse = ", ")
    )
  }

  group_shap_wide[
    ,
    (paste0("group__", group_name)) := rowSums(
      row_shap[, .SD, .SDcols = group_feature_columns]
    )
  ]
}

group_measure_columns <- grep(
  "^group__",
  names(group_shap_wide),
  value = TRUE
)
if (!length(group_measure_columns)) {
  stop("No biological-group SHAP columns were constructed.")
}
group_long <- melt(
  group_shap_wide,
  id.vars = group_id_columns,
  measure.vars = group_measure_columns,
  variable.name = "biological_group",
  value.name = "group_shap",
  variable.factor = FALSE
)
group_long[
  ,
  biological_group := sub(
    "^group__",
    "",
    biological_group
  )
]

group_unit <- group_long[, .(
  mean_abs_group_shap = mean(abs(group_shap)),
  mean_group_shap = mean(group_shap),
  n_observations = .N
), by = .(
  scheme,
  split_id,
  fold_id,
  env_id,
  site_id,
  season_year,
  env_zone,
  bootstrap_cluster,
  BBCH_stage,
  BBCH_phase,
  biological_group
)]

group_importance <- group_unit[, .(
  n_environment_stage_units = .N,
  mean_abs_group_shap = mean(mean_abs_group_shap),
  mean_group_shap = mean(mean_group_shap)
), by = .(scheme, biological_group)]
group_importance[
  ,
  importance_rank := frank(
    -mean_abs_group_shap,
    ties.method = "average"
  ),
  by = scheme
]
group_importance[
  ,
  importance_share := mean_abs_group_shap / sum(mean_abs_group_shap),
  by = scheme
]

bootstrap_group_importance <- rbindlist(lapply(
  seq_len(nrow(group_importance)),
  function(i) {
    group_row <- group_importance[i]
    unit_subset <- group_unit[
      scheme == group_row$scheme &
        biological_group == group_row$biological_group
    ]
    result <- cluster_bootstrap_mean(
      unit_table = unit_subset,
      value_column = "mean_abs_group_shap",
      cluster_column = "bootstrap_cluster",
      replicates = BOOTSTRAP_REPLICATES,
      seed = SEED + 25000L + i
    )
    result[, scheme := group_row$scheme]
    result[, biological_group := group_row$biological_group]
    result
  }
))
group_importance <- merge(
  group_importance,
  bootstrap_group_importance[, .(
    scheme,
    biological_group,
    ci_low,
    ci_high,
    n_bootstrap_clusters = n_clusters
  )],
  by = c("scheme", "biological_group"),
  all.x = TRUE
)
setorder(group_importance, scheme, importance_rank)

group_stage_importance <- group_unit[, .(
  n_environment_stage_units = .N,
  mean_abs_group_shap = mean(mean_abs_group_shap),
  mean_group_shap = mean(mean_group_shap)
), by = .(scheme, BBCH_stage, biological_group)]
group_stage_importance[
  ,
  stage_rank := frank(
    -mean_abs_group_shap,
    ties.method = "average"
  ),
  by = .(scheme, BBCH_stage)
]
group_stage_importance[
  ,
  importance_share := mean_abs_group_shap / sum(mean_abs_group_shap),
  by = .(scheme, BBCH_stage)
]

write_table(
  group_importance,
  "biological_group_shap_importance.csv"
)
write_table(
  group_stage_importance,
  "biological_group_shap_by_exact_bbch.csv"
)

# ------------------------------------------------------------------------------
# 11. CROSS-SCHEME ATTRIBUTION STABILITY
# ------------------------------------------------------------------------------

scheme_pairs <- if (length(schemes) >= 2L) {
  combn(sort(schemes), 2L, simplify = FALSE)
} else {
  list()
}

global_rank_stability <- rbindlist(lapply(
  seq_along(scheme_pairs),
  function(i) {
    pair <- scheme_pairs[[i]]
    left <- global_importance_macro[
      scheme == pair[1L],
      .(
        feature,
        importance_left = mean_abs_shap,
        rank_left = importance_rank
      )
    ]
    right <- global_importance_macro[
      scheme == pair[2L],
      .(
        feature,
        importance_right = mean_abs_shap,
        rank_right = importance_rank
      )
    ]
    paired <- merge(left, right, by = "feature")

    data.table(
      scheme_left = pair[1L],
      scheme_right = pair[2L],
      n_features = nrow(paired),
      spearman_importance = cor_or_na(
        paired$importance_left,
        paired$importance_right,
        method = "spearman"
      ),
      spearman_rank = cor_or_na(
        paired$rank_left,
        paired$rank_right,
        method = "spearman"
      ),
      top_3_overlap = length(intersect(
        head(left[order(rank_left), feature], 3L),
        head(right[order(rank_right), feature], 3L)
      ))
    )
  }
), fill = TRUE)

if (nrow(global_rank_stability)) {
  write_table(
    global_rank_stability,
    "shap_global_cross_scheme_stability.csv"
  )
}

observation_attribution_stability <- rbindlist(lapply(
  seq_along(scheme_pairs),
  function(i) {
    pair <- scheme_pairs[[i]]
    left <- shap_long[
      scheme == pair[1L],
      .(obs_id, feature, shap_left = shap_value)
    ]
    right <- shap_long[
      scheme == pair[2L],
      .(obs_id, feature, shap_right = shap_value)
    ]
    paired <- merge(
      left,
      right,
      by = c("obs_id", "feature"),
      all = FALSE
    )

    paired[, .(
      n_observations = .N,
      pearson_shap = cor_or_na(
        shap_left,
        shap_right,
        method = "pearson"
      ),
      spearman_shap = cor_or_na(
        shap_left,
        shap_right,
        method = "spearman"
      ),
      mean_absolute_difference =
        mean(abs(shap_left - shap_right)),
      sign_agreement = mean(
        sign(shap_left) == sign(shap_right)
      )
    ), by = feature][
      ,
      `:=`(
        scheme_left = pair[1L],
        scheme_right = pair[2L]
      )
    ]
  }
), fill = TRUE)

if (nrow(observation_attribution_stability)) {
  setcolorder(
    observation_attribution_stability,
    c(
      "scheme_left",
      "scheme_right",
      "feature",
      setdiff(
        names(observation_attribution_stability),
        c("scheme_left", "scheme_right", "feature")
      )
    )
  )
  write_table(
    observation_attribution_stability,
    "shap_observation_cross_scheme_stability.csv"
  )
}

# ------------------------------------------------------------------------------
# 12. PREDICTOR CORRELATION AUDIT
# ------------------------------------------------------------------------------

feature_matrix <- as.matrix(cohort[, ..features])
spearman_correlation <- cor(
  feature_matrix,
  method = "spearman",
  use = "pairwise.complete.obs"
)

correlation_matrix_dt <- as.data.table(
  spearman_correlation,
  keep.rownames = "feature"
)
write_table(
  correlation_matrix_dt,
  "predictor_spearman_correlation_matrix.csv"
)

if (length(features) >= 2L) {
  correlation_pairs_index <- which(
    upper.tri(spearman_correlation),
    arr.ind = TRUE
  )
  correlation_pairs <- data.table(
    feature_1 = rownames(spearman_correlation)[
      correlation_pairs_index[, 1L]
    ],
    feature_2 = colnames(spearman_correlation)[
      correlation_pairs_index[, 2L]
    ],
    spearman_rho = spearman_correlation[
      correlation_pairs_index
    ]
  )
  correlation_pairs[
    ,
    high_correlation :=
      abs(spearman_rho) >= CORRELATION_WARNING_THRESHOLD
  ]
  correlation_pairs[
    ,
    abs_spearman_rho := abs(spearman_rho)
  ]
  setorder(correlation_pairs, -abs_spearman_rho)
  write_table(
    correlation_pairs,
    "predictor_spearman_correlation_pairs.csv"
  )
}

# ------------------------------------------------------------------------------
# 13. EXACT TREESHAP INTERACTIONS
# ------------------------------------------------------------------------------

interaction_all <- rbindlist(
  lapply(shap_results, `[[`, "interaction_long"),
  use.names = TRUE,
  fill = TRUE
)

if (RUN_INTERACTIONS && nrow(interaction_all)) {
  atomic_save_rds(
    interaction_all,
    file.path(OUTPUT_DIR, "oof_shap_interactions_sample.rds"),
    compress = FALSE
  )

  interaction_unit <- interaction_all[, .(
    mean_abs_interaction = mean(abs(interaction_value)),
    mean_interaction = mean(interaction_value),
    n_sampled_observations = .N
  ), by = .(
    scheme,
    split_id,
    fold_id,
    env_id,
    site_id,
    season_year,
    env_zone,
    BBCH_stage,
    BBCH_phase,
    pair_id,
    feature_1,
    feature_2
  )]
  interaction_unit[
    ,
    bootstrap_cluster := scheme_bootstrap_cluster(
      scheme,
      site_id,
      season_year,
      env_zone,
      fold_id,
      env_id
    )
  ]

  interaction_importance <- interaction_unit[, .(
    n_environment_stage_units = .N,
    mean_abs_interaction = mean(mean_abs_interaction),
    mean_interaction = mean(mean_interaction),
    sd_unit_abs_interaction = sd(mean_abs_interaction)
  ), by = .(
    scheme,
    pair_id,
    feature_1,
    feature_2
  )]
  interaction_importance[
    ,
    interaction_rank := frank(
      -mean_abs_interaction,
      ties.method = "average"
    ),
    by = scheme
  ]
  interaction_importance[
    ,
    interaction_share :=
      mean_abs_interaction / sum(mean_abs_interaction),
    by = scheme
  ]

  interaction_importance <- merge(
    interaction_importance,
    predeclared_interactions[, .(pair_id, hypothesis)],
    by = "pair_id",
    all.x = TRUE
  )
  interaction_importance[
    ,
    interaction_status := fifelse(
      !is.na(hypothesis),
      "pre_specified",
      "exploratory"
    )
  ]

  bootstrap_interactions <- rbindlist(lapply(
    seq_len(nrow(interaction_importance)),
    function(i) {
      group_row <- interaction_importance[i]
      unit_subset <- interaction_unit[
        scheme == group_row$scheme &
          pair_id == group_row$pair_id
      ]
      result <- cluster_bootstrap_mean(
        unit_table = unit_subset,
        value_column = "mean_abs_interaction",
        cluster_column = "bootstrap_cluster",
        replicates = BOOTSTRAP_REPLICATES,
        seed = SEED + 30000L + i
      )
      result[, scheme := group_row$scheme]
      result[, pair_id := group_row$pair_id]
      result
    }
  ))

  interaction_importance <- merge(
    interaction_importance,
    bootstrap_interactions[
      ,
      .(
        scheme,
        pair_id,
        ci_low,
        ci_high,
        n_bootstrap_clusters = n_clusters
      )
    ],
    by = c("scheme", "pair_id"),
    all.x = TRUE
  )
  setorder(
    interaction_importance,
    scheme,
    interaction_rank
  )

  interaction_stage <- interaction_unit[, .(
    n_environment_stage_units = .N,
    mean_abs_interaction = mean(mean_abs_interaction),
    mean_interaction = mean(mean_interaction)
  ), by = .(
    scheme,
    BBCH_stage,
    pair_id,
    feature_1,
    feature_2
  )]
  interaction_stage[
    ,
    stage_rank := frank(
      -mean_abs_interaction,
      ties.method = "average"
    ),
    by = .(scheme, BBCH_stage)
  ]

  interaction_fold <- interaction_unit[, .(
    n_environment_stage_units = .N,
    mean_abs_interaction = mean(mean_abs_interaction),
    mean_interaction = mean(mean_interaction)
  ), by = .(
    scheme,
    split_id,
    fold_id,
    pair_id,
    feature_1,
    feature_2
  )]
  interaction_fold[
    ,
    fold_rank := frank(
      -mean_abs_interaction,
      ties.method = "average"
    ),
    by = .(scheme, split_id, fold_id)
  ]

  interaction_fold_stability <- interaction_fold[, .(
    n_folds = .N,
    mean_fold_importance = mean(mean_abs_interaction),
    sd_fold_importance = if (.N > 1L) {
      sd(mean_abs_interaction)
    } else {
      NA_real_
    },
    mean_fold_rank = mean(fold_rank),
    sd_fold_rank = if (.N > 1L) sd(fold_rank) else NA_real_,
    best_fold_rank = min(fold_rank),
    worst_fold_rank = max(fold_rank)
  ), by = .(
    scheme,
    pair_id,
    feature_1,
    feature_2
  )]

  write_table(
    interaction_importance,
    "shap_interaction_importance_all_pairs.csv"
  )
  write_table(
    interaction_importance[
      interaction_status == "pre_specified"
    ],
    "shap_interaction_hypotheses_pre_specified.csv"
  )
  write_table(
    interaction_stage,
    "shap_interactions_by_exact_bbch.csv"
  )
  write_table(
    interaction_fold,
    "shap_interactions_by_fold.csv"
  )
  write_table(
    interaction_fold_stability,
    "shap_interaction_fold_stability.csv"
  )
}

# ------------------------------------------------------------------------------
# 14. PUBLICATION-ORIENTED PLOTS
# ------------------------------------------------------------------------------

feature_label_lookup <- c(
  cumGDD_from_sow_Soltani = "Thermal time (GDD)",
  cumDTU_beta = "Cumulative beta thermal time",
  sowing_DOY = "Sowing day of year",
  photoperiod_civil = "Civil photoperiod",
  Daylength = "Day length",
  cumRR = "Cumulative precipitation",
  cumRad = "Cumulative radiation",
  cumCDD = "Cumulative cold\ndegree-days",
  cumHDD = "Cumulative heat\ndegree-days",
  paw_fc33_0_28_mm = "Plant-available water\ncapacity (0–28 cm)",
  paw_fc33_0_100_mm = "Plant-available water\ncapacity (0–100 cm)"
)
missing_feature_labels <- setdiff(features, names(feature_label_lookup))
if (length(missing_feature_labels)) {
  feature_label_lookup[missing_feature_labels] <- gsub(
    "_",
    " ",
    missing_feature_labels,
    fixed = TRUE
  )
}
feature_short_label_lookup <- c(
  cumGDD_from_sow_Soltani = "GDD",
  cumDTU_beta = "Beta thermal time",
  sowing_DOY = "Sowing DOY",
  photoperiod_civil = "Photoperiod",
  Daylength = "Day length",
  cumRR = "Cum. precipitation",
  cumRad = "Cum. radiation",
  cumCDD = "Cum. cold DD",
  cumHDD = "Cum. heat DD",
  paw_fc33_0_28_mm = "PAW capacity (0–28 cm)",
  paw_fc33_0_100_mm = "PAW capacity (0–100 cm)"
)
missing_short_labels <- setdiff(features, names(feature_short_label_lookup))
if (length(missing_short_labels)) {
  feature_short_label_lookup[missing_short_labels] <- gsub(
    "_",
    " ",
    missing_short_labels,
    fixed = TRUE
  )
}
group_label_lookup <- c(
  developmental_clock = "Developmental clock",
  sowing_context = "Sowing context",
  photoperiod = "Photoperiod",
  water_supply = "Water supply",
  radiation = "Radiation",
  cold_exposure = "Cold exposure",
  heat_exposure = "Heat exposure",
  soil_water_capacity = "Soil water capacity"
)
missing_group_labels <- setdiff(
  unique(feature_groups$biological_group),
  names(group_label_lookup)
)
if (length(missing_group_labels)) {
  group_label_lookup[missing_group_labels] <- gsub(
    "_",
    " ",
    missing_group_labels,
    fixed = TRUE
  )
}
scheme_label_lookup <- c(
  spatial_block = "Spatial block",
  leave_one_year_out = "Leave-one-year-out",
  expanding_origin = "Expanding origin",
  environmental_zone = "Environmental zone",
  joint_space_time = "Joint space–time"
)
phase_label_lookup <- c(
  establishment = "Establishment",
  tillering = "Tillering",
  stem_elongation = "Stem elongation",
  heading_flowering = "Heading/flowering",
  grain_development_ripening = "Grain development/ripening",
  maturity_post_maturity = "Maturity/post-maturity"
)

scheme_display_levels <- unname(scheme_label_lookup[schemes])
paper_facet_columns <- min(2L, length(schemes))
paper_facet_rows <- ceiling(length(schemes) / paper_facet_columns)

add_scheme_display <- function(table) {
  answer <- copy(table)
  answer[, scheme_display := factor(
    unname(scheme_label_lookup[scheme]),
    levels = scheme_display_levels
  )]
  answer
}

theme_pheno <- function() {
  theme_minimal(base_size = 10.5) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9.5, margin = margin(b = 8)),
      strip.text = element_text(face = "bold", size = 10),
      axis.title = element_text(size = 10),
      axis.text = element_text(color = "black"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8.5),
      plot.margin = margin(8, 12, 8, 8)
    )
}

save_publication_plot <- function(plot, stem, width, height) {
  ggsave(
    file.path(PLOT_DIR, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 450,
    bg = "white",
    limitsize = FALSE
  )
  ggsave(
    file.path(PLOT_DIR, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    units = "in",
    device = "pdf",
    bg = "white",
    limitsize = FALSE
  )
}

importance_order <- global_importance_macro[
  ,
  .(mean_importance_rank = mean(importance_rank)),
  by = feature
][order(mean_importance_rank), feature]
feature_factor_levels <- rev(importance_order)

global_plot_data <- add_scheme_display(global_importance_macro)
global_plot_data[, feature_display := factor(
  feature,
  levels = feature_factor_levels,
  labels = unname(feature_label_lookup[feature_factor_levels])
)]

global_plot <- ggplot(
  global_plot_data,
  aes(x = mean_abs_shap, y = feature_display)
) +
  geom_col(fill = "#2C7FB8", width = 0.68) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    width = 0.16,
    linewidth = 0.45
  ) +
  facet_wrap(~scheme_display, ncol = paper_facet_columns) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = "Mean absolute OOF SHAP (BBCH units)",
    y = NULL,
    title = "Out-of-fold feature importance",
    subtitle = paste0(
      "Environment × observed-stage macro-average; 95% cluster-bootstrap ",
      "intervals"
    )
  ) +
  theme_pheno()
save_publication_plot(
  global_plot,
  "figure_shap_global_importance",
  width = 11.2,
  height = 5.8 * paper_facet_rows
)

group_order <- group_importance[
  ,
  .(mean_importance_rank = mean(importance_rank)),
  by = biological_group
][order(mean_importance_rank), biological_group]
group_factor_levels <- rev(group_order)
group_plot_data <- add_scheme_display(group_importance)
group_plot_data[, group_display := factor(
  biological_group,
  levels = group_factor_levels,
  labels = unname(group_label_lookup[group_factor_levels])
)]

group_plot <- ggplot(
  group_plot_data,
  aes(x = mean_abs_group_shap, y = group_display)
) +
  geom_col(fill = "#3A9D5D", width = 0.68) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    width = 0.16,
    linewidth = 0.45
  ) +
  facet_wrap(~scheme_display, ncol = paper_facet_columns) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = "Mean absolute grouped OOF SHAP (BBCH units)",
    y = NULL,
    title = "Biological process-group attribution",
    subtitle = paste0(
      "Signed feature SHAP values are summed within group before taking\n",
      "the absolute value"
    )
  ) +
  theme_pheno()
save_publication_plot(
  group_plot,
  "figure_shap_biological_groups",
  width = 11.2,
  height = 5.8 * paper_facet_rows
)

stage_plot_data <- add_scheme_display(stage_importance)
stage_plot_data[, feature_display := factor(
  feature,
  levels = feature_factor_levels,
  labels = unname(feature_label_lookup[feature_factor_levels])
)]
stage_breaks <- seq(0, 100, by = 10)

stage_plot_absolute <- ggplot(
  stage_plot_data,
  aes(x = BBCH_stage, y = feature_display, fill = mean_abs_shap)
) +
  geom_tile(width = 0.95, height = 0.92) +
  facet_wrap(~scheme_display, ncol = 1) +
  scale_x_continuous(
    breaks = stage_breaks,
    minor_breaks = seq(0, 100, by = 5),
    expand = expansion(add = 0.5)
  ) +
  scale_fill_viridis_c(
    option = "C",
    guide = guide_colorbar(
      title.position = "top",
      barwidth = 14,
      barheight = 0.55
    )
  ) +
  labs(
    x = "Observed BBCH stage",
    y = NULL,
    fill = "Mean |SHAP| (BBCH units)",
    title = "Stage-specific out-of-fold SHAP importance",
    subtitle = "Absolute attribution conditional on the observed stage"
  ) +
  theme_pheno() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 8.5)
  )
save_publication_plot(
  stage_plot_absolute,
  "figure_shap_stage_importance_absolute",
  width = 10.8,
  height = 4.2 + 2.4 * length(schemes)
)

stage_plot_share <- ggplot(
  stage_plot_data,
  aes(x = BBCH_stage, y = feature_display, fill = importance_share)
) +
  geom_tile(width = 0.95, height = 0.92) +
  facet_wrap(~scheme_display, ncol = 1) +
  scale_x_continuous(
    breaks = stage_breaks,
    minor_breaks = seq(0, 100, by = 5),
    expand = expansion(add = 0.5)
  ) +
  scale_fill_viridis_c(
    option = "B",
    labels = function(x) paste0(round(100 * x), "%"),
    guide = guide_colorbar(
      title.position = "top",
      barwidth = 14,
      barheight = 0.55
    )
  ) +
  labs(
    x = "Observed BBCH stage",
    y = NULL,
    fill = "Share of total mean |SHAP|",
    title = "Relative stage-specific attribution",
    subtitle = "Within-stage shares reveal changes that absolute magnitude can mask"
  ) +
  theme_pheno() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 8.5)
  )
save_publication_plot(
  stage_plot_share,
  "figure_shap_stage_importance_share",
  width = 10.8,
  height = 4.2 + 2.4 * length(schemes)
)

if (nrow(focal_stage_unit)) {
  focal_plot_data <- add_scheme_display(focal_stage_unit)
  focal_plot_data[, feature_display := factor(
    feature,
    levels = feature_factor_levels,
    labels = unname(feature_label_lookup[feature_factor_levels])
  )]
  focal_plot_data[, stage_display := paste("Observed BBCH", BBCH_stage)]

  focal_plot <- ggplot(
    focal_plot_data,
    aes(x = mean_shap, y = feature_display)
  ) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.45) +
    geom_boxplot(
      width = 0.62,
      outlier.shape = NA,
      fill = "#8FBBD9",
      color = "#184E77",
      linewidth = 0.45,
      orientation = "y"
    ) +
    facet_grid(stage_display ~ scheme_display, scales = "free_x") +
    labs(
      x = "Mean signed OOF SHAP per environment × stage unit (BBCH units)",
      y = NULL,
      title = "Attribution at the weakest stage (BBCH 75)",
      subtitle = "Boxes show the median and interquartile range across environments"
    ) +
    theme_pheno()
  save_publication_plot(
    focal_plot,
    "figure_shap_focal_bbch75_attribution",
    width = 11.2,
    height = 6.0
  )
}

# Balanced SHAP summary plot with within-feature colour scaling.
plot_sample <- balanced_plot_sample(
  shap_table = row_shap,
  maximum_rows = SUMMARY_PLOT_MAX_ROWS_PER_SCHEME,
  seed = SEED + 40000L
)

summary_shap_long <- melt(
  plot_sample,
  id.vars = c("obs_id", "scheme", "BBCH_stage", "BBCH_phase"),
  measure.vars = shap_columns,
  variable.name = "feature",
  value.name = "shap_value",
  variable.factor = FALSE
)
summary_shap_long[, feature := sub("^shap__", "", feature)]

x_columns <- paste0("x__", features)
summary_x_long <- melt(
  plot_sample,
  id.vars = c("obs_id", "scheme", "BBCH_stage", "BBCH_phase"),
  measure.vars = x_columns,
  variable.name = "feature",
  value.name = "feature_value",
  variable.factor = FALSE
)
summary_x_long[, feature := sub("^x__", "", feature)]

summary_plot_data <- merge(
  summary_shap_long,
  summary_x_long,
  by = c(
    "obs_id",
    "scheme",
    "BBCH_stage",
    "BBCH_phase",
    "feature"
  )
)
summary_plot_data[
  ,
  feature_value_scaled := {
    limits <- quantile(feature_value, c(0.01, 0.99), na.rm = TRUE)
    if (!is.finite(diff(limits)) || diff(limits) == 0) {
      rep(0.5, .N)
    } else {
      pmin(1, pmax(0, (feature_value - limits[1L]) / diff(limits)))
    }
  },
  by = .(scheme, feature)
]
summary_plot_data <- add_scheme_display(summary_plot_data)
summary_plot_data[, feature_display := factor(
  feature,
  levels = feature_factor_levels,
  labels = unname(feature_label_lookup[feature_factor_levels])
)]

summary_plot <- ggplot(
  summary_plot_data,
  aes(x = shap_value, y = feature_display, color = feature_value_scaled)
) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_point(
    alpha = 0.48,
    size = 0.62,
    position = position_jitter(height = 0.15, width = 0)
  ) +
  facet_wrap(~scheme_display, ncol = paper_facet_columns) +
  scale_color_gradientn(
    colors = c("#2166AC", "#BDBDBD", "#B2182B"),
    values = c(0, 0.5, 1),
    breaks = c(0, 1),
    labels = c("Low", "High"),
    guide = guide_colorbar(
      title.position = "top",
      barwidth = 12,
      barheight = 0.5
    )
  ) +
  labs(
    x = "OOF SHAP value (BBCH units)",
    y = NULL,
    color = "Feature value (within-feature 1st–99th percentile)",
    title = "Out-of-fold TreeSHAP summary",
    subtitle = paste0(
      "Environment × observed-stage-balanced plotting sample; vertical ",
      "jitter displays point density"
    )
  ) +
  theme_pheno()
save_publication_plot(
  summary_plot,
  "figure_shap_summary",
  width = 11.2,
  height = 6.4 * paper_facet_rows
)

# Dependence plots for the union of the top features across schemes.
top_dependence_features <- global_importance_macro[
  importance_rank <= DEPENDENCE_TOP_N,
  unique(feature)
]

for (feature_name in top_dependence_features) {
  dependence_data <- plot_sample[, .(
    scheme,
    BBCH_phase,
    feature_value = get(paste0("x__", feature_name)),
    shap_value = get(paste0("shap__", feature_name))
  )]
  dependence_data <- add_scheme_display(dependence_data)
  dependence_data[, phase_display := factor(
    unname(phase_label_lookup[BBCH_phase]),
    levels = unname(phase_label_lookup)
  )]

  dependence_plot <- ggplot(
    dependence_data,
    aes(x = feature_value, y = shap_value, color = phase_display)
  ) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
    geom_point(alpha = 0.28, size = 0.58) +
    geom_smooth(
      aes(group = 1),
      method = "gam",
      formula = y ~ s(x, k = 8),
      se = TRUE,
      color = "black",
      fill = "grey75",
      linewidth = 0.7
    ) +
    facet_wrap(
      ~scheme_display,
      ncol = paper_facet_columns,
      scales = "free_x"
    ) +
    scale_color_brewer(palette = "Dark2", drop = FALSE) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(
      x = gsub("\n", " ", unname(feature_label_lookup[feature_name])),
      y = "OOF SHAP value (BBCH units)",
      color = "Observed phase",
      title = paste(
        "SHAP dependence:",
        gsub("\n", " ", unname(feature_label_lookup[feature_name]))
      ),
      subtitle = paste0(
        "The smooth is descriptive and must not be read as a causal ",
        "response curve"
      )
    ) +
    theme_pheno()

  save_publication_plot(
    dependence_plot,
    paste0("figure_shap_dependence__", safe_id(feature_name)),
    width = 10.8,
    height = 5.8 * paper_facet_rows
  )
}

# The conventional summary plot above deliberately displaces points along the
# categorical feature axis to display density. Interaction clues instead come
# from vertical SHAP dispersion at the same predictor value in a dependence
# scatter. For every leading feature and CV scheme, identify its strongest
# exact TreeSHAP partner and colour an unstacked dependence scatter by that
# partner's observed value. The colour scale is normalized within each panel so
# that predictors with different units can be combined in one publication plot.
if (
  RUN_INTERACTIONS &&
    exists("interaction_importance") &&
    nrow(interaction_importance)
) {
  interaction_partner_candidates <- rbindlist(list(
    interaction_importance[, .(
      scheme,
      feature = feature_1,
      partner_feature = feature_2,
      pair_id,
      pair_mean_abs_interaction = mean_abs_interaction,
      pair_interaction_rank = interaction_rank
    )],
    interaction_importance[, .(
      scheme,
      feature = feature_2,
      partner_feature = feature_1,
      pair_id,
      pair_mean_abs_interaction = mean_abs_interaction,
      pair_interaction_rank = interaction_rank
    )]
  ))
  setorder(
    interaction_partner_candidates,
    scheme,
    feature,
    pair_interaction_rank,
    -pair_mean_abs_interaction
  )
  strongest_interaction_partner <- interaction_partner_candidates[
    ,
    .SD[1L],
    by = .(scheme, feature)
  ]

  leading_dependence_features <- global_importance_macro[
    importance_rank <= DEPENDENCE_TOP_N,
    .(scheme, feature, feature_importance_rank = importance_rank)
  ]
  strongest_interaction_partner <- merge(
    strongest_interaction_partner,
    leading_dependence_features,
    by = c("scheme", "feature")
  )
  setorder(
    strongest_interaction_partner,
    scheme,
    feature_importance_rank
  )
  write_table(
    strongest_interaction_partner,
    "shap_strongest_interaction_partner_for_top_features.csv"
  )

  partner_value_long <- copy(summary_x_long)
  setnames(
    partner_value_long,
    c("feature", "feature_value"),
    c("partner_feature", "partner_value")
  )

  dependence_interaction_data <- merge(
    summary_plot_data,
    strongest_interaction_partner,
    by = c("scheme", "feature")
  )
  dependence_interaction_data <- merge(
    dependence_interaction_data,
    partner_value_long,
    by = c(
      "obs_id",
      "scheme",
      "BBCH_stage",
      "BBCH_phase",
      "partner_feature"
    )
  )
  dependence_interaction_data[
    ,
    partner_value_scaled := {
      limits <- quantile(partner_value, c(0.01, 0.99), na.rm = TRUE)
      if (!is.finite(diff(limits)) || diff(limits) == 0) {
        rep(0.5, .N)
      } else {
        pmin(
          1,
          pmax(0, (partner_value - limits[1L]) / diff(limits))
        )
      }
    },
    by = .(scheme, feature, partner_feature)
  ]
  dependence_interaction_data[, panel_display := paste0(
    gsub("\n", " ", unname(feature_label_lookup[feature])),
    "\ncolour: ",
    gsub("\n", " ", unname(feature_label_lookup[partner_feature]))
  )]

  for (scheme_name in schemes) {
    scheme_dependence_data <- dependence_interaction_data[
      scheme == scheme_name
    ]
    if (!nrow(scheme_dependence_data)) next

    panel_levels <- unique(
      scheme_dependence_data[
        order(feature_importance_rank),
        panel_display
      ]
    )
    scheme_dependence_data[, panel_display := factor(
      panel_display,
      levels = panel_levels
    )]

    interaction_coloured_dependence_plot <- ggplot(
      scheme_dependence_data,
      aes(
        x = feature_value,
        y = shap_value,
        color = partner_value_scaled
      )
    ) +
      geom_hline(yintercept = 0, color = "grey60", linewidth = 0.35) +
      geom_point(alpha = 0.30, size = 0.60) +
      facet_wrap(~panel_display, ncol = 2L, scales = "free_x") +
      scale_color_gradientn(
        colors = c("#2166AC", "#BDBDBD", "#B2182B"),
        values = c(0, 0.5, 1),
        breaks = c(0, 1),
        labels = c("Low", "High"),
        guide = guide_colorbar(
          title.position = "top",
          barwidth = 11,
          barheight = 0.5
        )
      ) +
      labs(
        x = "Observed feature value",
        y = "OOF SHAP value (BBCH units)",
        color = "Strongest partner value\n(within-panel 1st–99th percentile)",
        title = paste0(
          "SHAP dependence and strongest interaction partner — ",
          unname(scheme_label_lookup[scheme_name])
        ),
        subtitle = paste0(
          "Points are not stacked; vertical dispersion at similar predictor ",
          "values indicates possible interaction structure"
        )
      ) +
      theme_pheno()

    save_publication_plot(
      interaction_coloured_dependence_plot,
      paste0(
        "figure_shap_dependence_interaction__",
        safe_id(scheme_name)
      ),
      width = 11.6,
      height = 3.3 * ceiling(length(panel_levels) / 2L) + 1.4
    )
  }
}

if (
  RUN_INTERACTIONS &&
    exists("interaction_importance") &&
    nrow(interaction_importance)
) {
  interaction_plot_data <- interaction_importance[
    interaction_rank <= INTERACTION_TOP_N
  ]
  interaction_plot_data <- add_scheme_display(interaction_plot_data)
  interaction_plot_data[, pair_display := paste0(
    unname(feature_short_label_lookup[feature_1]),
    "\n× ",
    unname(feature_short_label_lookup[feature_2])
  )]
  interaction_plot_data[, pair_display := factor(
    pair_display,
    levels = rev(unique(interaction_plot_data[
      order(interaction_rank),
      pair_display
    ]))
  )]

  interaction_plot <- ggplot(
    interaction_plot_data,
    aes(
      x = mean_abs_interaction,
      y = pair_display,
      fill = interaction_status
    )
  ) +
    geom_col(width = 0.68) +
    geom_errorbar(
      aes(xmin = ci_low, xmax = ci_high),
      width = 0.14,
      linewidth = 0.4
    ) +
    facet_wrap(
      ~scheme_display,
      ncol = paper_facet_columns,
      scales = "free_y"
    ) +
    scale_fill_manual(
      values = c(
        pre_specified = "#D95F0E",
        exploratory = "#756BB1"
      ),
      breaks = c("pre_specified", "exploratory"),
      labels = c("Pre-specified", "Exploratory")
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      x = "Mean absolute pair interaction (BBCH units)",
      y = NULL,
      fill = "Analysis status",
      title = "Out-of-fold TreeSHAP interaction importance",
      subtitle = paste0(
        "Exact interactions on a balanced environment × stage sample;\n",
        "95% cluster-bootstrap intervals"
      )
    ) +
    theme_pheno()

  save_publication_plot(
    interaction_plot,
    "figure_shap_interaction_importance",
    width = 11.8,
    height = 7.8 * paper_facet_rows
  )

  # Direction and developmental timing of the strongest interactions. The
  # ranking above uses absolute magnitude, whereas this profile retains the
  # signed pair contribution. Positive values raise the predicted BBCH stage
  # relative to an additive representation of the two features; negative
  # values lower it. Each point is one environment x exact-stage unit from the
  # balanced interaction sample. Lines and vertical bars show the stage median
  # and interquartile range, respectively.
  top_interaction_keys <- interaction_plot_data[, .(
    scheme,
    pair_id,
    interaction_rank
  )]
  top_interaction_stage_units <- merge(
    interaction_unit,
    top_interaction_keys,
    by = c("scheme", "pair_id")
  )
  top_interaction_stage_units[, pair_display := paste0(
    unname(feature_short_label_lookup[feature_1]),
    " × ",
    unname(feature_short_label_lookup[feature_2])
  )]

  top_interaction_stage_profile <- top_interaction_stage_units[, .(
    n_environment_stage_units = .N,
    mean_signed_interaction = mean(mean_interaction),
    median_signed_interaction = median(mean_interaction),
    q25_signed_interaction = as.numeric(
      quantile(mean_interaction, 0.25, names = FALSE)
    ),
    q75_signed_interaction = as.numeric(
      quantile(mean_interaction, 0.75, names = FALSE)
    ),
    mean_abs_interaction = mean(mean_abs_interaction)
  ), by = .(
    scheme,
    pair_id,
    feature_1,
    feature_2,
    pair_display,
    interaction_rank,
    BBCH_stage
  )]
  setorder(
    top_interaction_stage_profile,
    scheme,
    interaction_rank,
    BBCH_stage
  )
  write_table(
    top_interaction_stage_profile,
    "shap_top_interactions_by_exact_bbch.csv"
  )

  for (scheme_name in schemes) {
    scheme_stage_units <- top_interaction_stage_units[
      scheme == scheme_name
    ]
    scheme_stage_profile <- top_interaction_stage_profile[
      scheme == scheme_name
    ]
    if (!nrow(scheme_stage_units)) next

    pair_levels <- unique(
      scheme_stage_units[
        order(interaction_rank),
        pair_display
      ]
    )
    scheme_stage_units[, pair_display := factor(
      pair_display,
      levels = pair_levels
    )]
    scheme_stage_profile[, pair_display := factor(
      pair_display,
      levels = pair_levels
    )]

    observed_stage_breaks <- sort(unique(scheme_stage_units$BBCH_stage))
    top_interaction_stage_plot <- ggplot(
      scheme_stage_units,
      aes(x = BBCH_stage, y = mean_interaction)
    ) +
      geom_hline(yintercept = 0, color = "grey35", linewidth = 0.45) +
      geom_point(
        color = "grey45",
        alpha = 0.18,
        size = 0.60,
        position = position_jitter(width = 0.45, height = 0)
      ) +
      geom_linerange(
        data = scheme_stage_profile,
        aes(
          x = BBCH_stage,
          ymin = q25_signed_interaction,
          ymax = q75_signed_interaction
        ),
        color = "#0072B2",
        linewidth = 0.65,
        inherit.aes = FALSE
      ) +
      geom_line(
        data = scheme_stage_profile,
        aes(
          x = BBCH_stage,
          y = median_signed_interaction,
          group = 1L
        ),
        color = "#0072B2",
        linewidth = 0.75,
        inherit.aes = FALSE
      ) +
      geom_point(
        data = scheme_stage_profile,
        aes(x = BBCH_stage, y = median_signed_interaction),
        shape = 21,
        fill = "white",
        color = "#0072B2",
        stroke = 0.65,
        size = 2.0,
        inherit.aes = FALSE
      ) +
      facet_wrap(~pair_display, ncol = 2L) +
      scale_x_continuous(
        breaks = observed_stage_breaks,
        expand = expansion(add = 1.5)
      ) +
      labs(
        x = "Observed BBCH stage",
        y = "Signed TreeSHAP interaction (BBCH units)",
        title = paste0(
          "Top ",
          length(pair_levels),
          " interactions across development — ",
          unname(scheme_label_lookup[scheme_name])
        ),
        subtitle = paste0(
          "Grey points: environment × exact-stage units; blue points/lines: ",
          "stage median; vertical bars: interquartile range. Descriptive ",
          "attribution conditional on observed stage"
        )
      ) +
      theme_pheno() +
      theme(panel.grid.major.x = element_line(color = "grey90"))

    save_publication_plot(
      top_interaction_stage_plot,
      paste0(
        "figure_shap_top",
        INTERACTION_TOP_N,
        "_interactions_by_stage__",
        safe_id(scheme_name)
      ),
      width = 11.6,
      height = 3.3 * ceiling(length(pair_levels) / 2L) + 1.7
    )
  }
}

figure_manifest <- data.table(
  figure_stem = c(
    "figure_shap_main_composite",
    "figure_shap_global_importance",
    "figure_shap_stage_importance_share",
    "figure_shap_stage_importance_absolute",
    "figure_shap_focal_bbch75_attribution",
    "figure_shap_summary",
    "figure_shap_biological_groups",
    "figure_shap_dependence__<feature>",
    "figure_shap_dependence_interaction__<scheme>",
    "figure_shap_interaction_importance",
    paste0(
      "figure_shap_top",
      INTERACTION_TOP_N,
      "_interactions_by_stage__<scheme>"
    )
  ),
  suggested_role = c(
    "main manuscript figure",
    "main or supplementary panel",
    "main stage-dynamics panel",
    "supplementary magnitude heatmap",
    "main BBCH 75 diagnostic",
    "supplementary SHAP distribution",
    "supplementary grouped attribution",
    "supplementary dependence series",
    "interaction-diagnostic dependence summary",
    "supplementary interaction ranking",
    "main or supplementary interaction direction panel"
  ),
  generated_when = c(
    "patchwork installed and BBCH 75 present",
    "always",
    "always",
    "always",
    "BBCH 75 present",
    "always",
    "always",
    "for each top-ranked feature",
    "interactions enabled; one plot per scheme",
    "interactions enabled and sampled",
    "interactions enabled and sampled"
  ),
  formats = "PNG (450 dpi) and vector PDF"
)
write_table(figure_manifest, "shap_figure_manifest.csv")

# A compact manuscript composite: global ranking, stage-specific relative
# attribution, and the pre-declared BBCH 75 diagnostic. Individual panels are
# also exported above at full size for supplementary or alternative layouts.
if (
  exists("focal_plot") &&
    requireNamespace("patchwork", quietly = TRUE)
) {
  figure_shap_main <- patchwork::wrap_plots(
    global_plot,
    stage_plot_share,
    focal_plot,
    ncol = 1,
    guides = "collect",
    heights = c(1, 1.35, 1)
  ) +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(
      legend.position = "bottom",
      plot.tag = element_text(face = "bold", size = 12)
    )
  save_publication_plot(
    figure_shap_main,
    "figure_shap_main_composite",
    width = 11.8,
    height = 16.0 + 3.0 * max(0L, paper_facet_rows - 1L)
  )
} else if (exists("focal_plot")) {
  warning(
    "Package 'patchwork' is not installed; individual SHAP figures were ",
    "created, but figure_shap_main_composite was skipped."
  )
}

# ------------------------------------------------------------------------------
# 15. PROVENANCE AND FINAL MESSAGES
# ------------------------------------------------------------------------------

analysis_configuration <- list(
  shap_run_id = SHAP_RUN_ID,
  benchmark_run_id = BENCHMARK_RUN_ID,
  benchmark_dir = normalizePath(BENCHMARK_DIR),
  input_file = normalizePath(input_file),
  input_md5 = input_md5,
  feature_set = feature_set,
  features = features,
  feature_groups = feature_groups,
  learner_id = LEARNER_ID,
  schemes = schemes,
  scheme_estimands = scheme_estimands,
  cross_scheme_pooling = FALSE,
  seed = SEED,
  reconstruction_engine = "mlr3learners::regr.xgboost",
  configured_features = features,
  model_feature_order = model_features,
  package_version_audit = package_version_audit,
  benchmark_xgboost_version = benchmark_xgboost_version,
  current_xgboost_version = current_xgboost_version,
  model_reconstruction_source = list(
    split_rules = normalizePath(split_specs_path),
    tuned_parameters = normalizePath(model_specifications_path),
    prediction_reference = normalizePath(oof_predictions_path)
  ),
  prediction_match_tolerance = PREDICTION_MATCH_TOLERANCE,
  shap_additivity_tolerance = SHAP_ADDITIVITY_TOLERANCE,
  shap_additivity_relative_tolerance =
    SHAP_ADDITIVITY_RELATIVE_TOLERANCE,
  interaction_additivity_tolerance =
    INTERACTION_ADDITIVITY_TOLERANCE,
  interaction_additivity_relative_tolerance =
    INTERACTION_ADDITIVITY_RELATIVE_TOLERANCE,
  run_interactions = RUN_INTERACTIONS,
  interaction_max_rows_per_scheme =
    INTERACTION_MAX_ROWS_PER_SCHEME,
  interaction_sample_definition = paste0(
    "one representative per environment x exact BBCH, round-robin ",
    "balanced over outer fold x stage, then capped per scheme"
  ),
  summary_plot_max_rows_per_scheme =
    SUMMARY_PLOT_MAX_ROWS_PER_SCHEME,
  dependence_top_n = DEPENDENCE_TOP_N,
  interaction_plot_top_n = INTERACTION_TOP_N,
  interaction_coloured_dependence_definition = paste0(
    "unstacked feature-value versus SHAP scatter coloured by the ",
    "strongest exact interaction partner within each CV scheme"
  ),
  interaction_stage_profile_definition = paste0(
    "signed pair contribution over observed exact BBCH; environment-stage ",
    "points with stage median and interquartile range"
  ),
  focal_bbch_stages = FOCAL_BBCH_STAGES,
  bootstrap_replicates = BOOTSTRAP_REPLICATES,
  bootstrap_cluster_rule = list(
    spatial_block = "physical site",
    leave_one_year_out = "season year",
    expanding_origin = "season year",
    environmental_zone = "environmental zone",
    joint_space_time = "joint outer space-time fold"
  ),
  primary_shap_weighting =
    "macro mean over environment x exact observed BBCH units",
  grouped_importance_definition =
    "absolute value after summing signed feature SHAP within group",
  interaction_pair_definition =
    "two times one symmetric off-diagonal interaction cell",
  checkpoint_mode = CHECKPOINT_MODE,
  checkpoint_version = CHECKPOINT_VERSION,
  workers = workers,
  n_cohort_rows = nrow(cohort),
  n_oof_rows_by_scheme = row_shap[, .N, by = scheme],
  n_environments = uniqueN(cohort$env_id),
  n_sites = uniqueN(cohort$site_id),
  n_years = uniqueN(cohort$season_year)
)

atomic_save_rds(
  analysis_configuration,
  file.path(OUTPUT_DIR, "shap_analysis_configuration.rds")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "session_info.txt")
)

message("Out-of-fold TreeSHAP analysis finished.")
message("Results written to: ", OUTPUT_DIR)
message(
  "Primary importance table: ",
  file.path(
    OUTPUT_DIR,
    "global_shap_importance_macro_environment_stage.csv"
  )
)
message(
  "Important: interpret SHAP as model attribution, not causal effect. ",
  "Any XAI-guided interaction constraints must be evaluated in a new ",
  "nested experiment or on independent BRIWECS data, not on these same ",
  "OOF explanations."
)

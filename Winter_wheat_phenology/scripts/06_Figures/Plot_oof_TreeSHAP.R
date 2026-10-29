# ==============================================================================
# PUBLICATION PLOTS FROM EXISTING OUT-OF-FOLD TREESHAP RESULTS
# MAIN: JOINT SPACE-TIME + EXPANDING ORIGIN; SUPPLEMENT: OTHER SAVED SCHEMES
# ==============================================================================
#
# Purpose
# -------
# Create four publication-oriented figure groups from the exact out-of-fold
# TreeSHAP objects written by:
#   shap_xgboost_phenology_enhanced_plots_check_points.R
#
# This is deliberately a plotting-only script. The source SHAP run used
# n = 10,000 for additive and interaction samples. It does not reconstruct folds,
# refit XGBoost models, or recalculate SHAP values.
#
# Required source files in SHAP_RESULTS_DIR
# -----------------------------------------
#   shap_analysis_configuration.rds
#   oof_shap_values_wide.rds
#   oof_shap_interactions_sample.rds
#
# The source SHAP run must contain the two main schemes. Any additional
# schemes with both SHAP values and interactions are plotted in the supplement.
#
# To plot supplemental schemes, the source SHAP run must include them in both
# saved SHAP values and interaction samples. Enable interactions in that run.
#
# Figures
# -------
#   1. SHAPforxgboost-style out-of-fold TreeSHAP summary. Supplemental
#      schemes are combined in one labelled figure.
#   2. Relative feature-attribution heatmap across observed BBCH stages.
#   3. Five strongest TreeSHAP pairs plus a dependence plot of the strongest
#      pair, shown separately for both cross-validation schemes.
#   4. Dependence plots for the six most important features in each scheme.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. USER CONFIGURATION
# ------------------------------------------------------------------------------

PROJECT_DIR <- Sys.getenv(
  "PHENOLOGY_PROJECT_DIR",
  unset = "/data/analysis3/milekj/UseCase_WW02_2024"
)


SHAP_RUN_ID <- Sys.getenv(
  "PHENOLOGY_SHAP_RUN_ID",
  unset = "xgboost_generalisation_oof_treeshap_paper_cum"
)

SHAP_RESULTS_DIR <- Sys.getenv(
  "PHENOLOGY_SHAP_RESULTS_DIR",
  unset = file.path(PROJECT_DIR, "shap_results", SHAP_RUN_ID)
)

OUTPUT_DIR <- Sys.getenv(
  "PHENOLOGY_SHAP_FIGURE_OUTPUT_DIR",
  unset = file.path(SHAP_RESULTS_DIR, "joint_expanding_and_supplement_plots")
)
PLOT_DIR <- file.path(OUTPUT_DIR, "plots")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")

# These schemes remain the main figures; other saved schemes form the supplement.
MAIN_SCHEMES <- c("joint_space_time", "expanding_origin")

SEED <- 123L
TOP_N_FEATURES <- 6L
TOP_N_INTERACTIONS <- 5L
WIDE_FIGURE_WIDTH <- 10.5

read_env_integer <- function(name, default, minimum = 1L) {
  value <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) return(as.integer(default))
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < minimum) {
    stop(name, " must be an integer >= ", minimum, ".")
  }
  parsed
}

SUMMARY_MAX_UNITS_PER_SCHEME <- read_env_integer(
  "PHENOLOGY_SHAP_SUMMARY_MAX_ROWS_PER_SCHEME",
  10000L
)
DEPENDENCE_MAX_UNITS_PER_SCHEME <- read_env_integer(
  "PHENOLOGY_SHAP_DEPENDENCE_MAX_ROWS_PER_SCHEME",
  10000L
)


# ------------------------------------------------------------------------------
# 2. PACKAGES AND REPRODUCIBILITY
# ------------------------------------------------------------------------------

SERVER_R_LIBRARY <- Sys.getenv(
  "PHENOLOGY_R_LIBRARY",
  unset = "/data/analysis3/milekj/Rlibs/4.5"
)
if (dir.exists(SERVER_R_LIBRARY)) {
  .libPaths(unique(c(SERVER_R_LIBRARY, .libPaths())))
}

if (requireNamespace("renv", quietly = TRUE)) {
  try(
    renv::load(project = PROJECT_DIR),
    silent = TRUE
  )
}

required_packages <- c(
  "data.table",
  "ggplot2",
  "patchwork",
  "scales",
  "mgcv",
  "systemfonts"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install these packages before plotting: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

HAS_SHAPFORXGBOOST <- requireNamespace(
  "SHAPforxgboost",
  quietly = TRUE
)

RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)
data.table::setDTthreads(1L)

if (!capabilities("cairo")) {
  stop("Cairo graphics support is required for Liberation Sans PDF and EPS export.")
}
if (!any(tolower(trimws(systemfonts::system_fonts()$family)) == "liberation sans")) {
  stop("Liberation Sans is not available to R on this system.")
}

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 3. GENERAL HELPERS
# ------------------------------------------------------------------------------

assert_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      label,
      " was not found: ",
      path,
      "\nRun shap_xgboost_phenology_enhanced_plots_check_points.R first with ",
      "PHENOLOGY_SHAP_SCHEMES=\"joint_space_time,expanding_origin\" and ",
      "PHENOLOGY_SHAP_RUN_INTERACTIONS=true."
    )
  }
  invisible(path)
}

normalise_scheme <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[- ]", "_", x)
  x[x %in% c("joint", "joint_spacetime", "joint_space_time_cv")] <-
    "joint_space_time"
  x[x %in% c(
    "expanding",
    "expanding_window",
    "expanding_origin_cv",
    "rolling_origin",
    "forward_chaining"
  )] <- "expanding_origin"
  x
}

scheme_label <- function(x) {
  labels <- c(
    joint_space_time = "Joint space-time",
    expanding_origin = "Expanding origin",
    spatial = "Spatial blocks",
    temporal_block = "Temporal blocks"
  )
  result <- unname(labels[as.character(x)])
  missing <- is.na(result)
  result[missing] <- gsub("_", " ", as.character(x)[missing])
  result
}

safe_id <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
}

write_table <- function(x, filename) {
  fwrite(
    as.data.table(x),
    file.path(TABLE_DIR, filename),
    na = "NA"
  )
}

atomic_save_rds <- function(object, path, compress = FALSE) {
  temporary <- paste0(path, ".tmp_", Sys.getpid())
  saveRDS(object, temporary, compress = compress)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not atomically write: ", path)
  }
  invisible(path)
}

save_plot <- function(plot, stem, width, height) {
  pdf_device <- grDevices::cairo_pdf
  
  ggsave(
    file.path(PLOT_DIR, paste0(stem, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = pdf_device,
    limitsize = FALSE
  )
  ggsave(
    file.path(PLOT_DIR, paste0(stem, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 450,
    bg = "white",
    limitsize = FALSE
  )
  eps_device <- grDevices::cairo_ps
  ggsave(
    file.path(PLOT_DIR, paste0(stem, ".eps")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = eps_device,
    onefile = FALSE,
    bg = "white",
    limitsize = FALSE
  )
  invisible(plot)
}

theme_pheno <- function(base_size = 10) {
  theme_minimal(base_size = base_size, base_family = "Liberation Sans") +
    theme(
      text = element_text(family = "Liberation Sans"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        colour = "grey90",
        linewidth = 0.25
      ),
      panel.grid.major.y = element_blank(),
      axis.text = element_text(colour = "black"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      
      
      plot.tag = element_text(family = "Liberation Sans", face = "bold"),
      plot.margin = margin(8, 10, 8, 8)
    )
}

clip_quantiles <- function(x, lower = 0.01, upper = 0.99) {
  limits <- quantile(
    x,
    probs = c(lower, upper),
    na.rm = TRUE,
    names = FALSE,
    type = 8
  )
  if (!all(is.finite(limits)) || limits[1L] == limits[2L]) return(x)
  pmin(pmax(x, limits[1L]), limits[2L])
}

scale_within_panel <- function(x) {
  limits <- quantile(
    x,
    probs = c(0.01, 0.99),
    na.rm = TRUE,
    names = FALSE,
    type = 8
  )
  if (!all(is.finite(limits)) || limits[1L] == limits[2L]) {
    return(rep(0.5, length(x)))
  }
  pmin(1, pmax(0, (x - limits[1L]) / diff(limits)))
}

# The source SHAP script defines the primary statistical unit as one
# environment x exact observed BBCH stage. This helper caps a plotting sample
# by cycling over fold x stage strata, preserving broad fold/stage coverage.
balanced_unit_sample <- function(data, maximum_rows, seed) {
  answer <- copy(as.data.table(data))
  if (!nrow(answer) || nrow(answer) <= maximum_rows) return(answer)
  
  set.seed(seed)
  answer[, draw_rank__ := sample.int(.N), by = .(fold_id, BBCH_stage)]
  setorder(answer, draw_rank__, fold_id, BBCH_stage)
  answer <- answer[seq_len(maximum_rows)]
  answer[, draw_rank__ := NULL]
  answer
}

feature_label_lookup <- c(
  cumGDD_from_sow_Soltani = "Growing degree-days (cumGDD)",
  cumDTU_beta = "Beta-response thermal time (cumDTU)",
  cumdevfun_phase = "Process-derived development",
  sowing_DOY = "Sowing day of year",
  Daylength = "Day length",
  photoperiod_civil = "Photoperiod",
  cumRR = "Cumulative precipitation",
  cumRad = "Cumulative radiation",
  cumCDD = "Cumulative cold degree-days",
  cumHDD = "Cumulative heat degree-days",
  paw_fc33_0_28_mm = "Plant-available water capacity (0-28 cm)",
  paw_fc33_0_100_mm = "Plant-available water capacity (0-100 cm)",
  FTSW_0_100_proxy = "Fraction of transpirable soil water",
  sumRR_30 = "Precipitation (30 d)",
  sumRad_30 = "Radiation (30 d)",
  sumCDD_30 = "Cold degree-days (30 d)",
  sumHDD_30 = "Heat degree-days (30 d)",
  sumPEvap_30 = "Potential evaporation (30 d)",
  climaticWaterBalance_30 = "Climatic water balance (30 d)",
  meanSoilWater_0_100_30 = "Mean soil water (0-100 cm, 30 d)"
)

feature_label <- function(x) {
  result <- unname(feature_label_lookup[as.character(x)])
  missing <- is.na(result)
  result[missing] <- gsub("_", " ", as.character(x)[missing])
  result
}

# Wrap only figure labels; keep full, single-line names for captions and tables.
feature_plot_label <- function(x, width = 26L) {
  vapply(feature_label(x), function(label) {
    wrapped <- strwrap(label, width = width)
    if (!length(wrapped)) return(label)
    paste(wrapped, collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

# Wrap a complete interaction pair to keep the bar-chart label column narrow.
pair_plot_label <- function(feature_1, feature_2, width = 26L) {
  vapply(seq_along(feature_1), function(i) {
    label <- paste(
      feature_label(feature_1[i]),
      feature_label(feature_2[i]),
      sep = " × "
    )
    paste(strwrap(label, width = width), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}


# ------------------------------------------------------------------------------
# 4. LOAD THE AUTHORITATIVE SHAP-RUN OUTPUTS
# ------------------------------------------------------------------------------

analysis_configuration_path <- file.path(
  SHAP_RESULTS_DIR,
  "shap_analysis_configuration.rds"
)
row_shap_path <- file.path(
  SHAP_RESULTS_DIR,
  "oof_shap_values_wide.rds"
)
interaction_path <- file.path(
  SHAP_RESULTS_DIR,
  "oof_shap_interactions_sample.rds"
)

assert_file(
  analysis_configuration_path,
  "SHAP analysis configuration"
)
assert_file(row_shap_path, "Out-of-fold SHAP values")
assert_file(interaction_path, "Out-of-fold SHAP interactions")

analysis_configuration <- readRDS(analysis_configuration_path)
row_shap <- as.data.table(readRDS(row_shap_path))
interaction_all <- as.data.table(readRDS(interaction_path))

if (!"scheme" %in% names(row_shap)) {
  stop("oof_shap_values_wide.rds has no scheme column.")
}
if (!"scheme" %in% names(interaction_all)) {
  stop("oof_shap_interactions_sample.rds has no scheme column.")
}

row_shap[, scheme := normalise_scheme(scheme)]
interaction_all[, scheme := normalise_scheme(scheme)]

missing_main_schemes <- setdiff(MAIN_SCHEMES, unique(row_shap$scheme))
missing_interaction_schemes <- setdiff(
  MAIN_SCHEMES,
  unique(interaction_all$scheme)
)
if (length(missing_main_schemes) || length(missing_interaction_schemes)) {
  stop(
    "The saved SHAP run is missing a required main scheme.\n",
    "Missing in oof_shap_values_wide.rds: ",
    if (length(missing_main_schemes)) paste(missing_main_schemes, collapse = ", ") else "none",
    "\nMissing in oof_shap_interactions_sample.rds: ",
    if (length(missing_interaction_schemes)) paste(missing_interaction_schemes, collapse = ", ") else "none"
  )
}

supplement_candidates <- sort(setdiff(
  union(unique(row_shap$scheme), unique(interaction_all$scheme)),
  MAIN_SCHEMES
))
missing_supplement_data <- supplement_candidates[
  !supplement_candidates %in% unique(row_shap$scheme) |
    !supplement_candidates %in% unique(interaction_all$scheme)
]
if (length(missing_supplement_data)) {
  warning(
    "Skipping supplemental schemes without both SHAP values and interactions: ",
    paste(missing_supplement_data, collapse = ", ")
  )
}
SUPPLEMENT_SCHEMES <- setdiff(supplement_candidates, missing_supplement_data)
SCHEMES <- c(MAIN_SCHEMES, SUPPLEMENT_SCHEMES)
if (!length(SUPPLEMENT_SCHEMES)) {
  warning(
    "No additional schemes with both SHAP values and interactions were found; ",
    "supplemental figures cannot be made from this saved run."
  )
}

row_shap <- row_shap[scheme %in% SCHEMES]
interaction_all <- interaction_all[scheme %in% SCHEMES]

shap_columns_found <- grep("^shap__", names(row_shap), value = TRUE)
x_columns_found <- grep("^x__", names(row_shap), value = TRUE)
features_from_shap <- sub("^shap__", "", shap_columns_found)
features_from_x <- sub("^x__", "", x_columns_found)

configured_features <- as.character(analysis_configuration$features)
if (!length(configured_features)) {
  configured_features <- features_from_shap
}

if (!setequal(configured_features, features_from_shap)) {
  stop(
    "Feature names in shap_analysis_configuration.rds and ",
    "oof_shap_values_wide.rds differ."
  )
}
if (!setequal(configured_features, features_from_x)) {
  stop(
    "Predictor-value and SHAP columns do not contain the same feature set."
  )
}

FEATURES <- configured_features
shap_columns <- paste0("shap__", FEATURES)
x_columns <- paste0("x__", FEATURES)

required_row_columns <- c(
  "obs_id",
  "env_id",
  "BBCH_stage",
  "scheme",
  "fold_id",
  shap_columns,
  x_columns
)
missing_row_columns <- setdiff(required_row_columns, names(row_shap))
if (length(missing_row_columns)) {
  stop(
    "oof_shap_values_wide.rds is missing: ",
    paste(missing_row_columns, collapse = ", ")
  )
}

required_interaction_columns <- c(
  "obs_id",
  "env_id",
  "BBCH_stage",
  "scheme",
  "fold_id",
  "pair_id",
  "feature_1",
  "feature_2",
  "feature_1_value",
  "feature_2_value",
  "interaction_value"
)
missing_interaction_columns <- setdiff(
  required_interaction_columns,
  names(interaction_all)
)
if (length(missing_interaction_columns)) {
  stop(
    "oof_shap_interactions_sample.rds is missing: ",
    paste(missing_interaction_columns, collapse = ", ")
  )
}

if (anyDuplicated(row_shap[, .(scheme, obs_id)])) {
  stop("An observation has duplicate out-of-fold SHAP rows within a scheme.")
}
if (any(!is.finite(row_shap$BBCH_stage))) {
  stop("BBCH_stage contains non-finite values.")
}
if (any(!is.finite(interaction_all$interaction_value))) {
  stop("interaction_value contains non-finite values.")
}

source_scheme_counts <- row_shap[, .(
  n_oof_rows = .N,
  n_environments = uniqueN(env_id),
  n_folds = uniqueN(fold_id),
  n_stages = uniqueN(BBCH_stage)
), by = scheme]
write_table(source_scheme_counts, "source_shap_scheme_counts.csv")

# Interaction SHAP was sampled upstream; report its actual observation count.
source_interaction_counts <- interaction_all[, .(
  n_interaction_rows = .N,
  n_sampled_observations = uniqueN(obs_id),
  n_pairs = uniqueN(pair_id)
), by = scheme]
write_table(source_interaction_counts, "source_interaction_scheme_counts.csv")


# ------------------------------------------------------------------------------
# 5. ENVIRONMENT x STAGE AGGREGATION
# ------------------------------------------------------------------------------

# Mean predictor and signed SHAP values for dependence and summary plots.
unit_values <- row_shap[, lapply(.SD, mean), by = .(
  scheme,
  fold_id,
  env_id,
  BBCH_stage
), .SDcols = c(x_columns, shap_columns)]

# Mean absolute attribution is calculated before aggregation, matching the
# enhanced SHAP script's environment x exact-stage macro weighting.
shap_long <- melt(
  row_shap,
  id.vars = c("scheme", "fold_id", "env_id", "BBCH_stage"),
  measure.vars = shap_columns,
  variable.name = "feature",
  value.name = "shap_value",
  variable.factor = FALSE
)
shap_long[, feature := sub("^shap__", "", feature)]

feature_unit <- shap_long[, .(
  mean_abs_shap = mean(abs(shap_value)),
  mean_shap = mean(shap_value),
  n_observations = .N
), by = .(
  scheme,
  fold_id,
  env_id,
  BBCH_stage,
  feature
)]

global_importance <- feature_unit[, .(
  mean_abs_shap = mean(mean_abs_shap),
  mean_shap = mean(mean_shap),
  n_environment_stage_units = .N
), by = .(scheme, feature)]
global_importance[, importance_rank := frank(
  -mean_abs_shap,
  ties.method = "first"
), by = scheme]
global_importance[, importance_share := mean_abs_shap / sum(mean_abs_shap),
                  by = scheme]
setorder(global_importance, scheme, importance_rank)

stage_attribution <- feature_unit[, .(
  mean_abs_shap = mean(mean_abs_shap),
  mean_shap = mean(mean_shap),
  n_environment_stage_units = .N
), by = .(scheme, BBCH_stage, feature)]
stage_attribution[, relative_attribution :=
                    mean_abs_shap / sum(mean_abs_shap),
                  by = .(scheme, BBCH_stage)]
stage_attribution[, stage_rank := frank(
  -mean_abs_shap,
  ties.method = "first"
), by = .(scheme, BBCH_stage)]
setorder(stage_attribution, scheme, BBCH_stage, stage_rank)

write_table(
  global_importance,
  "oof_shap_global_importance_environment_stage.csv"
)
write_table(
  stage_attribution,
  "oof_shap_relative_stage_attribution.csv"
)


# ------------------------------------------------------------------------------
# 6. INTERACTION AGGREGATION
# ------------------------------------------------------------------------------

interaction_unit <- interaction_all[, .(
  mean_abs_interaction = mean(abs(interaction_value)),
  mean_interaction = mean(interaction_value),
  feature_1_value = mean(feature_1_value),
  feature_2_value = mean(feature_2_value),
  n_sampled_observations = .N
), by = .(
  scheme,
  fold_id,
  env_id,
  BBCH_stage,
  pair_id,
  feature_1,
  feature_2
)]

interaction_importance <- interaction_unit[, .(
  mean_abs_interaction = mean(mean_abs_interaction),
  mean_interaction = mean(mean_interaction),
  n_environment_stage_units = .N
), by = .(
  scheme,
  pair_id,
  feature_1,
  feature_2
)]
interaction_importance[, interaction_rank := frank(
  -mean_abs_interaction,
  ties.method = "first"
), by = scheme]
interaction_importance[, interaction_share :=
                         mean_abs_interaction / sum(mean_abs_interaction),
                       by = scheme]
setorder(interaction_importance, scheme, interaction_rank)

write_table(
  interaction_importance,
  "oof_shap_interaction_importance_environment_stage.csv"
)


# ------------------------------------------------------------------------------
# 7. FIGURE 1: SHAPFORXGBOOST-STYLE OOF SUMMARY
# ------------------------------------------------------------------------------

make_fallback_summary <- function(sample_data, scheme_name) {
  scheme_id <- scheme_name
  shap_long_plot <- melt(
    sample_data,
    id.vars = c("fold_id", "env_id", "BBCH_stage"),
    measure.vars = shap_columns,
    variable.name = "feature",
    value.name = "shap_value",
    variable.factor = FALSE
  )
  shap_long_plot[, feature := sub("^shap__", "", feature)]
  
  x_long_plot <- melt(
    sample_data,
    id.vars = c("fold_id", "env_id", "BBCH_stage"),
    measure.vars = x_columns,
    variable.name = "feature",
    value.name = "feature_value",
    variable.factor = FALSE
  )
  x_long_plot[, feature := sub("^x__", "", feature)]
  
  plot_data <- merge(
    shap_long_plot,
    x_long_plot,
    by = c("fold_id", "env_id", "BBCH_stage", "feature")
  )
  plot_data[, feature_value_scaled := scale_within_panel(feature_value),
            by = feature]
  
  feature_order <- global_importance[
    scheme == scheme_id
  ][order(importance_rank), feature]
  plot_data[, feature_display := factor(
    feature_plot_label(feature),
    levels = rev(feature_plot_label(feature_order))
  )]
  
  ggplot(
    plot_data,
    aes(
      x = shap_value,
      y = feature_display,
      colour = feature_value_scaled
    )
  ) +
    geom_vline(
      xintercept = 0,
      colour = "grey55",
      linewidth = 0.35
    ) +
    geom_point(
      position = position_jitter(width = 0, height = 0.16),
      size = 0.65,
      alpha = 0.50
    ) +
    scale_colour_gradientn(
      colours = c("#2166AC", "#BDBDBD", "#B2182B"),
      values = c(0, 0.5, 1),
      limits = c(0, 1),
      breaks = c(0, 1),
      labels = c("Low", "High")
    ) +
    labs(
      x = "SHAP value (BBCH units)",
      y = NULL,
      colour = "Feature value"
    ) +
    theme_pheno(10)
}

make_summary_plot <- function(scheme_name, scheme_number) {
  scheme_id <- scheme_name
  sample_data <- balanced_unit_sample(
    unit_values[scheme == scheme_id],
    maximum_rows = SUMMARY_MAX_UNITS_PER_SCHEME,
    seed = SEED + scheme_number
  )
  
  if (!HAS_SHAPFORXGBOOST) {
    return(make_fallback_summary(sample_data, scheme_name))
  }
  
  S <- as.matrix(sample_data[, ..shap_columns])
  X_table <- sample_data[, ..x_columns]
  X_table <- X_table[, lapply(.SD, clip_quantiles)]
  X <- as.matrix(X_table)
  
  pretty_names <- make.unique(feature_plot_label(FEATURES), sep = " ")
  colnames(S) <- pretty_names
  colnames(X) <- pretty_names
  
  tryCatch(
    {
      shap_long_package <- SHAPforxgboost::shap.prep(
        shap_contrib = S,
        X_train = X
      )
      SHAPforxgboost::shap.plot.summary(
        shap_long_package,
        scientific = FALSE
      ) +
        labs(
          x = "SHAP value (BBCH units)",
          y = NULL
        ) +
        theme_pheno(10)
    },
    error = function(error) {
      warning(
        "SHAPforxgboost plotting failed for ",
        scheme_name,
        "; using the equivalent ggplot layout. Error: ",
        conditionMessage(error)
      )
      make_fallback_summary(sample_data, scheme_name)
    }
  )
}

summary_plots <- lapply(
  seq_along(MAIN_SCHEMES),
  function(i) make_summary_plot(MAIN_SCHEMES[i], i)
)

figure_summary <- wrap_plots(summary_plots, ncol = 1L) +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "(",
    tag_suffix = ")"
  )

save_plot(
  figure_summary,
  "figure_1_oof_treeshap_summary_joint_expanding",
  width = 7.2,
  height = 10.0
)

if (length(SUPPLEMENT_SCHEMES)) {
  supplement_summary_plots <- lapply(
    seq_along(SUPPLEMENT_SCHEMES),
    function(i) make_summary_plot(
      SUPPLEMENT_SCHEMES[i],
      i + length(MAIN_SCHEMES)
    )
  )
  figure_summary_supplement <- wrap_plots(
    supplement_summary_plots,
    ncol = 1L
  ) + plot_annotation(
    tag_levels = "A",
    tag_prefix = "(",
    tag_suffix = ")"
  )
  save_plot(
    figure_summary_supplement,
    "figure_1_oof_treeshap_summary_supplement",
    width = 7.2,
    height = 3.5 * length(SUPPLEMENT_SCHEMES)
  )
}


# ------------------------------------------------------------------------------
# 8. FIGURE 2: RELATIVE STAGE ATTRIBUTION
# ------------------------------------------------------------------------------

make_stage_attribution_plot <- function(schemes) {
  feature_order <- global_importance[scheme %in% schemes, .(
    mean_rank = mean(importance_rank)
  ), by = feature][order(mean_rank), feature]
  
  stage_plot_data <- copy(stage_attribution[scheme %in% schemes])
  stage_plot_data[, scheme_display := factor(
    scheme_label(scheme),
    levels = scheme_label(schemes)
  )]
  stage_plot_data[, feature_display := factor(
    feature_plot_label(feature),
    levels = rev(feature_plot_label(feature_order))
  )]
  stage_plot_data[, stage_display := factor(
    BBCH_stage,
    levels = sort(unique(BBCH_stage))
  )]
  
  maximum_relative_attribution <- max(
    stage_plot_data$relative_attribution,
    na.rm = TRUE
  )
  
  ggplot(
    stage_plot_data,
    aes(
      x = stage_display,
      y = feature_display,
      fill = relative_attribution
    )
  ) +
    geom_tile(
      width = 0.94,
      height = 0.90,
      colour = "white",
      linewidth = 0.18
    ) +
    facet_wrap(~scheme_display, ncol = 1L) +
    scale_fill_viridis_c(
      option = "H",
      direction = 1,
      limits = c(0, maximum_relative_attribution),
      labels = scales::label_percent(accuracy = 1),
      guide = guide_colorbar(
        title.position = "top",
        barwidth = 15,
        barheight = 0.55
      )
    ) +
    labs(
      x = "Observed BBCH stage",
      y = NULL,
      fill = "Relative absolute SHAP attribution"
    ) +
    theme_pheno(10) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      axis.text.x = element_text(size = 9)
    )
}

figure_stage <- make_stage_attribution_plot(MAIN_SCHEMES)
save_plot(
  figure_stage,
  "figure_2_relative_stage_attribution_joint_expanding",
  width = WIDE_FIGURE_WIDTH,
  height = 8.5
)

if (length(SUPPLEMENT_SCHEMES)) {
  figure_stage_supplement <- make_stage_attribution_plot(SUPPLEMENT_SCHEMES)
  save_plot(
    figure_stage_supplement,
    "figure_2_relative_stage_attribution_supplement",
    width = WIDE_FIGURE_WIDTH,
    height = 4.25 * length(SUPPLEMENT_SCHEMES)
  )
}


# ------------------------------------------------------------------------------
# 9. DEPENDENCE-PLOT HELPERS
# ------------------------------------------------------------------------------

strongest_partner <- function(feature_name, scheme_name) {
  scheme_id <- scheme_name
  feature_id <- feature_name
  candidates <- interaction_importance[
    scheme == scheme_id &
      (feature_1 == feature_id | feature_2 == feature_id)
  ][order(interaction_rank)]
  
  if (!nrow(candidates)) return(NA_character_)
  if (candidates$feature_1[1L] == feature_name) {
    candidates$feature_2[1L]
  } else {
    candidates$feature_1[1L]
  }
}

dependence_samples <- setNames(
  lapply(seq_along(SCHEMES), function(i) {
    scheme_id <- SCHEMES[i]
    balanced_unit_sample(
      unit_values[scheme == scheme_id],
      maximum_rows = DEPENDENCE_MAX_UNITS_PER_SCHEME,
      seed = SEED + 1000L + i
    )
  }),
  SCHEMES
)

make_feature_dependence_plot <- function(
    scheme_name,
    feature_name,
    partner_name
) {
  sample_data <- dependence_samples[[scheme_name]]
  plot_data <- data.table(
    feature_value = sample_data[[paste0("x__", feature_name)]],
    shap_value = sample_data[[paste0("shap__", feature_name)]],
    partner_value = sample_data[[paste0("x__", partner_name)]],
    BBCH_stage = sample_data$BBCH_stage
  )
  plot_data[, feature_value_plot := clip_quantiles(feature_value)]
  plot_data[, partner_value_plot := clip_quantiles(partner_value)]
  
  smooth_layer <- if (uniqueN(plot_data$feature_value_plot) >= 8L) {
    geom_smooth(
      aes(group = 1),
      method = "gam",
      formula = y ~ s(x, bs = "cs", k = 8),
      se = FALSE,
      colour = "black",
      linewidth = 0.7
    )
  } else {
    NULL
  }
  
  ggplot(
    plot_data,
    aes(
      x = feature_value_plot,
      y = shap_value,
      colour = partner_value_plot
    )
  ) +
    geom_hline(
      yintercept = 0,
      colour = "grey55",
      linewidth = 0.35
    ) +
    geom_point(size = 0.72, alpha = 0.40) +
    smooth_layer +
    scale_colour_viridis_c(option = "C", end = 0.92) +
    labs(
      x = feature_plot_label(feature_name),
      y = "OOF SHAP value (BBCH units)",
      colour = feature_plot_label(partner_name)
    ) +
    theme_pheno(9) +
    theme(
      legend.position = "right"
    )
}

make_interaction_dependence_plot <- function(
    scheme_name,
    strongest_pair,
    x_feature
) {
  scheme_id <- scheme_name
  pair_id_value <- strongest_pair$pair_id[1L]
  pair_data <- copy(interaction_unit[
    scheme == scheme_id & pair_id == pair_id_value
  ])
  
  if (x_feature == strongest_pair$feature_1[1L]) {
    colour_feature <- strongest_pair$feature_2[1L]
    pair_data[, x_value := feature_1_value]
    pair_data[, colour_value := feature_2_value]
  } else {
    colour_feature <- strongest_pair$feature_1[1L]
    pair_data[, x_value := feature_2_value]
    pair_data[, colour_value := feature_1_value]
  }
  
  pair_data[, x_value_plot := clip_quantiles(x_value)]
  pair_data[, colour_value_plot := clip_quantiles(colour_value)]
  
  smooth_layer <- if (uniqueN(pair_data$x_value_plot) >= 8L) {
    geom_smooth(
      aes(group = 1),
      method = "gam",
      formula = y ~ s(x, bs = "cs", k = 8),
      se = FALSE,
      colour = "black",
      linewidth = 0.75
    )
  } else {
    NULL
  }
  
  ggplot(
    pair_data,
    aes(
      x = x_value_plot,
      y = mean_interaction,
      colour = colour_value_plot
    )
  ) +
    geom_hline(
      yintercept = 0,
      colour = "grey55",
      linewidth = 0.35
    ) +
    geom_point(size = 0.78, alpha = 0.43) +
    smooth_layer +
    scale_colour_viridis_c(option = "C", end = 0.92) +
    labs(
      x = feature_plot_label(x_feature),
      y = "TreeSHAP interaction (BBCH units)",
      colour = feature_plot_label(colour_feature)
    ) +
    theme_pheno(9) +
    theme(
      legend.position = "right"
    )
}


# ------------------------------------------------------------------------------
# 10. FIGURE 3: TOP-FIVE INTERACTIONS + STRONGEST-PAIR DEPENDENCE
# ------------------------------------------------------------------------------

make_interaction_figure <- function(schemes) {
  interaction_panels <- list()
  
  for (scheme_number in seq_along(schemes)) {
    scheme_name <- schemes[scheme_number]
    scheme_id <- scheme_name
    
    top_pairs <- interaction_importance[
      scheme == scheme_id & interaction_rank <= TOP_N_INTERACTIONS
    ][order(interaction_rank)]
    if (!nrow(top_pairs)) {
      stop("No interaction pairs are available for ", scheme_name, ".")
    }
    
    top_pairs[, pair_display := pair_plot_label(feature_1, feature_2)]
    top_pairs[, pair_display := factor(
      pair_display,
      levels = rev(pair_display)
    )]
    
    bar_plot <- ggplot(
      top_pairs,
      aes(x = pair_display, y = mean_abs_interaction)
    ) +
      geom_col(fill = "#0072B2", width = 0.70) +
      coord_flip() +
      scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
      labs(
        x = NULL,
        y = "Mean absolute pair interaction (BBCH units)"
      ) +
      theme_pheno(9)
    
    strongest_pair <- top_pairs[1L]
    member_features <- c(
      strongest_pair$feature_1[1L],
      strongest_pair$feature_2[1L]
    )
    member_importance <- global_importance[
      scheme == scheme_id & feature %in% member_features
    ]
    setorder(member_importance, importance_rank)
    x_feature <- member_importance$feature[1L]
    
    interaction_dependence <- make_interaction_dependence_plot(
      scheme_name = scheme_name,
      strongest_pair = strongest_pair,
      x_feature = x_feature
    )
    
    interaction_panels[[length(interaction_panels) + 1L]] <- bar_plot
    interaction_panels[[length(interaction_panels) + 1L]] <-
      interaction_dependence
  }
  
  wrap_plots(
    interaction_panels,
    ncol = 2L,
    widths = c(0.90, 1.10)
  ) + plot_annotation(
    tag_levels = "A",
    tag_prefix = "(",
    tag_suffix = ")"
  )
}

figure_interactions <- make_interaction_figure(MAIN_SCHEMES)
save_plot(
  figure_interactions,
  "figure_3_oof_treeshap_interactions_joint_expanding",
  width = WIDE_FIGURE_WIDTH,
  height = 8.5
)

if (length(SUPPLEMENT_SCHEMES)) {
  figure_interactions_supplement <- make_interaction_figure(SUPPLEMENT_SCHEMES)
  save_plot(
    figure_interactions_supplement,
    "figure_3_oof_treeshap_interactions_supplement",
    width = WIDE_FIGURE_WIDTH,
    height = 4.25 * length(SUPPLEMENT_SCHEMES)
  )
}


# ------------------------------------------------------------------------------
# 11. FIGURE 4: TOP-SIX FEATURE DEPENDENCE PLOTS BY SCHEME
# ------------------------------------------------------------------------------

for (scheme_number in seq_along(SCHEMES)) {
  scheme_name <- SCHEMES[scheme_number]
  scheme_id <- scheme_name
  
  top_features <- global_importance[
    scheme == scheme_id & importance_rank <= TOP_N_FEATURES
  ][order(importance_rank), feature]
  
  if (length(top_features) < TOP_N_FEATURES) {
    stop("Fewer than six ranked features are available for ", scheme_name, ".")
  }
  
  dependence_plots <- lapply(seq_along(top_features), function(i) {
    feature_name <- top_features[i]
    partner_name <- strongest_partner(feature_name, scheme_name)
    if (is.na(partner_name)) {
      stop(
        "No TreeSHAP interaction partner is available for ",
        feature_name,
        " in ",
        scheme_name,
        "."
      )
    }
    
    make_feature_dependence_plot(
      scheme_name = scheme_name,
      feature_name = feature_name,
      partner_name = partner_name
    )
  })
  
  figure_dependence <- wrap_plots(
    dependence_plots,
    ncol = 2L
  ) +
    plot_annotation(
      tag_levels = "A",
      tag_prefix = "(",
      tag_suffix = ")"
    )
  
  save_plot(
    figure_dependence,
    paste0(
      "figure_4_dependence_top6_",
      safe_id(scheme_name)
    ),
    width = WIDE_FIGURE_WIDTH,
    height = 9.5
  )
}


# ------------------------------------------------------------------------------
# 12. PROVENANCE
# ------------------------------------------------------------------------------

panel_scheme_caption <- function(schemes) {
  sample_counts <- vapply(schemes, function(scheme_name) {
    min(
      nrow(unit_values[scheme == scheme_name]),
      SUMMARY_MAX_UNITS_PER_SCHEME
    )
  }, integer(1))
  paste(
    sprintf(
      "(%s) %s (n = %s)",
      LETTERS[seq_along(schemes)],
      scheme_label(schemes),
      sample_counts
    ),
    collapse = "; "
  )
}

interaction_panel_caption <- function(schemes) {
  paste(
    vapply(seq_along(schemes), function(i) {
      letters <- LETTERS[2L * i - c(1L, 0L)]
      scheme_id <- schemes[i]
      strongest <- interaction_importance[
        scheme == scheme_id
      ][order(interaction_rank)][1L]
      pair_name <- paste(
        feature_label(strongest$feature_1),
        feature_label(strongest$feature_2),
        sep = " x "
      )
      n_sampled <- source_interaction_counts[
        scheme == scheme_id, n_sampled_observations
      ][1L]
      paste0(
        "(", letters[1L], ") ", scheme_label(schemes[i]),
        " five strongest interaction pairs (n = ", n_sampled,
        " source observations); (", letters[2L],
        ") signed dependence of the strongest pair (", pair_name, ")"
      )
    }, character(1)),
    collapse = "; "
  )
}

figure_manifest <- data.table(
  figure = c(
    "figure_1_oof_treeshap_summary_joint_expanding",
    "figure_2_relative_stage_attribution_joint_expanding",
    "figure_3_oof_treeshap_interactions_joint_expanding"
  ),
  group = "main",
  description = c(
    "Out-of-fold TreeSHAP summary",
    "Relative attribution by observed BBCH stage",
    "Top-five exact interactions and strongest-pair dependence"
  ),
  caption = c(
    paste0(
      "Out-of-fold TreeSHAP summary. ",
      panel_scheme_caption(MAIN_SCHEMES),
      ". Each point represents a held-out environment-by-observed-stage unit; ",
      "colour shows the predictor value."
    ),
    paste0(
      "Relative absolute TreeSHAP attribution across observed BBCH stages. ",
      "Panels show ", paste(scheme_label(MAIN_SCHEMES), collapse = " and "),
      ". Each cell is the mean absolute attribution for a feature at a stage, ",
      "expressed as a share of total attribution at that stage."
    ),
    paste0(
      "Out-of-fold TreeSHAP interaction importance. ",
      interaction_panel_caption(MAIN_SCHEMES),
      ". Bars rank mean absolute pair contributions on the balanced interaction ",
      "sample; dependence points show signed pair contributions, with colour ",
      "indicating the other feature in each pair."
    )
  )
)

if (length(SUPPLEMENT_SCHEMES)) {
  figure_manifest <- rbind(
    figure_manifest,
    data.table(
      figure = c(
        "figure_1_oof_treeshap_summary_supplement",
        "figure_2_relative_stage_attribution_supplement",
        "figure_3_oof_treeshap_interactions_supplement"
      ),
      group = "supplement",
      description = c(
        "Out-of-fold TreeSHAP summary for other schemes",
        "Relative attribution by observed BBCH stage for other schemes",
        "Top-five exact interactions for other schemes"
      ),
      caption = c(
        paste0(
          "Out-of-fold TreeSHAP summary for supplemental validation schemes. ",
          panel_scheme_caption(SUPPLEMENT_SCHEMES),
          ". Each point represents a held-out environment-by-observed-stage ",
          "unit; colour shows the predictor value."
        ),
        paste0(
          "Relative absolute TreeSHAP attribution across observed BBCH stages ",
          "for ", paste(scheme_label(SUPPLEMENT_SCHEMES), collapse = ", "),
          ". Each cell is the mean absolute attribution for a feature at a ",
          "stage, expressed as a share of total attribution at that stage."
        ),
        paste0(
          "Out-of-fold TreeSHAP interaction importance for supplemental ",
          "validation schemes. ", interaction_panel_caption(SUPPLEMENT_SCHEMES),
          ". Bars rank mean absolute pair contributions on the balanced ",
          "interaction sample; dependence points show signed pair contributions, ",
          "with colour indicating the other feature in each pair."
        )
      )
    )
  )
}

dependence_manifest <- data.table(
  figure = paste0("figure_4_dependence_top6_", safe_id(SCHEMES)),
  group = ifelse(SCHEMES %in% MAIN_SCHEMES, "main", "supplement"),
  description = paste0(
    "Top-six feature dependence plots for ", scheme_label(SCHEMES)
  ),
  caption = paste0(
    "Out-of-fold TreeSHAP dependence for the six features with the highest ",
    "mean absolute attribution across environment-by-stage units under ",
    scheme_label(SCHEMES), " validation. Each panel shows the ",
    "predictor value against its signed SHAP contribution; point colour ",
    "shows that predictor's strongest exact interaction partner. The black ",
    "curve is a GAM smooth."
  )
)
figure_manifest <- rbind(figure_manifest, dependence_manifest)
write_table(figure_manifest, "figure_manifest.csv")

caption_lines <- unlist(lapply(seq_len(nrow(figure_manifest)), function(i) {
  c(
    paste0("## ", figure_manifest$figure[i]),
    "",
    figure_manifest$caption[i],
    ""
  )
}))
writeLines(
  c("# Draft figure captions", "", caption_lines),
  file.path(OUTPUT_DIR, "figure_captions_draft.md")
)

plot_configuration <- list(
  source_shap_results_dir = normalizePath(SHAP_RESULTS_DIR),
  source_shap_run_id = SHAP_RUN_ID,
  source_analysis_configuration = analysis_configuration_path,
  source_oof_shap = row_shap_path,
  source_oof_interactions = interaction_path,
  output_dir = normalizePath(OUTPUT_DIR),
  main_schemes = MAIN_SCHEMES,
  supplement_schemes = SUPPLEMENT_SCHEMES,
  schemes = SCHEMES,
  feature_set = analysis_configuration$feature_set,
  features = FEATURES,
  source_scheme_counts = source_scheme_counts,
  source_interaction_counts = source_interaction_counts,
  primary_weighting =
    "macro mean over environment x exact observed BBCH units",
  summary_layout = if (HAS_SHAPFORXGBOOST) {
    "SHAPforxgboost::shap.plot.summary"
  } else {
    "equivalent custom ggplot layout"
  },
  summary_max_units_per_scheme = SUMMARY_MAX_UNITS_PER_SCHEME,
  dependence_max_units_per_scheme = DEPENDENCE_MAX_UNITS_PER_SCHEME,
  top_n_features = TOP_N_FEATURES,
  top_n_interactions = TOP_N_INTERACTIONS,
  interaction_definition =
    "two times one symmetric off-diagonal TreeSHAP interaction cell",
  created_at_utc = format(Sys.time(), tz = "UTC")
)
atomic_save_rds(
  plot_configuration,
  file.path(OUTPUT_DIR, "plot_configuration.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "session_info.txt")
)

message("Main and available supplemental SHAP plots finished.")
message("Source SHAP results: ", SHAP_RESULTS_DIR)
message("Plots written to: ", PLOT_DIR)
message("Tables written to: ", TABLE_DIR)
message(
  "SHAP values are model attributions, not causal effects; stage-specific ",
  "patterns are conditional on the observed BBCH stage."
)

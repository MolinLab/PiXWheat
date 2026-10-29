# ==============================================================================
# PAPER FIGURES AND DIAGNOSTICS FOR THE MULTI-SCHEME GENERALISATION BENCHMARK
# ==============================================================================
#
# This script consumes the outputs of multischeme_generalisation_benchmark.R.
# It does not refit or retune any model.
#
# Main paper outputs
# ------------------
#   * figure_cv_schemes_paper.{pdf,png,eps}
#       Data-specific visualization of all five outer validation schemes.
#   * figure_fold_performance_selected_learner.{pdf,png,eps}
#       Outer-fold MAE and RMSE distributions for the selected learner.
#   * figure_fold_performance_vs_gdd.{pdf,png,eps}
#       The same distributions with the calibrated GDD baseline.
#
# Diagnostic outputs
# ------------------
#   * fold_performance_quantiles.csv
#   * stage_specific_mae_across_cv_schemes.{pdf,png,eps}
#   * stage_specific_rmse_across_cv_schemes.{pdf,png,eps}
#   * stage_of_interest_diagnostics_*.csv
#   * tuning_budget_diagnostics_*.csv
#   * tuning_best_parameter_summary.csv
#   * figure_tuning_convergence_diagnostic.{pdf,png,eps}
#   * tuning_budget_interpretation.txt
#
# Boxplot definition
# ------------------
# The box spans the 25th--75th percentiles of OUTER-FOLD scores, the central
# line is the median, whiskers extend to 1.5 x IQR, and points are individual
# outer folds. A white diamond marks the pooled out-of-fold score reported by
# the benchmark. Fold distributions are descriptive and do not replace the
# scheme-specific cluster-bootstrap confidence intervals.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. USER CONFIGURATION
# ------------------------------------------------------------------------------

PROJECT_DIR <- Sys.getenv(
  "PHENOLOGY_PROJECT_DIR",
  unset = "/data/analysis3/milekj/UseCase_WW02_2024"
)

SELECTED_LEARNER <- tolower(trimws(Sys.getenv(
  "PHENOLOGY_GENERALISATION_LEARNER",
  unset = "xgboost"
)))

RUN_ID <- Sys.getenv(
  "PHENOLOGY_GENERALISATION_RUN_ID",
  unset = paste0("multischeme_", SELECTED_LEARNER, "_paper_cum")
)

RESULT_DIR <- Sys.getenv(
  "PHENOLOGY_GENERALISATION_RESULT_DIR",
  unset = file.path(PROJECT_DIR, "benchmark_results", RUN_ID)
)

FIGURE_DIR <- Sys.getenv(
  "PHENOLOGY_GENERALISATION_FIGURE_DIR",
  unset = file.path(RESULT_DIR, "new_figures")
)

BASELINE_LEARNER <- "gdd_threshold"
STAGE_OF_INTEREST <- as.numeric(Sys.getenv(
  "PHENOLOGY_STAGE_OF_INTEREST",
  unset = "75"
))

# A fold is called plateaued when the best inner-CV MAE improves by less than
# this fraction during the final quarter of the MBO evaluations. This is a
# transparent heuristic diagnostic, not a formal convergence test.
MBO_PLATEAU_RELATIVE_THRESHOLD <- 0.005
MBO_BOUNDARY_FRACTION <- 0.05
MBO_SYSTEMATIC_BOUNDARY_PROPORTION <- 0.50

PNG_DPI <- 600

parse_boolean <- function(x, default = FALSE) {
  if (length(x) == 0L || is.na(x) || !nzchar(trimws(x))) return(default)
  normalized <- tolower(trimws(x))
  if (normalized %in% c("1", "true", "yes", "y")) return(TRUE)
  if (normalized %in% c("0", "false", "no", "n")) return(FALSE)
  stop("Cannot parse logical value: ", x)
}

RUN_TUNING_AUDIT <- parse_boolean(
  Sys.getenv("PHENOLOGY_RUN_TUNING_AUDIT", unset = "TRUE"),
  default = TRUE
)

RUN_STAGE_DIAGNOSTICS <- parse_boolean(
  Sys.getenv("PHENOLOGY_RUN_STAGE_DIAGNOSTICS", unset = "TRUE"),
  default = TRUE
)

# ------------------------------------------------------------------------------
# 2. PACKAGES, LABELS AND GRAPHICAL CONSTANTS
# ------------------------------------------------------------------------------

required_packages <- c("data.table", "ggplot2", "patchwork", "scales", "systemfonts")
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

if (!capabilities("cairo")) {
  stop("Cairo graphics support is required for Liberation Sans PDF and EPS export.")
}
if (!any(tolower(trimws(systemfonts::system_fonts()$family)) == "liberation sans")) {
  stop("Liberation Sans is not available to R on this system.")
}

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

scheme_levels <- c(
  "spatial_block",
  "leave_one_year_out",
  "expanding_origin",
  "environmental_zone",
  "joint_space_time"
)

scheme_labels <- c(
  spatial_block = "Spatial\nblocks",
  leave_one_year_out = "Leave-one-\nseason-year-out",
  expanding_origin = "Expanding-\norigin",
  environmental_zone = "Environmental\nzones",
  joint_space_time = "Joint\nspace–time"
)

model_labels <- c(
  random_forest = "Random forest",
  xgboost = "XGBoost",
  gdd_threshold = "Calibrated thermal time (cumDTU)"
)

model_colours <- c(
  random_forest = "#009E73",
  xgboost = "#0072B2",
  gdd_threshold = "#777777"
)

partition_colours <- c(
  Training = "#D9D9D9",
  Test = "#0072B2",
  `Not used` = "#FFFFFF"
)

theme_paper <- function(base_size = 9.5) {
  theme_bw(base_size = base_size, base_family = "Liberation Sans") +
    theme(
      text = element_text(family = "Liberation Sans"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.20, colour = "grey90"),
      axis.title = element_text(colour = "black"),
      axis.text = element_text(colour = "black"),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      plot.tag = element_text(family = "Liberation Sans", face = "bold", size = rel(1.15))
    )
}

factor_scheme <- function(x) {
  factor(x, levels = scheme_levels, labels = unname(scheme_labels[scheme_levels]))
}

model_display <- function(x) {
  answer <- unname(model_labels[x])
  answer[is.na(answer)] <- x[is.na(answer)]
  answer
}

thin_breaks <- function(values, maximum = 9L) {
  values <- unique(as.character(values))
  if (length(values) <= maximum) return(values)
  positions <- unique(as.integer(round(seq(1, length(values), length.out = maximum))))
  values[positions]
}

wrap_text <- function(x, width = 80L) {
  vapply(
    as.character(x),
    function(value) paste(strwrap(value, width = width), collapse = "\n"),
    character(1),
    USE.NAMES = FALSE
  )
}

save_publication_plot <- function(plot, stem, width, height) {
  pdf_path <- file.path(FIGURE_DIR, paste0(stem, ".pdf"))
  png_path <- file.path(FIGURE_DIR, paste0(stem, ".png"))
  eps_path  <- file.path(FIGURE_DIR, paste0(stem, ".eps"))
  
  ggsave(
    filename = pdf_path,
    plot = plot,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
  ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = PNG_DPI,
    bg = "white",
    limitsize = FALSE
  )
  # EPS for journal submission
  ggsave(
    eps_path, plot,
    device = grDevices::cairo_ps,
    width = width, height = height,
    units = "in",
    bg = "white",
    limitsize = FALSE,
    onefile = FALSE
  )
  
  invisible(c(pdf = pdf_path, png = png_path, eps = eps_path))
}

read_required <- function(filename, select = NULL) {
  path <- file.path(RESULT_DIR, filename)
  if (!file.exists(path)) stop("Required benchmark output is missing: ", path)
  fread(path, select = select, na.strings = c("NA", "NaN", ""))
}

assert_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(
      label, " is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }
}

# ------------------------------------------------------------------------------
# 3. LOAD AND VALIDATE THE COMPLETED BENCHMARK
# ------------------------------------------------------------------------------

if (!dir.exists(RESULT_DIR)) {
  stop("Generalisation result directory does not exist: ", RESULT_DIR)
}

manifest <- read_required("outer_fold_manifest.csv")
split_specs <- read_required("outer_split_specs.csv")
split_sizes <- read_required("outer_split_sizes.csv")
fold_metrics <- read_required("outer_fold_metrics.csv")
performance <- read_required("generalisation_performance_by_scheme.csv")
stage_metrics <- read_required("stage_specific_metrics.csv")

# Keep block identifiers categorical even when fread infers all observed
# values as integers. This avoids type mismatches when joining display labels.
split_specs[, holdout_spatial_fold := as.character(holdout_spatial_fold)]
split_specs[, holdout_temporal_block := as.character(holdout_temporal_block)]
manifest[, spatial_block := as.character(spatial_block)]
manifest[, joint_temporal_block := as.character(joint_temporal_block)]

assert_columns(
  manifest,
  c(
    "obs_id", "site_id", "season_year", "env_zone", "lon", "lat",
    "spatial_block", "joint_temporal_block"
  ),
  "outer_fold_manifest.csv"
)
assert_columns(
  split_specs,
  c(
    "split_id", "scheme", "fold_id", "holdout_spatial_fold",
    "holdout_year", "holdout_zone", "holdout_temporal_block"
  ),
  "outer_split_specs.csv"
)
assert_columns(
  split_sizes,
  c("split_id", "scheme", "fold_id", "n_train", "n_test", "n_excluded"),
  "outer_split_sizes.csv"
)
assert_columns(
  fold_metrics,
  c("learner_id", "scheme", "split_id", "fold_id", "mae", "rmse"),
  "outer_fold_metrics.csv"
)
assert_columns(
  performance,
  c("learner_id", "scheme", "mae", "rmse"),
  "generalisation_performance_by_scheme.csv"
)
assert_columns(
  stage_metrics,
  c("learner_id", "scheme", "BBCH_stage", "n", "mae", "rmse", "bias"),
  "stage_specific_metrics.csv"
)

available_schemes <- intersect(scheme_levels, unique(split_specs$scheme))
if (!length(available_schemes)) stop("No recognized validation schemes found.")

if (!SELECTED_LEARNER %in% fold_metrics$learner_id) {
  stop(
    "Selected learner '", SELECTED_LEARNER,
    "' is absent from outer_fold_metrics.csv. Available learners: ",
    paste(sort(unique(fold_metrics$learner_id)), collapse = ", ")
  )
}

job_status_path <- file.path(RESULT_DIR, "job_status.csv")
if (file.exists(job_status_path)) {
  job_status <- fread(job_status_path, na.strings = c("NA", "NaN", ""))
  if ("status" %in% names(job_status) && any(job_status$status != "ok")) {
    stop(
      "At least one outer job did not finish successfully. Inspect ",
      job_status_path, " before making paper figures."
    )
  }
}

# ------------------------------------------------------------------------------
# 4. FIGURE: DATA-SPECIFIC OUTER RESAMPLING SCHEMES
# ------------------------------------------------------------------------------

site_qc <- manifest[, .(
  n_coordinates = uniqueN(paste(round(lon, 6), round(lat, 6))),
  n_spatial_blocks = uniqueN(spatial_block),
  n_environmental_zones = uniqueN(env_zone)
), by = site_id]

if (any(site_qc$n_coordinates != 1L)) {
  stop("At least one site has multiple coordinate pairs in the fold manifest.")
}
if (any(site_qc$n_spatial_blocks != 1L)) {
  stop("At least one site was assigned to multiple spatial blocks.")
}
if (any(site_qc$n_environmental_zones != 1L)) {
  stop("At least one site was assigned to multiple environmental zones.")
}

sites <- manifest[, .(
  lon = lon[1L],
  lat = lat[1L],
  spatial_block = as.character(spatial_block[1L]),
  env_zone = as.character(env_zone[1L])
), by = site_id]

map_extent <- sites[, .(
  xmin = min(lon), xmax = max(lon),
  ymin = min(lat), ymax = max(lat)
)]
x_padding <- max(0.15, 0.05 * (map_extent$xmax - map_extent$xmin))
y_padding <- max(0.10, 0.05 * (map_extent$ymax - map_extent$ymin))

country_map <- NULL
if (requireNamespace("maps", quietly = TRUE)) {
  country_map <- ggplot2::map_data("world")
  country_map <- country_map[country_map$region %in% c("Germany", "Austria"), ]
}

geographic_base <- function() {
  p <- ggplot()
  if (!is.null(country_map) && nrow(country_map)) {
    p <- p + geom_polygon(
      data = country_map,
      aes(x = long, y = lat, group = group),
      inherit.aes = FALSE,
      fill = "grey97",
      colour = "grey70",
      linewidth = 0.25
    )
  }
  p
}

map_theme <- theme_paper() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.key.width = grid::unit(0.70, "lines"),
    legend.key.height = grid::unit(0.70, "lines"),
    legend.text = element_text(size = 6.5),
    legend.title = element_text(size = 7.5, face = "bold"),
    legend.spacing.x = grid::unit(0.15, "lines"),
    legend.margin = margin(t = 1, r = 1, b = 1, l = 1)
  )

p_spatial <- geographic_base() +
  geom_point(
    data = sites,
    aes(x = lon, y = lat, colour = factor(spatial_block)),
    size = 1.45,
    alpha = 0.80
  ) +
  scale_colour_viridis_d(
    option = "D",
    end = 0.90,
    name = "Spatial block",
    guide = guide_legend(
      ncol = 1,
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  coord_quickmap(
    xlim = c(map_extent$xmin - x_padding, map_extent$xmax + x_padding),
    ylim = c(map_extent$ymin - y_padding, map_extent$ymax + y_padding),
    expand = FALSE
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  map_theme

p_zone <- geographic_base() +
  geom_point(
    data = sites,
    aes(x = lon, y = lat, colour = factor(env_zone)),
    size = 1.45,
    alpha = 0.80
  ) +
  scale_colour_viridis_d(
    option = "C",
    end = 0.92,
    name = "Environmental zone",
    labels = function(x) wrap_text(x, width = 14L),
    guide = guide_legend(
      ncol = 1,
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  coord_quickmap(
    xlim = c(map_extent$xmin - x_padding, map_extent$xmax + x_padding),
    ylim = c(map_extent$ymin - y_padding, map_extent$ymax + y_padding),
    expand = FALSE
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  map_theme

build_year_role_matrix <- function(scheme_id) {
  specs <- copy(split_specs[scheme == scheme_id])
  if (!nrow(specs)) return(NULL)
  
  years <- sort(unique(manifest$season_year))
  grid <- CJ(split_id = specs$split_id, season_year = years, unique = TRUE)
  grid <- merge(
    grid,
    specs[, .(split_id, fold_id, holdout_year)],
    by = "split_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  if (identical(scheme_id, "leave_one_year_out")) {
    grid[, partition := fifelse(season_year == holdout_year, "Test", "Training")]
  } else if (identical(scheme_id, "expanding_origin")) {
    grid[, partition := fcase(
      season_year < holdout_year, "Training",
      season_year == holdout_year, "Test",
      default = "Not used"
    )]
  } else {
    stop("Unsupported temporal scheme: ", scheme_id)
  }
  
  split_order <- specs[order(holdout_year), as.character(holdout_year)]
  grid[, split_label := factor(
    as.character(holdout_year),
    levels = rev(split_order)
  )]
  grid[, partition := factor(
    partition,
    levels = c("Training", "Test", "Not used")
  )]
  grid
}

plot_temporal_matrix <- function(role_data) {
  y_levels <- levels(role_data$split_label)
  ggplot(role_data, aes(x = season_year, y = split_label, fill = partition)) +
    geom_tile(colour = "grey85", linewidth = 0.10) +
    scale_fill_manual(
      values = partition_colours,
      breaks = c("Training", "Test", "Not used"),
      drop = FALSE,
      name = "Partition"
    ) +
    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 8),
      expand = expansion(mult = 0)
    ) +
    scale_y_discrete(
      breaks = thin_breaks(y_levels, maximum = 8L),
      drop = FALSE
    ) +
    labs(
      x = "Season year",
      y = "Outer test year"
    ) +
    theme_paper() +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      axis.text.y = element_text(size = 7)
    )
}

loyo_roles <- build_year_role_matrix("leave_one_year_out")
expanding_roles <- build_year_role_matrix("expanding_origin")

if (is.null(loyo_roles) || is.null(expanding_roles)) {
  stop(
    "Both leave_one_year_out and expanding_origin outputs are required for ",
    "the complete validation-scheme figure."
  )
}

p_loyo <- plot_temporal_matrix(loyo_roles)
p_expanding <- plot_temporal_matrix(expanding_roles)

temporal_row <- (
  (p_loyo + theme(legend.position = "none")| p_expanding + theme(legend.position = "right"))) 


joint_specs <- merge(
  split_specs[scheme == "joint_space_time"],
  split_sizes[
    scheme == "joint_space_time",
    .(split_id, n_train, n_test, n_excluded)
  ],
  by = "split_id",
  all.x = TRUE,
  sort = FALSE
)

if (!nrow(joint_specs)) {
  stop("joint_space_time split specifications are missing.")
}

time_block_labels <- manifest[
  !is.na(joint_temporal_block),
  .(
    year_min = min(season_year),
    year_max = max(season_year)
  ),
  by = .(holdout_temporal_block = as.character(joint_temporal_block))
]
time_block_labels[, temporal_label := ifelse(
  year_min == year_max,
  paste0("T", holdout_temporal_block, "\n", year_min),
  paste0("T", holdout_temporal_block, "\n", year_min, "–", year_max)
)]

joint_specs <- merge(
  joint_specs,
  time_block_labels,
  by = "holdout_temporal_block",
  all.x = TRUE,
  sort = FALSE
)

spatial_order <- sort(unique(as.character(joint_specs$holdout_spatial_fold)))
temporal_order <- time_block_labels[
  order(year_min, year_max),
  temporal_label
]
joint_specs[, spatial_label := factor(
  paste0("S", holdout_spatial_fold),
  levels = rev(paste0("S", spatial_order))
)]
joint_specs[, temporal_label := factor(
  temporal_label,
  levels = temporal_order
)]

p_joint <- ggplot(
  joint_specs,
  aes(x = temporal_label, y = spatial_label, fill = n_test)
) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(
    aes(label = scales::comma(n_test)),
    size = 2.5,
    family = "Liberation Sans",
    fontface = "bold",
    colour = "black"
  ) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#6BAED6",
    name = "Test rows",
    labels = scales::label_comma()
  ) +
  labs(
    x = "Temporal block (season-year range)",
    y = "Spatial block"
  ) +
  theme_paper() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(size = 7.5)
  )

cv_figure <- (
  (p_spatial | p_zone) /
    temporal_row /
    p_joint
) +
  # Maps have a fixed aspect ratio; a taller first row adds unused space.
  plot_layout(heights = c(1.15, 1.1, 1)) +
  plot_annotation(
    tag_levels = "A", tag_prefix = "(", tag_suffix = ")",
    theme = theme(
      plot.tag = element_text(family = "Liberation Sans", face = "bold", size = 11),
      plot.margin = margin(t = 1, r = 2, b = 2, l = 2, unit = "pt")
    )
  )

save_publication_plot(
  cv_figure,
  "figure_cv_schemes_paper",
  width = 7.4,
  height = 8.9
)

# ------------------------------------------------------------------------------
# 5. FIGURES: OUTER-FOLD PERFORMANCE DISTRIBUTIONS
# ------------------------------------------------------------------------------

models_to_plot <- intersect(
  c(SELECTED_LEARNER, BASELINE_LEARNER),
  unique(fold_metrics$learner_id)
)

fold_metrics_plot <- fold_metrics[
  learner_id %in% models_to_plot & is.finite(mae) & is.finite(rmse)
]
fold_metrics_plot[, scheme_display := factor_scheme(scheme)]
fold_metrics_plot[, model_display := factor(
  model_display(learner_id),
  levels = model_display(models_to_plot)
)]

fold_long <- melt(
  fold_metrics_plot,
  id.vars = c(
    "learner_id", "model_display", "scheme", "scheme_display",
    "split_id", "fold_id", "n_train", "n_test"
  ),
  measure.vars = c("mae", "rmse"),
  variable.name = "metric",
  value.name = "error"
)
fold_long[, metric := factor(
  metric,
  levels = c("mae", "rmse"),
  labels = c("MAE", "RMSE")
)]

performance_plot <- performance[learner_id %in% models_to_plot]
performance_plot[, scheme_display := factor_scheme(scheme)]
performance_plot[, model_display := factor(
  model_display(learner_id),
  levels = model_display(models_to_plot)
)]
performance_long <- melt(
  performance_plot,
  id.vars = c("learner_id", "model_display", "scheme", "scheme_display"),
  measure.vars = c("mae", "rmse"),
  variable.name = "metric",
  value.name = "error"
)
performance_long[, metric := factor(
  metric,
  levels = c("mae", "rmse"),
  labels = c("MAE", "RMSE")
)]

fold_quantiles <- fold_long[, .(
  n_outer_folds = .N,
  minimum = min(error),
  q25 = as.numeric(quantile(error, 0.25)),
  median = median(error),
  q75 = as.numeric(quantile(error, 0.75)),
  maximum = max(error),
  mean = mean(error),
  sd = sd(error)
), by = .(learner_id, scheme, metric)]

pooled_values <- performance_long[, .(
  learner_id,
  scheme,
  metric,
  pooled_oof_error = error
)]
fold_quantiles <- merge(
  fold_quantiles,
  pooled_values,
  by = c("learner_id", "scheme", "metric"),
  all.x = TRUE
)
setorder(fold_quantiles, learner_id, scheme, metric)
fwrite(
  fold_quantiles,
  file.path(FIGURE_DIR, "fold_performance_quantiles.csv"),
  na = "NA"
)

selected_colour <- unname(model_colours[SELECTED_LEARNER])
if (is.na(selected_colour)) selected_colour <- "#0072B2"

selected_fold_long <- fold_long[learner_id == SELECTED_LEARNER]
selected_performance_long <- performance_long[learner_id == SELECTED_LEARNER]

p_selected_performance <- ggplot(
  selected_fold_long,
  aes(x = scheme_display, y = error)
) +
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    fill = scales::alpha(selected_colour, 0.25),
    colour = selected_colour,
    linewidth = 0.50
  ) +
  geom_jitter(
    width = 0.10,
    height = 0,
    size = 1.25,
    alpha = 0.50,
    colour = selected_colour
  ) +
  geom_point(
    data = selected_performance_long,
    aes(x = scheme_display, y = error),
    inherit.aes = FALSE,
    shape = 23,
    size = 2.6,
    stroke = 0.65,
    fill = "white",
    colour = "black"
  ) +
  facet_wrap(vars(metric), scales = "free_y", nrow = 1) +
  labs(
    x = NULL,
    y = "Prediction error (BBCH stages)"
  ) +
  theme_paper(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8, lineheight = 0.90),
    legend.position = "none"
  )

save_publication_plot(
  p_selected_performance,
  "figure_fold_performance_selected_learner",
  width = 7.4,
  height = 4.7
)

if (BASELINE_LEARNER %in% models_to_plot) {
  comparison_colours <- model_colours[models_to_plot]
  names(comparison_colours) <- model_display(models_to_plot)
  dodge <- position_dodge(width = 0.72)
  
  p_performance_comparison <- ggplot(
    fold_long,
    aes(x = scheme_display, y = error, fill = model_display)
  ) +
    geom_boxplot(
      position = dodge,
      width = 0.62,
      outlier.shape = NA,
      linewidth = 0.42,
      alpha = 0.75
    ) +
    geom_point(
      aes(colour = model_display),
      position = position_jitterdodge(
        jitter.width = 0.08,
        jitter.height = 0,
        dodge.width = 0.72
      ),
      size = 0.95,
      alpha = 0.42,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = comparison_colours, name = "Model") +
    scale_colour_manual(values = comparison_colours, name = "Model") +
    facet_wrap(vars(metric), scales = "free_y", nrow = 1) +
    labs(
      x = NULL,
      y = "Prediction error (BBCH stages)"
    ) +
    theme_paper(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0.9, size = 8, lineheight = 0.90),
      legend.position = "bottom"
    )
  
  save_publication_plot(
    p_performance_comparison,
    "figure_fold_performance_vs_gdd",
    width = 7.4,
    height = 4.8
  )
}

# ------------------------------------------------------------------------------
# 6. STAGE-SPECIFIC PERFORMANCE AND BBCH-50 DIAGNOSTICS
# ------------------------------------------------------------------------------
CV_SCHEMES <- c(
  "spatial_block",
  "leave_one_year_out",
  "environmental_zone",
  "joint_space_time",
  "expanding_origin"
)

stage_selected <- copy(
  stage_metrics[
    learner_id == SELECTED_LEARNER &
      scheme %in% CV_SCHEMES
  ]
)

# Robust conversion, including values such as "BBCH 31"
stage_selected[, BBCH_stage_numeric := as.numeric(
  gsub("[^0-9.]+", "", as.character(BBCH_stage))
)]

if (anyNA(stage_selected$BBCH_stage_numeric)) {
  stop("Some BBCH stages could not be converted to numeric values.")
}

stage_levels <- sort(unique(stage_selected$BBCH_stage_numeric))

# Ordered plotting position
stage_selected[, stage_index := match(
  BBCH_stage_numeric,
  stage_levels
)]

setorder(
  stage_selected,
  scheme,
  BBCH_stage_numeric
)

stage_long <- melt(
  stage_selected,
  id.vars = c(
    "scheme",
    "BBCH_stage_numeric",
    "stage_index",
    "learner_id"
  ),
  measure.vars = c("mae", "rmse"),
  variable.name = "metric",
  value.name = "error"
)

stage_summary <- stage_long[, .(
  median_error = median(error, na.rm = TRUE),
  minimum_error = min(error, na.rm = TRUE),
  maximum_error = max(error, na.rm = TRUE),
  n_schemes = uniqueN(scheme)
), by = .(
  metric,
  BBCH_stage_numeric,
  stage_index
)]
make_stage_scheme_summary_plot <- function(metric_name, metric_label) {
  
  plot_data <- stage_long[metric == metric_name]
  summary_data <- stage_summary[metric == metric_name]
  
  ggplot() +
    geom_ribbon(
      data = summary_data,
      aes(
        x = stage_index,
        ymin = minimum_error,
        ymax = maximum_error
      ),
      fill = "#56B4E9",
      alpha = 0.20
    ) +
    geom_line(
      data = plot_data,
      aes(
        x = stage_index,
        y = error,
        group = scheme
      ),
      colour = "grey55",
      linewidth = 0.45,
      alpha = 0.55
    ) +
    geom_point(
      data = plot_data,
      aes(
        x = stage_index,
        y = error,
        group = scheme
      ),
      colour = "grey55",
      size = 1.3,
      alpha = 0.55
    ) +
    geom_line(
      data = summary_data,
      aes(
        x = stage_index,
        y = median_error,
        group = 1
      ),
      colour = "#0072B2",
      linewidth = 1.05
    ) +
    geom_point(
      data = summary_data,
      aes(
        x = stage_index,
        y = median_error
      ),
      colour = "#0072B2",
      size = 2.5
    ) +
    scale_x_continuous(
      breaks = seq_along(stage_levels),
      labels = stage_levels,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
      x = "Observed BBCH stage",
      y = paste0(metric_label, " (BBCH units)")
    ) +
    theme_paper(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        size = 10,
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 11)
    )
}

p_stage_mae <- make_stage_scheme_summary_plot("mae", "MAE")
p_stage_rmse <- make_stage_scheme_summary_plot("rmse", "RMSE")

save_publication_plot(
  p_stage_mae,
  "stage_specific_mae_across_cv_schemes",
  width = 7.4,
  height = 6.8
)

save_publication_plot(
  p_stage_rmse,
  "stage_specific_rmse_across_cv_schemes",
  width = 7.4,
  height = 6.8
)
##########################################################################

stage_interest_from_summary <- stage_selected[BBCH_stage == STAGE_OF_INTEREST]
fwrite(
  stage_interest_from_summary,
  file.path(
    FIGURE_DIR,
    paste0("stage_", STAGE_OF_INTEREST, "_metrics_by_scheme.csv")
  ),
  na = "NA"
)


if (RUN_STAGE_DIAGNOSTICS) {
  oof_path <- file.path(RESULT_DIR, "oof_predictions.csv")
  if (!file.exists(oof_path)) {
    warning("Stage diagnostics skipped because oof_predictions.csv is missing.")
  } else {
    oof_columns <- c(
      "learner_id", "scheme", "source_row_id", "env_id", "site_id", "season_year",
      "env_zone", "truth", "response", "error", "absolute_error",
      "squared_error"
    )
    oof <- fread(
      oof_path,
      select = oof_columns,
      na.strings = c("NA", "NaN", "")
    )
    oof <- oof[
      learner_id == SELECTED_LEARNER & truth == STAGE_OF_INTEREST
    ]
    
    if (!nrow(oof)) {
      warning(
        "No selected-learner predictions found for BBCH ",
        STAGE_OF_INTEREST, "."
      )
    } else {
      summarize_stage_subset <- function(data, grouping_columns) {
        data[, .(
          n = .N,
          n_environments = uniqueN(env_id),
          n_sites = uniqueN(site_id),
          mae = mean(absolute_error),
          rmse = sqrt(mean(squared_error)),
          bias = mean(error),
          median_prediction = median(response),
          prediction_q25 = as.numeric(quantile(response, 0.25)),
          prediction_q75 = as.numeric(quantile(response, 0.75)),
          proportion_underpredicted = mean(error < 0),
          proportion_overpredicted = mean(error > 0)
        ), by = grouping_columns]
      }
      
      stage_by_scheme <- summarize_stage_subset(oof, "scheme")
      stage_by_scheme_zone <- summarize_stage_subset(
        oof,
        c("scheme", "env_zone")
      )
      stage_by_scheme_year <- summarize_stage_subset(
        oof,
        c("scheme", "season_year")
      )
      
      # If the original input and a source/provenance column are available,
      # attach them through the source-row identifier saved by the benchmark.
      # This helps distinguish a genuinely difficult stage from a source-specific
      # BBCH coding or observation-protocol issue.
      stage_by_scheme_source <- NULL
      configuration_path <- file.path(RESULT_DIR, "benchmark_configuration.rds")
      if (file.exists(configuration_path)) {
        configuration <- readRDS(configuration_path)
        source_input <- configuration$input_file
        if (
          length(source_input) == 1L &&
          !is.na(source_input) &&
          file.exists(source_input)
        ) {
          input_columns <- names(fread(source_input, nrows = 0L))
          source_candidates <- c(
            "Source_name", "source_name", "Source", "source",
            "data_source", "Datenquelle"
          )
          source_column <- source_candidates[
            source_candidates %in% input_columns
          ][1L]
          if (length(source_column) && !is.na(source_column)) {
            provenance <- fread(
              source_input,
              select = source_column,
              na.strings = c("NA", "NaN", "")
            )
            provenance[, source_row_id := .I]
            setnames(provenance, source_column, "source_name")
            oof <- merge(
              oof,
              provenance,
              by = "source_row_id",
              all.x = TRUE,
              sort = FALSE
            )
            stage_by_scheme_source <- summarize_stage_subset(
              oof,
              c("scheme", "source_name")
            )
            rm(provenance)
          }
        }
      }
      
      fwrite(
        stage_by_scheme,
        file.path(
          FIGURE_DIR,
          paste0("stage_", STAGE_OF_INTEREST, "_diagnostics_by_scheme.csv")
        ),
        na = "NA"
      )
      fwrite(
        stage_by_scheme_zone,
        file.path(
          FIGURE_DIR,
          paste0(
            "stage_", STAGE_OF_INTEREST,
            "_diagnostics_by_scheme_and_zone.csv"
          )
        ),
        na = "NA"
      )
      fwrite(
        stage_by_scheme_year,
        file.path(
          FIGURE_DIR,
          paste0(
            "stage_", STAGE_OF_INTEREST,
            "_diagnostics_by_scheme_and_year.csv"
          )
        ),
        na = "NA"
      )
      if (!is.null(stage_by_scheme_source)) {
        fwrite(
          stage_by_scheme_source,
          file.path(
            FIGURE_DIR,
            paste0(
              "stage_", STAGE_OF_INTEREST,
              "_diagnostics_by_scheme_and_source.csv"
            )
          ),
          na = "NA"
        )
      }
    }
    rm(oof)
    invisible(gc())
  }
}

# ------------------------------------------------------------------------------
# 7. LIGHTWEIGHT AUDIT OF THE COMPLETED MBO SEARCHES
# ------------------------------------------------------------------------------

if (RUN_TUNING_AUDIT) {
  tuning_path <- file.path(RESULT_DIR, "tuning_archives_scalar.csv")
  
  if (!file.exists(tuning_path)) {
    warning("Tuning audit skipped because tuning_archives_scalar.csv is missing.")
  } else {
    tuning <- fread(tuning_path, na.strings = c("NA", "NaN", ""))
    assert_columns(
      tuning,
      c("job_id", "learner_id", "scheme", "split_id", "fold_id"),
      "tuning_archives_scalar.csv"
    )
    tuning <- tuning[learner_id == SELECTED_LEARNER]
    
    objective_candidates <- c("regr.mae", "regr.rmse")
    objective_column <- objective_candidates[
      objective_candidates %in% names(tuning)
    ][1L]
    if (is.na(objective_column) || !length(objective_column)) {
      measure_columns <- grep("^regr\\.", names(tuning), value = TRUE)
      if (!length(measure_columns)) {
        stop("Could not identify the tuning objective column in the archive.")
      }
      objective_column <- measure_columns[1L]
    }
    
    tuning[, archive_row_order := .I]
    order_columns <- c(
      "job_id",
      intersect(c("batch_nr", "iteration", "timestamp"), names(tuning)),
      "archive_row_order"
    )
    setorderv(tuning, order_columns)
    tuning[, evaluation := seq_len(.N), by = job_id]
    tuning <- tuning[is.finite(get(objective_column))]
    
    tuning_budget <- tuning[, {
      objective <- get(objective_column)
      evaluations <- seq_along(objective)
      cutoff <- max(1L, floor(0.75 * length(objective)))
      best_full <- min(objective)
      best_first_75 <- min(objective[evaluations <= cutoff])
      denominator <- max(abs(best_first_75), .Machine$double.eps)
      relative_gain <- (best_first_75 - best_full) / denominator
      
      .(
        n_evaluations = length(objective),
        final_best_inner_mae = best_full,
        best_inner_mae_after_first_75_percent = best_first_75,
        relative_gain_final_quarter = relative_gain,
        best_evaluation = evaluations[which.min(objective)][1L],
        best_found_in_final_quarter =
          evaluations[which.min(objective)][1L] > cutoff,
        plateaued = relative_gain <= MBO_PLATEAU_RELATIVE_THRESHOLD
      )
    }, by = .(job_id, learner_id, scheme, split_id, fold_id)]
    
    tuning_budget_summary <- tuning_budget[, .(
      n_outer_folds = .N,
      median_evaluations = median(n_evaluations),
      proportion_plateaued = mean(plateaued),
      proportion_best_found_in_final_quarter = mean(best_found_in_final_quarter),
      median_relative_gain_final_quarter = median(relative_gain_final_quarter),
      q75_relative_gain_final_quarter = as.numeric(quantile(
        relative_gain_final_quarter,
        0.75
      )),
      maximum_relative_gain_final_quarter = max(relative_gain_final_quarter)
    ), by = scheme]
    
    tuning_budget_overall <- tuning_budget[, .(
      scheme = "all_schemes",
      n_outer_folds = .N,
      median_evaluations = median(n_evaluations),
      proportion_plateaued = mean(plateaued),
      proportion_best_found_in_final_quarter = mean(best_found_in_final_quarter),
      median_relative_gain_final_quarter = median(relative_gain_final_quarter),
      q75_relative_gain_final_quarter = as.numeric(quantile(
        relative_gain_final_quarter,
        0.75
      )),
      maximum_relative_gain_final_quarter = max(relative_gain_final_quarter)
    )]
    tuning_budget_summary <- rbindlist(
      list(tuning_budget_summary, tuning_budget_overall),
      use.names = TRUE
    )
    
    fwrite(
      tuning_budget,
      file.path(FIGURE_DIR, "tuning_budget_diagnostics_by_fold.csv"),
      na = "NA"
    )
    fwrite(
      tuning_budget_summary,
      file.path(FIGURE_DIR, "tuning_budget_diagnostics_summary.csv"),
      na = "NA"
    )
    
    # Parameter ranges exactly match multischeme_generalisation_benchmark.R.
    if (identical(SELECTED_LEARNER, "xgboost")) {
      parameter_space <- data.table(
        parameter = c(
          "nrounds", "eta", "gamma", "max_depth", "min_child_weight",
          "subsample", "colsample_bytree", "lambda", "alpha"
        ),
        lower_natural = c(200, 0.01, 0, 1, 0.1, 0.5, 0.5, 1e-4, 1e-4),
        upper_natural = c(2500, 0.3, 5, 10, 30, 1, 1, 20, 20),
        log_scale = c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE),
        integer_log_scale = c(
          TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
        )
      )
    } else {
      parameter_space <- data.table(
        parameter = c("mtry", "min.node.size", "sample.fraction", "max.depth"),
        lower_natural = c(1, 1, 0.5, 0),
        upper_natural = c(8, 100, 1, 20),
        log_scale = c(FALSE, TRUE, FALSE, FALSE),
        integer_log_scale = c(FALSE, TRUE, FALSE, FALSE)
      )
    }
    # paradox archives retain transformed search-space coordinates. Convert
    # log-scaled parameters back to their natural learner scale for reporting,
    # while using the transformed coordinates for boundary detection.
    parameter_space[, lower_archive := lower_natural]
    parameter_space[, upper_archive := upper_natural]
    parameter_space[
      log_scale == TRUE & integer_log_scale == FALSE,
      `:=`(
        lower_archive = log(lower_natural),
        upper_archive = log(upper_natural)
      )
    ]
    parameter_space[
      integer_log_scale == TRUE,
      `:=`(
        lower_archive = log(pmax(lower_natural, 0.5)),
        upper_archive = log(upper_natural + 1)
      )
    ]
    
    finite_tuning <- tuning[is.finite(get(objective_column))]
    best_row_indices <- finite_tuning[, .I[which.min(get(objective_column))][1L], by = job_id]$V1
    best_configurations <- finite_tuning[best_row_indices]
    available_parameters <- intersect(
      parameter_space$parameter,
      names(best_configurations)
    )
    
    if (length(available_parameters)) {
      best_parameters <- melt(
        best_configurations,
        id.vars = c("job_id", "scheme", "split_id", "fold_id"),
        measure.vars = available_parameters,
        variable.name = "parameter",
        value.name = "archive_value",
        variable.factor = FALSE
      )
      best_parameters <- merge(
        best_parameters,
        parameter_space,
        by = "parameter",
        all.x = TRUE
      )
      best_parameters[, value := archive_value]
      best_parameters[
        log_scale == TRUE & integer_log_scale == FALSE,
        value := exp(archive_value)
      ]
      best_parameters[
        integer_log_scale == TRUE,
        value := pmin(
          upper_natural,
          pmax(lower_natural, floor(exp(archive_value)))
        )
      ]
      best_parameters[, normalized_position :=
                        (archive_value - lower_archive) /
                        (upper_archive - lower_archive)]
      best_parameters[, boundary := fcase(
        normalized_position <= MBO_BOUNDARY_FRACTION, "lower",
        normalized_position >= 1 - MBO_BOUNDARY_FRACTION, "upper",
        default = "interior"
      )]
      
      parameter_summary <- best_parameters[, .(
        n_outer_folds = .N,
        minimum = min(value),
        q25 = as.numeric(quantile(value, 0.25)),
        median = median(value),
        q75 = as.numeric(quantile(value, 0.75)),
        maximum = max(value),
        proportion_near_lower_boundary = mean(boundary == "lower"),
        proportion_near_upper_boundary = mean(boundary == "upper"),
        proportion_near_any_boundary = mean(boundary != "interior")
      ), by = parameter]
      parameter_summary[, systematic_boundary_hit :=
                          pmax(
                            proportion_near_lower_boundary,
                            proportion_near_upper_boundary
                          ) >= MBO_SYSTEMATIC_BOUNDARY_PROPORTION]
      
      fwrite(
        best_parameters,
        file.path(FIGURE_DIR, "tuning_best_parameters_by_fold.csv"),
        na = "NA"
      )
      fwrite(
        parameter_summary,
        file.path(FIGURE_DIR, "tuning_best_parameter_summary.csv"),
        na = "NA"
      )
    } else {
      parameter_summary <- data.table()
      warning("No selected-learner parameter columns were found in the scalar archive.")
    }
    
    tuning[, incumbent := cummin(get(objective_column)), by = job_id]
    tuning[, final_best := min(get(objective_column)), by = job_id]
    tuning[, remaining_improvement_percent := 100 * pmax(
      0,
      (incumbent - final_best) / pmax(abs(final_best), .Machine$double.eps)
    )]
    
    convergence <- tuning[, .(
      n_outer_folds = uniqueN(job_id),
      median_remaining_improvement_percent = median(remaining_improvement_percent),
      q25_remaining_improvement_percent = as.numeric(quantile(
        remaining_improvement_percent,
        0.25
      )),
      q75_remaining_improvement_percent = as.numeric(quantile(
        remaining_improvement_percent,
        0.75
      ))
    ), by = .(scheme, evaluation)]
    convergence[, scheme_display := factor_scheme(scheme)]
    
    convergence_ylim <- as.numeric(quantile(
      convergence$q75_remaining_improvement_percent,
      0.98,
      na.rm = TRUE
    ))
    convergence_ylim <- max(convergence_ylim, 0.1)
    
    p_tuning <- ggplot(
      convergence,
      aes(x = evaluation, y = median_remaining_improvement_percent)
    ) +
      geom_ribbon(
        aes(
          ymin = q25_remaining_improvement_percent,
          ymax = q75_remaining_improvement_percent
        ),
        fill = "grey75",
        alpha = 0.55
      ) +
      geom_line(colour = "#0072B2", linewidth = 0.65) +
      coord_cartesian(ylim = c(0, convergence_ylim)) +
      facet_wrap(vars(scheme_display), scales = "free_x") +
      labs(
        x = "MBO evaluation",
        y = "Remaining improvement\nto final best (%)"
      ) +
      theme_paper(base_size = 9.5) +
      theme(legend.position = "none")
    
    save_publication_plot(
      p_tuning,
      "figure_tuning_convergence_diagnostic",
      width = 7.4,
      height = 5.3
    )
    
    overall_tuning <- tuning_budget_summary[scheme == "all_schemes"]
    systematic_parameters <- if (nrow(parameter_summary)) {
      parameter_summary[systematic_boundary_hit == TRUE, parameter]
    } else {
      character()
    }
    
    budget_supported <-
      overall_tuning$proportion_plateaued >= 0.80 &&
      length(systematic_parameters) == 0L
    
    interpretation <- c(
      "MBO tuning-budget audit",
      "=======================",
      paste0("Selected learner: ", SELECTED_LEARNER),
      paste0("Outer-fold searches audited: ", overall_tuning$n_outer_folds),
      paste0(
        "Proportion plateaued in final quarter: ",
        sprintf("%.3f", overall_tuning$proportion_plateaued)
      ),
      paste0(
        "Median relative improvement in final quarter: ",
        sprintf("%.4f", overall_tuning$median_relative_gain_final_quarter)
      ),
      paste0(
        "Parameters with a systematic boundary hit: ",
        if (length(systematic_parameters)) {
          paste(systematic_parameters, collapse = ", ")
        } else {
          "none"
        }
      ),
      "",
      if (budget_supported) {
        paste(
          "Heuristic conclusion: the archives show no clear evidence that the",
          "MBO budget constrained performance. Additional outer-CV retuning is",
          "not indicated before SHAP."
        )
      } else {
        paste(
          "Heuristic conclusion: inspect the flagged folds or parameter",
          "boundaries before declaring the MBO budget sufficient. A targeted",
          "budget extension on representative folds is preferable to silently",
          "changing the completed benchmark."
        )
      },
      "",
      paste(
        "These rules are diagnostics only: plateau means <=",
        100 * MBO_PLATEAU_RELATIVE_THRESHOLD,
        "% improvement in the final quarter; systematic boundary means >=",
        100 * MBO_SYSTEMATIC_BOUNDARY_PROPORTION,
        "% of fold optima within",
        100 * MBO_BOUNDARY_FRACTION,
        "% of the same search-space boundary."
      )
    )
    writeLines(
      interpretation,
      file.path(FIGURE_DIR, "tuning_budget_interpretation.txt")
    )
  }
}

# ------------------------------------------------------------------------------
# 8. OUTPUT MANIFEST
# ------------------------------------------------------------------------------

output_files <- list.files(
  FIGURE_DIR,
  full.names = TRUE,
  recursive = FALSE
)
output_manifest <- data.table(
  file = basename(output_files),
  path = normalizePath(output_files),
  size_bytes = file.info(output_files)$size,
  modified_time = as.character(file.info(output_files)$mtime)
)
setorder(output_manifest, file)
fwrite(
  output_manifest,
  file.path(FIGURE_DIR, "figure_output_manifest.csv"),
  na = "NA"
)

message("Generalisation figures and diagnostics written to: ", FIGURE_DIR)
message(
  "Primary figures: figure_cv_schemes_paper, ",
  "figure_fold_performance_selected_learner, figure_fold_performance_vs_gdd"
)

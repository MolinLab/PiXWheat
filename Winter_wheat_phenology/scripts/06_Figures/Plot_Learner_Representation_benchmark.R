# ==============================================================================
# PUBLICATION FIGURES FOR POINT 4 AND LEARNER SELECTION
# ==============================================================================
#
# This script consumes the completed outputs of:
#   1. benchmark_features_point4_hypotheses.R
#   2. benchmark_learner_selection.R
#
# It does not refit or retune any learner.
#
# Usage
# -----
# Rscript plot_point4_and_learner_selection.R \
#   [point4_result_dir] [learner_result_dir] [figure_dir]
#
# The same paths can be set with:
#   PHENOLOGY_POINT4_RESULT_DIR
#   PHENOLOGY_LEARNER_RESULT_DIR
#   PHENOLOGY_BENCHMARK_FIGURE_DIR
#
# Point 4 terminology
# -------------------
# The beta-response thermal clock is cumDTU_beta. It is deliberately labelled
# as beta thermal time, not as Soltani GDD.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. USER CONFIGURATION
# ------------------------------------------------------------------------------

parse_boolean <- function(value, default = FALSE) {
  if (!length(value) || is.na(value) || !nzchar(trimws(value))) return(default)
  normalized <- tolower(trimws(value))
  if (normalized %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (normalized %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop("Could not parse a logical value from: ", value)
}

command_line_arguments <- commandArgs(trailingOnly = TRUE)

PROJECT_DIR <- Sys.getenv(
  "PHENOLOGY_PROJECT_DIR",
  unset = getwd()
)

POINT4_RUN_ID <- Sys.getenv(
  "PHENOLOGY_POINT4_RUN_ID",
  unset = "point4_hypothesis_benchmark_paper"
)
LEARNER_RUN_ID <- Sys.getenv(
  "PHENOLOGY_LEARNER_BENCHMARK_RUN_ID",
  unset = "learner_selection_cumulative_hybrid_paper"
)

POINT4_RESULT_DIR <- if (length(command_line_arguments) >= 1L) {
  command_line_arguments[1L]
} else {
  Sys.getenv(
    "PHENOLOGY_POINT4_RESULT_DIR",
    unset = file.path(PROJECT_DIR, "benchmark_results", POINT4_RUN_ID)
  )
}

LEARNER_RESULT_DIR <- if (length(command_line_arguments) >= 2L) {
  command_line_arguments[2L]
} else {
  Sys.getenv(
    "PHENOLOGY_LEARNER_RESULT_DIR",
    unset = file.path(PROJECT_DIR, "benchmark_results", LEARNER_RUN_ID)
  )
}

FIGURE_DIR <- if (length(command_line_arguments) >= 3L) {
  command_line_arguments[3L]
} else {
  Sys.getenv(
    "PHENOLOGY_BENCHMARK_FIGURE_DIR",
    unset = file.path(
      PROJECT_DIR,
      "benchmark_results",
      "paper_figures_point4_learner_selection",
      "new_plots"
    )
  )
}

PRIMARY_FEATURE_SET <- Sys.getenv(
  "PHENOLOGY_PRIMARY_FEATURE_SET",
  unset = "cumulative_weather"
)

ALLOW_INCOMPLETE_BENCHMARKS <- parse_boolean(
  Sys.getenv("PHENOLOGY_ALLOW_INCOMPLETE_BENCHMARKS", unset = "FALSE")
)

PNG_DPI <- as.integer(Sys.getenv("PHENOLOGY_FIGURE_DPI", unset = "600"))
if (is.na(PNG_DPI) || PNG_DPI < 72L) {
  stop("PHENOLOGY_FIGURE_DPI must be an integer >= 72.")
}


# ------------------------------------------------------------------------------
# 2. PACKAGES AND OUTPUT DIRECTORY
# ------------------------------------------------------------------------------

required_packages <- c("data.table", "ggplot2", "patchwork", "systemfonts")
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
})

if (!capabilities("cairo")) {
  stop("Cairo graphics support is required for Liberation Sans PDF and EPS export.")
}
if (!any(tolower(trimws(systemfonts::system_fonts()$family)) == "liberation sans")) {
  stop("Liberation Sans is not available to R on this system.")
}

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 3. LABELS, ORDERS AND GRAPHICAL CONSTANTS
# ------------------------------------------------------------------------------

scheme_order <- c("spatial", "temporal_block")
scheme_labels <- c(
  spatial = "Spatial blocks",
  temporal_block = "Temporal blocks"
)
scheme_colours <- c(
  spatial = "#0072B2",
  temporal_block = "#D55E00"
)

learner_order_default <- c("xgboost", "random_forest", "mars", "gam")
learner_labels <- c(
  xgboost = "XGBoost",
  random_forest = "Random forest",
  mars = "MARS",
  gam = "GAM",
  gdd_threshold = "Calibrated thermal time (DTU)"
)
learner_colours <- c(
  xgboost = "#0072B2",
  random_forest = "#009E73",
  mars = "#CC79A7",
  gam = "#D55E00",
  gdd_threshold = "#777777"
)

feature_order <- c(
  "elapsed_time",
  "gdd_linear",
  "thermal_beta",
  "process_phase",
  "process_phase_water",
  "cumulative_weather",
  "calendar_recent_weather",
  "process_components",
  "process_context",
  "hybrid_process_recent"
)
feature_labels <- c(
  elapsed_time = "Elapsed calendar time",
  gdd_linear = "Linear GDD",
  thermal_beta = "Beta-response thermal time (cumDTU)",
  process_phase = "Phase-aware process clock",
  process_phase_water = "Water-modified process clock",
  cumulative_weather = "Cumulative weather",
  calendar_recent_weather = "Calendar + recent weather",
  process_components = "Process components",
  process_context = "Integrated process index",
  hybrid_process_recent = "Process + recent weather"
)

hypothesis_order <- paste0("H", 1:8)
hypothesis_labels <- c(
  H1 = "H1: Linear GDD vs elapsed time",
  H2 = "H2: Beta-response thermal time vs linear GDD",
  H3 = "H3: Process clock vs beta-response thermal time",
  H4 = "H4: Water modifier vs process clock",
  H5 = "H5: Integrated process vs components",
  H6 = "H6: Hybrid vs integrated process",
  H7 = "H7: Recent weather vs cumulative weather",
  H8 = "H8: Hybrid vs cumulative weather"
)

metric_labels <- c(mae = "MAE", rmse = "RMSE")

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

lookup_labels <- function(x, lookup) {
  original <- as.character(x)
  answer <- unname(lookup[original])
  missing <- is.na(answer)
  answer[missing] <- gsub("_", " ", original[missing], fixed = TRUE)
  answer
}

ordered_display_factor <- function(
    x,
    canonical_order,
    lookup,
    reverse = FALSE
) {
  observed <- unique(as.character(x))
  ordered_values <- c(
    intersect(canonical_order, observed),
    sort(setdiff(observed, canonical_order))
  )
  display_levels <- lookup_labels(ordered_values, lookup)
  if (reverse) display_levels <- rev(display_levels)
  factor(lookup_labels(x, lookup), levels = display_levels)
}

scheme_factor <- function(x) {
  ordered_display_factor(x, scheme_order, scheme_labels)
}

feature_factor <- function(x, reverse = FALSE) {
  ordered_display_factor(x, feature_order, feature_labels, reverse)
}

hypothesis_factor <- function(x, reverse = FALSE) {
  ordered_display_factor(x, hypothesis_order, hypothesis_labels, reverse)
}


# ------------------------------------------------------------------------------
# 4. INPUT AND OUTPUT HELPERS
# ------------------------------------------------------------------------------

assert_directory <- function(path, label) {
  if (!dir.exists(path)) {
    benchmark_root <- file.path(PROJECT_DIR, "benchmark_results")
    candidates <- if (dir.exists(benchmark_root)) {
      list.dirs(benchmark_root, recursive = FALSE, full.names = TRUE)
    } else {
      character()
    }
    candidate_text <- if (length(candidates)) {
      paste0(
        "\nAvailable benchmark directories:\n  ",
        paste(candidates, collapse = "\n  ")
      )
    } else {
      ""
    }
    stop(label, " does not exist: ", path, candidate_text)
  }
}

read_result_table <- function(
    directory,
    filename,
    required_columns,
    label,
    optional = FALSE
) {
  path <- file.path(directory, filename)
  if (!file.exists(path)) {
    if (optional) return(NULL)
    stop("Missing ", label, " file: ", path)
  }
  
  data <- fread(path)
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns)) {
    stop(
      label,
      " is missing columns: ",
      paste(missing_columns, collapse = ", "),
      ". File: ",
      path
    )
  }
  data
}

validate_job_status <- function(directory, benchmark_label) {
  status_path <- file.path(directory, "job_status.csv")
  if (!file.exists(status_path)) {
    warning(
      benchmark_label,
      " has no job_status.csv; completion could not be verified."
    )
    return(invisible(NULL))
  }
  
  status <- fread(status_path)
  if (!"status" %in% names(status)) {
    stop("job_status.csv has no status column: ", status_path)
  }
  
  failed <- status[is.na(status) | status != "ok"]
  if (nrow(failed)) {
    message_text <- paste0(
      benchmark_label,
      " contains ",
      nrow(failed),
      " unfinished or failed outer jobs."
    )
    if (!ALLOW_INCOMPLETE_BENCHMARKS) {
      stop(
        message_text,
        " Set PHENOLOGY_ALLOW_INCOMPLETE_BENCHMARKS=TRUE only for diagnostics."
      )
    }
    warning(message_text)
  }
  
  invisible(status)
}

figure_manifest <- data.table(
  stem = character(),
  analysis = character(),
  description = character(),
  width_in = numeric(),
  height_in = numeric()
)

save_publication_plot <- function(
    plot,
    stem,
    width,
    height,
    analysis,
    description
) {
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
  
  figure_manifest <<- rbind(
    figure_manifest,
    data.table(
      stem = stem,
      analysis = analysis,
      description = description,
      width_in = width,
      height_in = height
    )
  )
  
  invisible(c(pdf = pdf_path, png = png_path, eps = eps_path))
}

write_plot_data <- function(data, filename) {
  fwrite(data, file.path(FIGURE_DIR, filename), na = "NA")
}

assert_directory(POINT4_RESULT_DIR, "Point 4 result directory")
assert_directory(LEARNER_RESULT_DIR, "Learner-selection result directory")
validate_job_status(POINT4_RESULT_DIR, "Point 4 benchmark")
validate_job_status(LEARNER_RESULT_DIR, "Learner-selection benchmark")


# ------------------------------------------------------------------------------
# 5. READ POINT 4 RESULTS
# ------------------------------------------------------------------------------

point4_hypotheses <- read_result_table(
  POINT4_RESULT_DIR,
  "hypothesis_summary.csv",
  c(
    "hypothesis_id", "analysis_family", "learner_id", "scheme",
    "reference", "alternative", "delta_macro_mae",
    "delta_macro_mae_ci_low", "delta_macro_mae_ci_high"
  ),
  "Point 4 hypothesis summary"
)

point4_stage_hypotheses <- read_result_table(
  POINT4_RESULT_DIR,
  "hypothesis_stage_summary.csv",
  c(
    "hypothesis_id", "analysis_family", "learner_id", "scheme",
    "BBCH_stage", "delta_macro_mae"
  ),
  "Point 4 stage-specific hypothesis summary"
)

point4_predictor_summary <- read_result_table(
  POINT4_RESULT_DIR,
  "predictor_set_summary.csv",
  c(
    "analysis_family", "feature_set", "learner_id", "scheme",
    "macro_mae", "rank_macro_mae"
  ),
  "Point 4 predictor-set summary"
)

point4_hypotheses[, scheme_display := scheme_factor(scheme)]
point4_hypotheses[
  ,
  hypothesis_display := hypothesis_factor(hypothesis_id, reverse = TRUE)
]

point4_stage_hypotheses[, BBCH_stage := as.numeric(BBCH_stage)]
point4_stage_hypotheses[, scheme_display := scheme_factor(scheme)]
point4_stage_hypotheses[
  ,
  hypothesis_display := hypothesis_factor(hypothesis_id)
]

write_plot_data(
  point4_hypotheses,
  "plot_data_point4_hypothesis_contrasts.csv"
)


# ------------------------------------------------------------------------------
# 6. POINT 4 FIGURE: MECHANISM HYPOTHESES
# ------------------------------------------------------------------------------

mechanism_hypotheses <- point4_hypotheses[analysis_family == "mechanism"]
if (!nrow(mechanism_hypotheses)) {
  stop("No mechanism rows were found in hypothesis_summary.csv.")
}

mechanism_ci <- mechanism_hypotheses[
  is.finite(delta_macro_mae_ci_low) &
    is.finite(delta_macro_mae_ci_high)
]
mechanism_dodge <- position_dodge(width = 0.55)

p_point4_mechanism <- ggplot(
  mechanism_hypotheses,
  aes(
    x = hypothesis_display,
    y = delta_macro_mae,
    colour = scheme_display,
    group = scheme_display
  )
) +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.45) +
  geom_errorbar(
    data = mechanism_ci,
    aes(
      ymin = delta_macro_mae_ci_low,
      ymax = delta_macro_mae_ci_high
    ),
    position = mechanism_dodge,
    width = 0.16,
    linewidth = 0.55
  ) +
  geom_point(position = mechanism_dodge, size = 2.3) +
  coord_flip() +
  scale_colour_manual(
    values = unname(scheme_colours[scheme_order]),
    breaks = unname(scheme_labels[scheme_order]),
    name = "Outer CV"
  ) +
  labs(
    
    x = NULL,
    y = expression(Delta * " macro MAE (BBCH stages)")
  ) +
  theme_paper(base_size = 10) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 8.5)
  )

save_publication_plot(
  p_point4_mechanism,
  "figure_point4_mechanism_hypotheses",
  width = 7.4,
  height = 4.4,
  analysis = "Point 4",
  description = "Paired mechanism-ablation macro-MAE contrasts with 95% CIs."
)


# ------------------------------------------------------------------------------
# 7. POINT 4 FIGURE: REPRESENTATION HYPOTHESES
# ------------------------------------------------------------------------------

representation_hypotheses <- point4_hypotheses[
  analysis_family == "representation"
]
if (!nrow(representation_hypotheses)) {
  stop("No representation rows were found in hypothesis_summary.csv.")
}

point4_learners <- c(
  intersect(learner_order_default, unique(representation_hypotheses$learner_id)),
  sort(setdiff(
    unique(representation_hypotheses$learner_id),
    learner_order_default
  ))
)
representation_hypotheses[
  ,
  learner_display := ordered_display_factor(
    learner_id,
    point4_learners,
    learner_labels,
    reverse = TRUE
  )
]
representation_hypotheses[
  ,
  hypothesis_display := hypothesis_factor(hypothesis_id)
]

representation_ci <- representation_hypotheses[
  is.finite(delta_macro_mae_ci_low) &
    is.finite(delta_macro_mae_ci_high)
]
representation_dodge <- position_dodge(width = 0.58)

p_point4_representation <- ggplot(
  representation_hypotheses,
  aes(
    x = learner_display,
    y = delta_macro_mae,
    colour = scheme_display,
    group = scheme_display
  )
) +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.42) +
  geom_errorbar(
    data = representation_ci,
    aes(
      ymin = delta_macro_mae_ci_low,
      ymax = delta_macro_mae_ci_high
    ),
    position = representation_dodge,
    width = 0.15,
    linewidth = 0.48
  ) +
  geom_point(position = representation_dodge, size = 1.9) +
  coord_flip() +
  facet_wrap(
    vars(hypothesis_display),
    ncol = 2,
    scales = "free_y"
  ) +
  scale_colour_manual(
    values = unname(scheme_colours[scheme_order]),
    breaks = unname(scheme_labels[scheme_order]),
    name = "Outer CV"
  ) +
  labs(
    
    x = NULL,
    y = expression(Delta * " macro MAE (BBCH stages)")
  ) +
  theme_paper(base_size = 9.5) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 8.2)
  )

save_publication_plot(
  p_point4_representation,
  "figure_point4_representation_hypotheses",
  width = 7.4,
  height = 6.7,
  analysis = "Point 4",
  description = "Paired representation contrasts by learner and outer CV."
)


# ------------------------------------------------------------------------------
# 8. POINT 4 FIGURE: REPRESENTATION SCREENING HEATMAP
# ------------------------------------------------------------------------------

representation_summary <- point4_predictor_summary[
  analysis_family == "representation"
]
if (!nrow(representation_summary)) {
  stop("No representation rows were found in predictor_set_summary.csv.")
}

representation_summary[, scheme_display := scheme_factor(scheme)]
representation_summary[
  ,
  feature_display := feature_factor(feature_set, reverse = TRUE)
]
representation_summary[
  ,
  learner_display := ordered_display_factor(
    learner_id,
    point4_learners,
    learner_labels
  )
]
representation_summary[
  ,
  column_display := paste(learner_display, scheme_display, sep = "\n")
]

column_order <- unlist(lapply(point4_learners, function(learner_id_value) {
  learner_label_value <- lookup_labels(learner_id_value, learner_labels)
  paste(
    learner_label_value,
    unname(scheme_labels[scheme_order]),
    sep = "\n"
  )
}))
representation_summary[
  ,
  column_display := factor(column_display, levels = column_order)
]
representation_summary[
  ,
  cell_label := sprintf("%.2f\n(rank %.1f)", macro_mae, rank_macro_mae)
]

p_point4_ranks <- ggplot(
  representation_summary,
  aes(x = column_display, y = feature_display, fill = rank_macro_mae)
) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = cell_label), family = "Liberation Sans", size = 2.50, lineheight = 0.90) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#9ECAE1",
    name = "MAE rank"
  ) +
  labs(
    
    x = NULL,
    y = NULL
  ) +
  theme_paper(base_size = 9.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 7.0, lineheight = 0.90, angle = 25, hjust= 1, vjust= 1),
    axis.text.y = element_text(size = 8.2),
    legend.position = "right"
  )

save_publication_plot(
  p_point4_ranks,
  "figure_point4_representation_ranks",
  width = 7.4,
  height = 4.7,
  analysis = "Point 4",
  description = "Macro-MAE values and ranks for representation screening."
)

write_plot_data(
  representation_summary,
  "plot_data_point4_representation_ranks.csv"
)

library(patchwork)
combined_plot <- (
  p_point4_mechanism /
    p_point4_representation /
    p_point4_ranks
) + plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")
combined_plot

save_publication_plot(
  combined_plot,
  "figure_combinend",
  width = 9,
  height = 11,
  analysis = "Point 4",
  description = ""
)

# ------------------------------------------------------------------------------
# 9. POINT 4 SUPPLEMENT: STAGE-SPECIFIC HYPOTHESIS DELTAS
# ------------------------------------------------------------------------------

stage_breaks <- sort(unique(point4_stage_hypotheses$BBCH_stage))

mechanism_stage <- point4_stage_hypotheses[analysis_family == "mechanism"]
if (nrow(mechanism_stage)) {
  p_point4_mechanism_stage <- ggplot(
    mechanism_stage,
    aes(
      x = BBCH_stage,
      y = delta_macro_mae,
      colour = scheme_display,
      group = scheme_display
    )
  ) +
    geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.35) +
    geom_line(linewidth = 0.60) +
    geom_point(size = 1.45) +
    facet_wrap(
      vars(hypothesis_display),
      ncol = 2,
      scales = "free_y"
    ) +
    scale_x_continuous(breaks = stage_breaks) +
    scale_colour_manual(
      values = unname(scheme_colours[scheme_order]),
      breaks = unname(scheme_labels[scheme_order]),
      name = "Outer CV"
    ) +
    labs(
      
      x = "Observed BBCH stage",
      y = expression(Delta * " macro MAE (BBCH stages)")
    ) +
    theme_paper(base_size = 9.3) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 8.1)
    )
  
  save_publication_plot(
    p_point4_mechanism_stage,
    "figure_point4_stage_deltas_mechanism",
    width = 7.4,
    height = 6.1,
    analysis = "Point 4 supplement",
    description = "Stage-specific mechanism-ablation MAE contrasts."
  )
}

representation_stage <- point4_stage_hypotheses[
  analysis_family == "representation"
]
if (nrow(representation_stage)) {
  representation_stage[
    ,
    learner_display := ordered_display_factor(
      learner_id,
      point4_learners,
      learner_labels
    )
  ]
  point4_learner_palette <- learner_colours[point4_learners]
  names(point4_learner_palette) <- lookup_labels(
    point4_learners,
    learner_labels
  )
  
  p_point4_representation_stage <- ggplot(
    representation_stage,
    aes(
      x = BBCH_stage,
      y = delta_macro_mae,
      colour = learner_display,
      group = learner_display
    )
  ) +
    geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.35) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.15) +
    facet_grid(
      rows = vars(hypothesis_display),
      cols = vars(scheme_display),
      scales = "free_y"
    ) +
    scale_x_continuous(breaks = stage_breaks) +
    scale_colour_manual(
      values = point4_learner_palette,
      name = "Learner"
    ) +
    labs(
      
      x = "Observed BBCH stage",
      y = expression(Delta * " macro MAE (BBCH stages)")
    ) +
    theme_paper(base_size = 8.8) +
    theme(
      legend.position = "bottom",
      strip.text.y = element_text(size = 7.5),
      strip.text.x = element_text(size = 8.0),
      axis.text.x = element_text(size = 7.0)
    )
  
  save_publication_plot(
    p_point4_representation_stage,
    "figure_point4_stage_deltas_representation",
    width = 7.4,
    height = 8.4,
    analysis = "Point 4 supplement",
    description = "Stage-specific representation MAE contrasts."
  )
}


# ------------------------------------------------------------------------------
# 10. READ LEARNER-SELECTION RESULTS
# ------------------------------------------------------------------------------

learner_environment_stage <- read_result_table(
  LEARNER_RESULT_DIR,
  "environment_stage_metrics.csv",
  c("feature_set", "learner_id", "scheme", "mae", "rmse"),
  "Learner-selection environment-stage metrics"
)

learner_stage_metrics <- read_result_table(
  LEARNER_RESULT_DIR,
  "stage_specific_metrics.csv",
  c(
    "feature_set", "learner_id", "scheme", "BBCH_stage",
    "mae", "rmse"
  ),
  "Learner-selection stage metrics"
)

learner_selection <- read_result_table(
  LEARNER_RESULT_DIR,
  "learner_selection_primary_cumulative.csv",
  c("learner_id", "selection_order",  "mean_macro_mae"),
  "Primary learner-selection summary"
)

paired_learners <- read_result_table(
  LEARNER_RESULT_DIR,
  "paired_learner_bootstrap_primary_cumulative.csv",
  c(
    "scheme", "learner_a", "learner_b", "observed_delta",
    "ci_low", "ci_high"
  ),
  "Paired learner bootstrap",
  optional = TRUE
)

representation_sensitivity <- read_result_table(
  LEARNER_RESULT_DIR,
  "hybrid_sensitivity_vs_cumulative.csv",
  c(
    "learner_id", "scheme",
    "delta_mae_hybrid_minus_cumulative",
    "delta_rmse_hybrid_minus_cumulative"
  ),
  "Learner representation sensitivity",
  optional = TRUE
)

learner_ids_present <- unique(c(
  learner_environment_stage$learner_id,
  learner_stage_metrics$learner_id
))
selection_order_ids <- learner_selection[
  order(selection_order),
  as.character(learner_id)
]
learner_order <- c(
  intersect(selection_order_ids, learner_ids_present),
  intersect(learner_order_default, learner_ids_present),
  sort(setdiff(
    learner_ids_present,
    c(selection_order_ids, learner_order_default)
  ))
)
learner_order <- unique(learner_order)

unknown_learner_ids <- setdiff(learner_order, names(learner_colours))
if (length(unknown_learner_ids)) {
  additional_colours <- grDevices::hcl.colors(
    length(unknown_learner_ids),
    palette = "Dark 3"
  )
  names(additional_colours) <- unknown_learner_ids
  learner_colours <- c(learner_colours, additional_colours)
}

learner_palette <- learner_colours[learner_order]
names(learner_palette) <- lookup_labels(learner_order, learner_labels)


# ------------------------------------------------------------------------------
# 11. LEARNER FIGURE: ENVIRONMENT-STAGE ERROR DISTRIBUTIONS
# ------------------------------------------------------------------------------

primary_environment_stage <- learner_environment_stage[
  feature_set == PRIMARY_FEATURE_SET
]
if (!nrow(primary_environment_stage)) {
  stop(
    "No environment-stage rows were found for PRIMARY_FEATURE_SET = ",
    PRIMARY_FEATURE_SET
  )
}

learner_error_long <- melt(
  primary_environment_stage,
  id.vars = c("learner_id", "scheme"),
  measure.vars = c("mae", "rmse"),
  variable.name = "metric",
  value.name = "error"
)
learner_error_long <- learner_error_long[is.finite(error)]

learner_box_summary <- learner_error_long[, {
  quartiles <- as.numeric(quantile(error, c(0.25, 0.50, 0.75)))
  interquartile_range <- quartiles[3L] - quartiles[1L]
  lower_fence <- quartiles[1L] - 1.5 * interquartile_range
  upper_fence <- quartiles[3L] + 1.5 * interquartile_range
  
  list(
    n_environment_stage_units = .N,
    whisker_low = min(error[error >= lower_fence]),
    q25 = quartiles[1L],
    median = quartiles[2L],
    q75 = quartiles[3L],
    whisker_high = max(error[error <= upper_fence]),
    mean_error = mean(error)
  )
}, by = .(learner_id, scheme, metric)]

learner_box_summary[, scheme_display := scheme_factor(scheme)]
learner_box_summary[
  ,
  learner_display := ordered_display_factor(
    learner_id,
    learner_order,
    learner_labels,
    reverse = TRUE
  )
]
learner_box_summary[
  ,
  metric_display := factor(
    lookup_labels(metric, metric_labels),
    levels = unname(metric_labels[c("mae", "rmse")])
  )
]

p_learner_boxplots <- ggplot(
  learner_box_summary,
  aes(x = learner_display, fill = learner_display)
) +
  geom_boxplot(
    aes(
      ymin = whisker_low,
      lower = q25,
      middle = median,
      upper = q75,
      ymax = whisker_high
    ),
    stat = "identity",
    width = 0.64,
    linewidth = 0.46,
    alpha = 0.80
  ) +
  #geom_point(
  #  aes(y = mean_error),
  #  shape = 23,
  #  size = 2.2,
  #  stroke = 0.55,
  #  fill = "white",
  # colour = "black"
  #) +
  coord_flip() +
  facet_grid(
    rows = vars(metric_display),
    cols = vars(scheme_display),
    scales = "free_y"
  ) +
  scale_fill_manual(values = learner_palette, guide = "none") +
  labs(
    
    x = NULL,
    y = "Prediction error (BBCH stages)"
  ) +
  theme_paper(base_size = 9.7) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 8.2),
    strip.text = element_text(size = 8.5)
  )

save_publication_plot(
  p_learner_boxplots,
  "figure_learner_selection_boxplots",
  width = 7.4,
  height = 6.0,
  analysis = "Learner selection",
  description = paste(
    "Median, IQR, whiskers and macro means for environment-stage",
    "MAE and RMSE."
  )
)

write_plot_data(
  learner_box_summary,
  "learner_environment_stage_error_quantiles.csv"
)


# ------------------------------------------------------------------------------
# 12. LEARNER FIGURES: STAGE-SPECIFIC MAE AND RMSE
# ------------------------------------------------------------------------------

primary_stage_metrics <- learner_stage_metrics[
  feature_set == PRIMARY_FEATURE_SET
]
primary_stage_metrics[, BBCH_stage := as.numeric(BBCH_stage)]
primary_stage_metrics[, scheme_display := scheme_factor(scheme)]
primary_stage_metrics[
  ,
  learner_display := ordered_display_factor(
    learner_id,
    learner_order,
    learner_labels
  )
]

learner_stage_breaks <- sort(unique(primary_stage_metrics$BBCH_stage))

make_stage_performance_plot <- function(metric_column, metric_label) {
  plot_data <- copy(primary_stage_metrics)
  plot_data[, error := get(metric_column)]
  
  ggplot(
    plot_data,
    aes(
      x = BBCH_stage,
      y = error,
      colour = learner_display,
      group = learner_display
    )
  ) +
    geom_line(linewidth = 0.68, alpha = 0.90) +
    geom_point(size = 1.60) +
    facet_wrap(vars(scheme_display), nrow = 1, scales = "free_y") +
    scale_x_continuous(breaks = learner_stage_breaks) +
    scale_colour_manual(values = learner_palette, name = "Learner") +
    labs(
      
      x = "Observed BBCH stage",
      y = paste(metric_label, "(BBCH stages)")
    ) +
    theme_paper(base_size = 9.7) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.text.x = element_text(size = 8)
    )
}

p_learner_stage_mae <- make_stage_performance_plot("mae", "MAE")
save_publication_plot(
  p_learner_stage_mae,
  "figure_learner_stage_mae",
  width = 7.4,
  height = 4.2,
  analysis = "Learner selection",
  description = "Stage-specific out-of-fold MAE by learner and scheme."
)

p_learner_stage_rmse <- make_stage_performance_plot("rmse", "RMSE")
save_publication_plot(
  p_learner_stage_rmse,
  "figure_learner_stage_rmse",
  width = 7.4,
  height = 4.2,
  analysis = "Learner selection supplement",
  description = "Stage-specific out-of-fold RMSE by learner and scheme."
)


# ------------------------------------------------------------------------------
# 13. LEARNER FIGURE: PAIRED BOOTSTRAP CONTRASTS
# ------------------------------------------------------------------------------

if (!is.null(paired_learners) && nrow(paired_learners)) {
  paired_learners[, scheme_display := scheme_factor(scheme)]
  paired_learners[
    ,
    pair_display := paste(
      lookup_labels(learner_a, learner_labels),
      lookup_labels(learner_b, learner_labels),
      sep = " - "
    )
  ]
  pair_order <- unique(paired_learners$pair_display)
  paired_learners[
    ,
    pair_display := factor(pair_display, levels = rev(pair_order))
  ]
  
  paired_ci <- paired_learners[is.finite(ci_low) & is.finite(ci_high)]
  paired_dodge <- position_dodge(width = 0.58)
  
  p_paired_learners <- ggplot(
    paired_learners,
    aes(
      x = pair_display,
      y = observed_delta,
      colour = scheme_display,
      group = scheme_display
    )
  ) +
    geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.45) +
    geom_errorbar(
      data = paired_ci,
      aes(ymin = ci_low, ymax = ci_high),
      position = paired_dodge,
      width = 0.16,
      linewidth = 0.52
    ) +
    geom_point(position = paired_dodge, size = 2.0) +
    coord_flip() +
    scale_colour_manual(
      values = unname(scheme_colours[scheme_order]),
      breaks = unname(scheme_labels[scheme_order]),
      name = "Outer CV"
    ) +
    labs(
      
      x = "Learner A - learner B",
      y = expression(Delta * " macro MAE: A - B (BBCH stages)")
    ) +
    theme_paper(base_size = 9.7) +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 8.1)
    )
  
  save_publication_plot(
    p_paired_learners,
    "figure_learner_paired_contrasts",
    width = 7.4,
    height = 5.2,
    analysis = "Learner selection",
    description = "Paired macro-MAE differences between learners with 95% CIs."
  )
  
  write_plot_data(
    paired_learners,
    "plot_data_learner_paired_contrasts.csv"
  )
} else {
  message(
    "Skipping paired learner figure: ",
    "paired_learner_bootstrap_primary_cumulative.csv is unavailable."
  )
}


# ------------------------------------------------------------------------------
# 14. LEARNER SUPPLEMENT: REPRESENTATION SENSITIVITY
# ------------------------------------------------------------------------------

if (!is.null(representation_sensitivity) && nrow(representation_sensitivity)) {
  sensitivity_long <- melt(
    representation_sensitivity,
    id.vars = c("learner_id", "scheme"),
    measure.vars = c(
      "delta_mae_hybrid_minus_cumulative",
      "delta_rmse_hybrid_minus_cumulative"
    ),
    variable.name = "metric",
    value.name = "delta"
  )
  sensitivity_long[
    ,
    metric_display := factor(
      fifelse(
        metric == "delta_mae_hybrid_minus_cumulative",
        "MAE",
        "RMSE"
      ),
      levels = c("MAE", "RMSE")
    )
  ]
  sensitivity_long[, scheme_display := scheme_factor(scheme)]
  sensitivity_long[
    ,
    learner_display := ordered_display_factor(
      learner_id,
      learner_order,
      learner_labels,
      reverse = TRUE
    )
  ]
  
  sensitivity_dodge <- position_dodge(width = 0.55)
  p_representation_sensitivity <- ggplot(
    sensitivity_long,
    aes(
      x = learner_display,
      y = delta,
      colour = scheme_display,
      group = scheme_display
    )
  ) +
    geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.45) +
    geom_point(position = sensitivity_dodge, size = 2.2) +
    coord_flip() +
    facet_wrap(vars(metric_display), nrow = 1, scales = "free_y") +
    scale_colour_manual(
      values = unname(scheme_colours[scheme_order]),
      breaks = unname(scheme_labels[scheme_order]),
      name = "Outer CV"
    ) +
    labs(
      
      x = NULL,
      y = expression(Delta * " error (BBCH stages)")
    ) +
    theme_paper(base_size = 9.7) +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank()
    )
  
  save_publication_plot(
    p_representation_sensitivity,
    "figure_learner_representation_sensitivity",
    width = 7.4,
    height = 4.4,
    analysis = "Learner selection supplement",
    description = "Hybrid-minus-cumulative sensitivity by learner and scheme."
  )
}


# ------------------------------------------------------------------------------
# 15. AUDIT TABLES AND COMPLETION MESSAGE
# ------------------------------------------------------------------------------

learner_selection_export <- copy(learner_selection)
learner_selection_export[
  ,
  learner_label := lookup_labels(learner_id, learner_labels)
]
write_plot_data(
  learner_selection_export,
  "learner_selection_order_for_figures.csv"
)

write_plot_data(figure_manifest, "benchmark_figure_manifest.csv")

configuration <- data.table(
  setting = c(
    "point4_result_dir",
    "learner_result_dir",
    "figure_dir",
    "primary_feature_set",
    "beta_thermal_clock_label",
    "allow_incomplete_benchmarks",
    "png_dpi"
  ),
  value = c(
    normalizePath(POINT4_RESULT_DIR),
    normalizePath(LEARNER_RESULT_DIR),
    normalizePath(FIGURE_DIR),
    PRIMARY_FEATURE_SET,
    "Beta-response thermal time (cumDTU)",
    as.character(ALLOW_INCOMPLETE_BENCHMARKS),
    as.character(PNG_DPI)
  )
)
write_plot_data(configuration, "benchmark_plot_configuration.csv")

message("Benchmark plotting finished successfully.")
message("Figures written to: ", normalizePath(FIGURE_DIR))
message("Number of figure stems: ", nrow(figure_manifest))

# ==============================================================================
# PUBLICATION FIGURE: OBSERVATION SUPPORT ACROSS BBCH STAGES
# ==============================================================================
#
# This script describes the frozen development cohort used by
# multischeme_generalisation_benchmark.R. It does not use model predictions.
#
# Main outputs
# ------------
#   * figure_observation_support_by_bbch.pdf
#   * figure_observation_support_by_bbch.png
#   * figure_observation_support_by_bbch.eps
#   * observation_support_by_exact_bbch.csv
#   * observation_support_by_principal_stage.csv
#   * figure_observation_support_by_bbch_caption.txt
#
# Panel A shows the number of observations at each exact BBCH code. A logarithmic
# axis is used because the scheduled stages contain tens of thousands of records,
# whereas several exact stages contain fewer than ten observations.
#
# Panel B separates raw sample size from independent support by showing the
# numbers of environments, sites, and season years represented at each stage.
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

MANIFEST_FILE <- Sys.getenv(
  "PHENOLOGY_OBSERVATION_MANIFEST",
  unset = file.path(RESULT_DIR, "outer_fold_manifest.csv")
)

FIGURE_DIR <- Sys.getenv(
  "PHENOLOGY_GENERALISATION_FIGURE_DIR",
  unset = file.path(RESULT_DIR, "paper_figures")
)

PNG_DPI <- 600

# Stages labelled directly in panel A. BBCH 50 and 59 are included because their
# sparse support is important for interpreting stage-specific model results.
LABELLED_STAGES <- c(0, 10, 31, 51, 60, 75, 87, 89)

# This flag is descriptive only; it is not used to exclude observations.
SPARSE_STAGE_N <- 30L

# ------------------------------------------------------------------------------
# 2. PACKAGES AND GRAPHICAL CONSTANTS
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

stage_names <- c(
  "0" = "Germination",
  "1" = "Leaf development",
  "2" = "Tillering",
  "3" = "Stem elongation",
  "4" = "Booting",
  "5" = "Heading",
  "6" = "Flowering",
  "7" = "Milk development",
  "8" = "Ripening",
  "9" = "Senescence"
)

support_colours <- c(
  "Environments" = "#0072B2",
  "Sites" = "#009E73",
  "Season years" = "#CC79A7"
)

theme_publication <- function(base_size = 9.5) {
  theme_bw(base_size = base_size, base_family = "Liberation Sans") +
    theme(
      text = element_text(family = "Liberation Sans"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title = element_text(colour = "black"),
      axis.text = element_text(colour = "black"),
      legend.position = "bottom",
      legend.justification = "center",
      legend.title = element_blank(),
      plot.tag = element_text(family = "Liberation Sans", face = "bold", size = 11),
      plot.margin = margin(4, 6, 4, 4)
    )
}

# ------------------------------------------------------------------------------
# 3. LOAD AND VALIDATE THE FROZEN ANALYSIS COHORT
# ------------------------------------------------------------------------------

if (!file.exists(MANIFEST_FILE)) {
  stop(
    "Frozen outer-fold manifest not found: ", MANIFEST_FILE, "\n",
    "Set PHENOLOGY_OBSERVATION_MANIFEST to outer_fold_manifest.csv from the ",
    "completed multischeme benchmark."
  )
}

manifest <- fread(MANIFEST_FILE)

required_columns <- c(
  "obs_id", "env_id", "site_id", "season_year", "BBCH"
)
missing_columns <- setdiff(required_columns, names(manifest))
if (length(missing_columns)) {
  stop(
    "outer_fold_manifest.csv is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (anyDuplicated(manifest$obs_id)) {
  stop("outer_fold_manifest.csv contains duplicated obs_id values.")
}

manifest[, BBCH := suppressWarnings(as.numeric(BBCH))]
manifest[, season_year := suppressWarnings(as.integer(season_year))]

if (manifest[, any(!is.finite(BBCH))]) {
  stop("BBCH contains missing or non-numeric values in the frozen manifest.")
}
if (manifest[, any(BBCH < 0 | BBCH > 99)]) {
  stop("BBCH values outside the expected 0--99 range were found.")
}
if (manifest[, anyNA(env_id) || anyNA(site_id) || anyNA(season_year)]) {
  stop("env_id, site_id, or season_year contains missing values in the manifest.")
}

# ------------------------------------------------------------------------------
# 4. SUMMARIZE EXACT AND PRINCIPAL BBCH-STAGE SUPPORT
# ------------------------------------------------------------------------------

stage_support <- manifest[
  , .(
    n_observations = .N,
    n_environments = uniqueN(env_id),
    n_sites = uniqueN(site_id),
    n_season_years = uniqueN(season_year)
  ),
  by = BBCH
][order(BBCH)]

stage_support[, principal_stage := pmin(9L, floor(BBCH / 10))]
stage_support[, principal_stage_name := unname(stage_names[as.character(principal_stage)])]
stage_support[, sparse_exact_stage := n_observations < SPARSE_STAGE_N]
stage_support[, labelled_stage := BBCH %in% LABELLED_STAGES]
stage_support[, observation_label := fifelse(
  labelled_stage,
  paste0("BBCH ", BBCH, "\n", scales::comma(n_observations)),
  NA_character_
)]
stage_support[, observation_label_y := n_observations * 1.42]
# Separate the adjacent late-stage labels at journal figure width.
stage_support[, observation_label_x := fifelse(
  BBCH == 87, BBCH - 4,
  fifelse(BBCH == 89, BBCH + 1, BBCH)
)]

principal_support <- manifest[
  , principal_stage := pmin(9L, floor(BBCH / 10))
][
  , .(
    n_observations = .N,
    n_exact_stages = uniqueN(BBCH),
    n_environments = uniqueN(env_id),
    n_sites = uniqueN(site_id),
    n_season_years = uniqueN(season_year)
  ),
  by = principal_stage
][order(principal_stage)]
principal_support[, principal_stage_name := unname(
  stage_names[as.character(principal_stage)]
)]
principal_support[, observation_percentage := 100 * n_observations / sum(n_observations)]

fwrite(
  stage_support,
  file.path(FIGURE_DIR, "observation_support_by_exact_bbch.csv")
)
fwrite(
  principal_support,
  file.path(FIGURE_DIR, "observation_support_by_principal_stage.csv")
)

# ------------------------------------------------------------------------------
# 5. BUILD THE TWO-PANEL PUBLICATION FIGURE
# ------------------------------------------------------------------------------

major_breaks <- seq(0, 90, by = 5)
count_breaks <- 10^(0:5)
max_count <- max(stage_support$n_observations)

p_observations <- ggplot(
  stage_support,
  aes(x = BBCH, y = n_observations)
) +
  geom_segment(
    aes(xend = BBCH, y = 1, yend = n_observations),
    linewidth = 1.25,
    colour = "#7A8B99",
    lineend = "butt"
  ) +
  geom_point(
    aes(
      fill = sparse_exact_stage,
      size = labelled_stage
    ),
    shape = 21,
    stroke = 0.45,
    colour = "white"
  ) +
  geom_text(
    data = stage_support[labelled_stage == TRUE],
    aes(x = observation_label_x, y = observation_label_y * 1.5,
        label = observation_label),
    size = 2.55,
    lineheight = .92,
    check_overlap = TRUE,
    colour = "grey15"
  ) +
  scale_fill_manual(
    values = c("FALSE" = "#365F79", "TRUE" = "#D55E00"),
    breaks = c(FALSE, TRUE),
    labels = c(
      paste0("At least ", SPARSE_STAGE_N, " observations"),
      paste0("Fewer than ", SPARSE_STAGE_N, " observations")
    )
  ) +
  scale_size_manual(values = c("FALSE" = 1.9, "TRUE" = 2.8), guide = "none") +
  scale_x_continuous(
    breaks = major_breaks,
    minor_breaks = NULL,
    limits = c(-5, 95),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_log10(
    breaks = count_breaks[count_breaks <= 10^ceiling(log10(max_count))],
    labels = scales::label_comma(),
    limits = c(1, max_count * 5.5),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "Number of observations (log scale)",
    fill = NULL
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3.2))) +
  theme_publication()

support_long <- melt(
  stage_support,
  id.vars = c("BBCH", "sparse_exact_stage"),
  measure.vars = c("n_environments", "n_sites", "n_season_years"),
  variable.name = "support_unit",
  value.name = "n_units"
)
support_long[, support_unit := factor(
  support_unit,
  levels = c("n_environments", "n_sites", "n_season_years"),
  labels = c("Environments", "Sites", "Season years")
)]

support_breaks <- 10^(0:5)
max_support <- max(support_long$n_units)

p_independence <- ggplot(
  support_long,
  aes(x = BBCH, y = n_units, colour = support_unit)
) +
  geom_line(linewidth = 0.55, alpha = 0.75) +
  geom_point(size = 1.65, alpha = 0.95) +
  scale_colour_manual(values = support_colours) +
  scale_x_continuous(
    breaks = major_breaks,
    minor_breaks = NULL,
    limits = c(-5, 95),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_log10(
    breaks = support_breaks[support_breaks <= 10^ceiling(log10(max_support))],
    labels = scales::label_comma(),
    limits = c(1, max_support * 1.65),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Exact BBCH stage",
    y = "Number of independent units (log scale)",
    colour = NULL
  ) +
  theme_publication()

combined_figure <- p_observations / p_independence +
  plot_layout(heights = c(1.08, 1)) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")

figure_pdf <- file.path(FIGURE_DIR, "figure_observation_support_by_bbch.pdf")
figure_png <- file.path(FIGURE_DIR, "figure_observation_support_by_bbch.png")
figure_eps <- file.path(FIGURE_DIR, "figure_observation_support_by_bbch.eps")

ggsave(
  figure_pdf,
  combined_figure,
  width = 7.4,
  height = 8.2,
  units = "in",
  device = cairo_pdf
)
ggsave(
  figure_png,
  combined_figure,
  width = 7.4,
  height = 8.2,
  units = "in",
  dpi = PNG_DPI,
  bg = "white"
)
ggsave(
  figure_eps,
  combined_figure,
  width = 7.4,
  height = 8.2,
  units = "in",
  device = grDevices::cairo_ps,
  onefile = FALSE,
  bg = "white"
)

# ------------------------------------------------------------------------------
# 6. WRITE A READY-TO-EDIT FIGURE CAPTION
# ------------------------------------------------------------------------------

bbch50_n <- stage_support[BBCH == 50, n_observations]
bbch59_n <- stage_support[BBCH == 59, n_observations]
if (!length(bbch50_n)) bbch50_n <- 0L
if (!length(bbch59_n)) bbch59_n <- 0L

caption <- paste0(
  "Figure X. Distribution and independent support of phenological observations ",
  "in the frozen development cohort. (A) Number of observations at each exact ",
  "BBCH code on a logarithmic scale. Orange points identify exact stages with ",
  "fewer than ", SPARSE_STAGE_N, " observations; this threshold is descriptive ",
  "and was not used for filtering. (B) Numbers of environments, sites, and ",
  "season years represented at each exact stage, also on a logarithmic scale. ",
  "The cohort comprised ", scales::comma(nrow(manifest)), " observations across ",
  nrow(stage_support), " exact BBCH codes. Exact BBCH 50 and BBCH 59 were ",
  "represented by ", scales::comma(bbch50_n), " and ", scales::comma(bbch59_n),
  " observations, respectively, and their stage-specific estimates should ",
  "therefore be interpreted cautiously."
)

writeLines(
  caption,
  file.path(FIGURE_DIR, "figure_observation_support_by_bbch_caption.txt")
)

message("Observation-support figure and summaries written to: ", FIGURE_DIR)


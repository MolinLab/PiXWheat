#!/usr/bin/env Rscript

# ==============================================================================
# Environment-level BRIWECS milestone targets
# ==============================================================================
# Estimate BBCH10, BBCH59 and BBCH87 dates for the reference treatment.
# Average cultivar BLUEs with equal weight after adjusting for field row and
# column. If the mixed model cannot be fitted, use equal-weight cultivar means.
# Both routes create one target per location, season and milestone.
# Run from the repository root or set PHENOLOGY_PROJECT_DIR to its path.

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(emmeans)
})

# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------

PROJECT_DIR <- Sys.getenv("PHENOLOGY_PROJECT_DIR", unset = getwd())
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = TRUE)

BRIWECS_DIR <- Sys.getenv(
  "PHENOLOGY_BRIWECS_DIR",
  unset = file.path(PROJECT_DIR, "Data", "BRIWECS")
)

INPUT_FILE <- Sys.getenv(
  "PHENOLOGY_BRIWECS_PHENOLOGY_FILE",
  unset = file.path(BRIWECS_DIR, "BRIWECS_data_publication.csv")
)

OUTPUT_DIR <- Sys.getenv(
  "PHENOLOGY_BRIWECS_TARGET_OUTPUT_DIR",
  unset = file.path(
    PROJECT_DIR,
    "external_validation",
    "briwecs_environment_targets"
  )
)

REFERENCE_TREATMENT <- Sys.getenv(
  "PHENOLOGY_BRIWECS_REFERENCE_TREATMENT",
  unset = "HN_WF_RF"
)

SINGULAR_TOLERANCE <- 1e-5
MAX_DAYS_AFTER_SOWING <- 400L

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2. Helpers
# ------------------------------------------------------------------------------

resolve_column <- function(data, candidates, label, required = TRUE) {
  hits <- candidates[candidates %in% names(data)]
  if (!length(hits)) {
    if (!required) return(NA_character_)
    stop(
      "Could not identify ", label, ". Expected one of: ",
      paste(candidates, collapse = ", ")
    )
  }
  hits[[1L]]
}

numeric_locale <- function(x) {
  suppressWarnings(as.numeric(gsub(
    ",", ".", trimws(as.character(x)), fixed = TRUE
  )))
}

clean_identifier <- function(x) {
  value <- trimws(as.character(x))
  value[!nzchar(value) | value %in% c("NA", "NaN")] <- NA_character_
  value
}

date_from_year_doy <- function(year, doy) {
  result <- as.IDate(rep(NA_character_, length(year)))
  keep <- !is.na(year) & !is.na(doy)
  result[keep] <- as.IDate(sprintf("%d-01-01", year[keep])) + doy[keep] - 1L
  result
}

median_date <- function(x) {
  value <- as.numeric(x)
  value <- value[is.finite(value)]
  if (!length(value)) return(as.IDate(NA))
  as.IDate(round(stats::median(value)), origin = "1970-01-01")
}

safe_sd <- function(x) {
  value <- as.numeric(x)
  value <- value[is.finite(value)]
  if (length(value) < 2L) return(NA_real_)
  as.numeric(stats::sd(value))
}

safe_iqr <- function(x) {
  value <- as.numeric(x)
  value <- value[is.finite(value)]
  if (!length(value)) return(NA_real_)
  as.numeric(stats::IQR(value))
}

fit_paper_target <- function(data) {
  model_data <- copy(data[
    is.finite(days_after_sowing) &
      !is.na(genotype) &
      !is.na(row) &
      !is.na(column)
  ])

  result <- list(
    status = "not_attempted",
    estimate = NA_real_,
    standard_error = NA_real_,
    lower_95 = NA_real_,
    upper_95 = NA_real_,
    singular = NA,
    n_model_rows = nrow(model_data),
    message = NA_character_
  )

  design_is_estimable <-
    nrow(model_data) >= 4L &&
    uniqueN(model_data$genotype) >= 2L &&
    uniqueN(model_data$row) >= 2L &&
    uniqueN(model_data$column) >= 2L &&
    uniqueN(model_data$row) < nrow(model_data) &&
    uniqueN(model_data$column) < nrow(model_data)

  if (!design_is_estimable) {
    result$status <- "insufficient_paper_design"
    result$message <- paste0(
      "Need at least four complete plots, two cultivars, two rows and two ",
      "columns, with replicated row and column levels."
    )
    return(result)
  }

  model_data[, `:=`(
    genotype = factor(genotype),
    row = factor(row),
    column = factor(column)
  )]

  captured_warnings <- character()
  fit <- withCallingHandlers(
    tryCatch(
      lmer(
        days_after_sowing ~ genotype + (1 | row) + (1 | column),
        data = model_data,
        REML = TRUE,
        control = lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 2e5)
        )
      ),
      error = function(error) error
    ),
    warning = function(warning) {
      captured_warnings <<- c(captured_warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )

  if (inherits(fit, "error")) {
    result$status <- "paper_model_failed"
    result$message <- conditionMessage(fit)
    return(result)
  }

  cultivar_blue <- tryCatch(
    emmeans(fit, specs = ~ genotype),
    error = function(error) error
  )
  if (inherits(cultivar_blue, "error")) {
    result$status <- "paper_emmeans_failed"
    result$message <- conditionMessage(cultivar_blue)
    return(result)
  }

  n_cultivars <- length(levels(model_data$genotype))
  equal_average <- tryCatch(
    contrast(
      cultivar_blue,
      method = list(
        equal_cultivar_average = rep(1 / n_cultivars, n_cultivars)
      )
    ),
    error = function(error) error
  )
  if (inherits(equal_average, "error")) {
    result$status <- "paper_contrast_failed"
    result$message <- conditionMessage(equal_average)
    return(result)
  }

  estimate <- as.data.table(summary(
    equal_average,
    infer = c(TRUE, TRUE),
    level = 0.95
  ))
  estimate_column <- intersect(c("estimate", "emmean", "response"), names(estimate))
  lower_column <- grep("^lower", names(estimate), value = TRUE)
  upper_column <- grep("^upper", names(estimate), value = TRUE)
  if (
    nrow(estimate) != 1L ||
    !length(estimate_column) ||
    !length(lower_column) ||
    !length(upper_column)
  ) {
    result$status <- "unexpected_paper_estimate"
    result$message <- paste(names(estimate), collapse = ", ")
    return(result)
  }

  result$singular <- lme4::isSingular(fit, tol = SINGULAR_TOLERANCE)
  result$status <- if (result$singular) "paper_model_singular" else "paper_model_ok"
  result$estimate <- as.numeric(estimate[[estimate_column[[1L]]]][[1L]])
  result$standard_error <- as.numeric(estimate$SE[[1L]])
  result$lower_95 <- as.numeric(estimate[[lower_column[[1L]]]][[1L]])
  result$upper_95 <- as.numeric(estimate[[upper_column[[1L]]]][[1L]])
  if (length(captured_warnings)) {
    result$message <- paste(unique(captured_warnings), collapse = " | ")
  }
  result
}

calculate_environment_target <- function(data) {
  usable <- copy(data[is.finite(days_after_sowing) & !is.na(genotype)])
  if (!nrow(usable)) {
    stop("Internal error: target group contains no usable cultivar observations.")
  }

  sowing_references <- unique(usable$reference_sowing_date)
  sowing_references <- sowing_references[!is.na(sowing_references)]
  if (length(sowing_references) != 1L) {
    stop(
      "Internal error: expected one environment-level reference sowing date, ",
      "but found ", length(sowing_references), "."
    )
  }
  sowing_reference <- sowing_references[[1L]]

  cultivar_means <- usable[, .(
    cultivar_mean_das = mean(days_after_sowing),
    n_cultivar_plots = .N
  ), by = genotype]

  simple_estimate <- as.numeric(mean(cultivar_means$cultivar_mean_das))
  simple_se <- if (nrow(cultivar_means) >= 2L) {
    safe_sd(cultivar_means$cultivar_mean_das) / sqrt(nrow(cultivar_means))
  } else {
    NA_real_
  }
  simple_critical <- if (nrow(cultivar_means) >= 2L) {
    stats::qt(0.975, df = nrow(cultivar_means) - 1L)
  } else {
    NA_real_
  }

  paper <- fit_paper_target(usable)
  paper_available <- paper$status %in% c(
    "paper_model_ok",
    "paper_model_singular"
  ) && is.finite(paper$estimate)

  target_estimate <- if (paper_available) paper$estimate else simple_estimate
  target_se <- if (paper_available) paper$standard_error else simple_se
  target_lower <- if (paper_available) {
    paper$lower_95
  } else if (is.finite(simple_se) && is.finite(simple_critical)) {
    simple_estimate - simple_critical * simple_se
  } else {
    NA_real_
  }
  target_upper <- if (paper_available) {
    paper$upper_95
  } else if (is.finite(simple_se) && is.finite(simple_critical)) {
    simple_estimate + simple_critical * simple_se
  } else {
    NA_real_
  }

  target_source <- if (paper_available) {
    "paper_BLUE_equal_cultivar_average"
  } else {
    "equal_weight_cultivar_mean_fallback"
  }

  data.table(
    target_available = TRUE,
    observed_environment_days_after_sowing = as.numeric(target_estimate),
    observed_environment_milestone_date = sowing_reference +
      as.integer(round(target_estimate)),
    reference_sowing_date = sowing_reference,
    reference_treatment = REFERENCE_TREATMENT,
    target_source = target_source,
    target_standard_error_days = as.numeric(target_se),
    target_lower_95_das = as.numeric(target_lower),
    target_upper_95_das = as.numeric(target_upper),
    adjustment_singular = paper$singular,
    paper_model_status = paper$status,
    paper_model_message = paper$message,
    n_plots = nrow(usable),
    n_genotypes = nrow(cultivar_means),
    n_blocks = uniqueN(usable$block, na.rm = TRUE),
    n_rows = uniqueN(usable$row, na.rm = TRUE),
    n_columns = uniqueN(usable$column, na.rm = TRUE),
    n_sowing_dates = uniqueN(usable$sowing_date, na.rm = TRUE),
    n_environment_sowing_dates = first(
      usable$n_environment_sowing_dates
    ),
    earliest_environment_sowing_date = first(
      usable$earliest_environment_sowing_date
    ),
    latest_environment_sowing_date = first(
      usable$latest_environment_sowing_date
    ),
    environment_sowing_range_days = first(
      usable$environment_sowing_range_days
    ),
    n_paper_model_rows = paper$n_model_rows,
    raw_plot_mean_das = as.numeric(mean(usable$days_after_sowing)),
    raw_plot_median_das = as.numeric(stats::median(usable$days_after_sowing)),
    plot_sd_days = safe_sd(usable$days_after_sowing),
    plot_iqr_days = safe_iqr(usable$days_after_sowing),
    equal_weight_cultivar_mean_das = simple_estimate,
    paper_minus_simple_days = if (paper_available) {
      as.numeric(paper$estimate - simple_estimate)
    } else {
      NA_real_
    }
  )
}

# ------------------------------------------------------------------------------
# 3. Load and standardise plot-level BRIWECS observations
# ------------------------------------------------------------------------------

if (!file.exists(INPUT_FILE)) {
  stop("BRIWECS input file does not exist: ", INPUT_FILE)
}

raw <- fread(
  INPUT_FILE,
  sep = ";",
  na.strings = c("", "NA", "NaN"),
  encoding = "UTF-8"
)
raw[, record_id := .I]

column_map <- list(
  Location = c("Location", "location"),
  sowing_year = c("Year", "year"),
  sowing_doy = c("Sowing_date", "sowing_DOY", "Sowing_DOY"),
  bbch10_doy = c("Emergence_date", "BBCH10"),
  bbch59_doy = c("BBCH59"),
  bbch87_doy = c("BBCH87"),
  genotype = c("BRISONr", "Genotype", "genotype"),
  treatment = c("Treatment", "treatment"),
  block = c("Block", "block"),
  row = c("Row", "row"),
  column = c("Column", "column")
)

resolved <- vapply(
  names(column_map),
  function(name) resolve_column(
    raw,
    column_map[[name]],
    name,
    required = !name %in% c("block", "row", "column")
  ),
  character(1)
)

optional_identifier <- function(name) {
  column <- resolved[[name]]
  if (is.na(column)) return(rep(NA_character_, nrow(raw)))
  clean_identifier(raw[[column]])
}

bri <- data.table(
  record_id = raw$record_id,
  Location = clean_identifier(raw[[resolved[["Location"]]]]),
  sowing_year = as.integer(numeric_locale(raw[[resolved[["sowing_year"]]]])),
  sowing_doy = as.integer(numeric_locale(raw[[resolved[["sowing_doy"]]]])),
  bbch10_doy = as.integer(numeric_locale(raw[[resolved[["bbch10_doy"]]]])),
  bbch59_doy = as.integer(numeric_locale(raw[[resolved[["bbch59_doy"]]]])),
  bbch87_doy = as.integer(numeric_locale(raw[[resolved[["bbch87_doy"]]]])),
  genotype = clean_identifier(raw[[resolved[["genotype"]]]]),
  treatment = clean_identifier(raw[[resolved[["treatment"]]]]),
  block = optional_identifier("block"),
  row = optional_identifier("row"),
  column = optional_identifier("column")
)

bri <- bri[
  !is.na(Location) &
    !is.na(sowing_year) &
    !is.na(sowing_doy) &
    treatment == REFERENCE_TREATMENT
]
if (!nrow(bri)) {
  stop("No records remain for reference treatment: ", REFERENCE_TREATMENT)
}
if (bri[, any(sowing_doy < 1L | sowing_doy > 366L)]) {
  stop("Sowing day-of-year lies outside 1-366.")
}

bri[, `:=`(
  SeasonYear = sowing_year + 1L,
  sowing_date = date_from_year_doy(sowing_year, sowing_doy),
  BBCH10_date = date_from_year_doy(sowing_year, bbch10_doy),
  BBCH59_date = date_from_year_doy(sowing_year + 1L, bbch59_doy),
  BBCH87_date = date_from_year_doy(sowing_year + 1L, bbch87_doy)
)]
bri[, environment_id := paste(Location, SeasonYear, sep = "__")]

# Define the sowing-date reference once per environment, before milestone-
# specific missing observations can change the contributing plot set. The same
# reference date is therefore used for BBCH10, BBCH59 and BBCH87. The median is
# calculated across all records of the selected reference treatment.
environment_sowing_reference <- bri[
  !is.na(sowing_date),
  .(
    reference_sowing_date = median_date(sowing_date),
    n_environment_sowing_dates = uniqueN(sowing_date),
    earliest_environment_sowing_date = min(sowing_date),
    latest_environment_sowing_date = max(sowing_date),
    environment_sowing_range_days = as.integer(
      max(sowing_date) - min(sowing_date)
    )
  ),
  by = .(environment_id, Location, SeasonYear)
]

if (environment_sowing_reference[, anyNA(reference_sowing_date)]) {
  stop("At least one environment has no valid reference sowing date.")
}

milestone_columns <- c("BBCH10_date", "BBCH59_date", "BBCH87_date")
milestones <- melt(
  bri,
  id.vars = setdiff(names(bri), milestone_columns),
  measure.vars = milestone_columns,
  variable.name = "milestone",
  value.name = "observed_date",
  variable.factor = FALSE
)
milestones[, milestone := sub("_date$", "", milestone)]
milestones <- milestones[!is.na(observed_date) & !is.na(sowing_date)]
milestones[, days_after_sowing := as.numeric(observed_date - sowing_date)]
milestones <- merge(
  milestones,
  environment_sowing_reference,
  by = c("environment_id", "Location", "SeasonYear"),
  all.x = TRUE,
  sort = FALSE
)

if (milestones[, anyNA(reference_sowing_date)]) {
  stop("An environment milestone could not be assigned a reference sowing date.")
}

n_missing_genotype_records <- milestones[is.na(genotype), .N]
if (n_missing_genotype_records) {
  message(
    "Excluding ", n_missing_genotype_records,
    " plot milestone records without a cultivar identifier."
  )
}
milestones <- milestones[!is.na(genotype)]
if (!nrow(milestones)) {
  stop("No milestone observations with a cultivar identifier remain.")
}

invalid <- milestones[
  !is.finite(days_after_sowing) |
    days_after_sowing < 0 |
    days_after_sowing > MAX_DAYS_AFTER_SOWING
]
if (nrow(invalid)) {
  stop(
    "Invalid milestone chronology for ", nrow(invalid),
    " plot records. Check dates before calculating targets."
  )
}

# ------------------------------------------------------------------------------
# 4. Calculate and validate one target per environment x milestone
# ------------------------------------------------------------------------------

targets <- milestones[, calculate_environment_target(.SD), by = .(
  environment_id,
  Location,
  SeasonYear,
  milestone
)]
targets[, BBCH := as.integer(sub("^BBCH", "", milestone))]

inconsistent_sowing_references <- targets[, .(
  n_reference_sowing_dates = uniqueN(reference_sowing_date, na.rm = TRUE)
), by = .(
  environment_id,
  Location,
  SeasonYear
)][n_reference_sowing_dates != 1L]

if (nrow(inconsistent_sowing_references)) {
  stop(
    "Reference sowing date is not unique within environment(s): ",
    paste(inconsistent_sowing_references$environment_id, collapse = ", ")
  )
}

setcolorder(targets, c(
  "environment_id",
  "Location",
  "SeasonYear",
  "milestone",
  "BBCH",
  "reference_treatment",
  "observed_environment_milestone_date",
  "observed_environment_days_after_sowing",
  "target_source",
  "target_standard_error_days",
  "target_lower_95_das",
  "target_upper_95_das",
  setdiff(names(targets), c(
    "environment_id",
    "Location",
    "SeasonYear",
    "milestone",
    "BBCH",
    "reference_treatment",
    "observed_environment_milestone_date",
    "observed_environment_days_after_sowing",
    "target_source",
    "target_standard_error_days",
    "target_lower_95_das",
    "target_upper_95_das"
  ))
))
setorder(targets, Location, SeasonYear, BBCH)

duplicate_targets <- targets[, .N, by = .(
  environment_id,
  milestone
)][N != 1L]
if (nrow(duplicate_targets)) {
  stop("Environment targets are not unique by environment and milestone.")
}

chronology <- dcast(
  targets[target_available == TRUE],
  environment_id + Location + SeasonYear ~ milestone,
  value.var = "observed_environment_days_after_sowing"
)
for (column in c("BBCH10", "BBCH59", "BBCH87")) {
  if (!column %in% names(chronology)) chronology[, (column) := NA_real_]
}
invalid_chronology <- chronology[
  (!is.na(BBCH10) & !is.na(BBCH59) & BBCH10 > BBCH59) |
    (!is.na(BBCH59) & !is.na(BBCH87) & BBCH59 > BBCH87) |
    (!is.na(BBCH10) & !is.na(BBCH87) & BBCH10 > BBCH87)
]
if (nrow(invalid_chronology)) {
  stop(
    "Aggregated milestone chronology is invalid for: ",
    paste(invalid_chronology$environment_id, collapse = ", ")
  )
}

output_file <- file.path(
  OUTPUT_DIR,
  "briwecs_environment_milestone_targets_paper_style.csv"
)
fwrite(targets, output_file, na = "NA")

message("Reference treatment: ", REFERENCE_TREATMENT)
message("Environment x milestone targets: ", nrow(targets))
print(targets[, .N, by = .(milestone, target_source)][
  order(milestone, target_source)
])
message("Saved: ", output_file)

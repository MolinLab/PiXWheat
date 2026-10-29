#!/usr/bin/env Rscript

# ==============================================================================
# EXTRACT HYRAS-DE WEATHER FOR THE SIX BRIWECS LOCATIONS, 2014-2022
# ==============================================================================
#
# The script reads the already downloaded HYRAS-DE v6-1 daily NetCDF files for
# minimum temperature, maximum temperature and precipitation, extracts the
# nearest 1-km grid cell for each BRIWECS location, and writes one combined
# semicolon-separated file.
#
# Output columns follow the weather contract of the BRIWECS validation script.
# Global.Radiation_MJ_m2_d is deliberately NA: radiation predictors are built
# from ERA5-Land ssrd in the validation pipeline, not from this HYRAS file.
#
# Required input layout (the layout produced by get_dwd_data.sh):
#   <HYRAS_DIR>/tasmin/tasmin_hyras_1_2014_v6-1_de.nc
#   <HYRAS_DIR>/tasmax/tasmax_hyras_1_2014_v6-1_de.nc
#   <HYRAS_DIR>/pr/pr_hyras_1_2014_v6-1_de.nc
#   ... through 2022
#
# Run from the repository root, or set PHENOLOGY_PROJECT_DIR.
# The BRIWECS and HYRAS input directories can be overridden separately.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(terra)
})

# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------

PROJECT_DIR <- Sys.getenv(
  "PHENOLOGY_PROJECT_DIR",
  unset = "."
)
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = TRUE)

BRIWECS_DIR <- Sys.getenv(
  "PHENOLOGY_BRIWECS_DIR",
  unset = file.path(PROJECT_DIR, "Data", "BRIWECS")
)

LOCATION_FILE <- Sys.getenv(
  "PHENOLOGY_BRIWECS_LOCATION_FILE",
  unset = file.path(BRIWECS_DIR, "BRIWECS_location_coordinates.csv")
)

HYRAS_DIR <- Sys.getenv(
  "PHENOLOGY_BRIWECS_HYRAS_DIR",
  unset = file.path(PROJECT_DIR, "Data", "DWD")
)

OUTPUT_FILE <- Sys.getenv(
  "PHENOLOGY_BRIWECS_HYRAS_OUTPUT_FILE",
  unset = file.path(
    BRIWECS_DIR,
    "DWD_HYRAS_2014_2022_all_locations.csv"
  )
)

START_YEAR <- 2014L
END_YEAR <- 2022L
EXPECTED_N_LOCATIONS <- 6L

VARIABLE_SPECS <- list(
  TN = list(
    directory = "tasmin",
    filename = "tasmin_hyras_1_%d_v6-1_de.nc"
  ),
  TX = list(
    directory = "tasmax",
    filename = "tasmax_hyras_1_%d_v6-1_de.nc"
  ),
  RR = list(
    directory = "pr",
    filename = "pr_hyras_1_%d_v6-1_de.nc"
  )
)

# ------------------------------------------------------------------------------
# 2. Small helpers
# ------------------------------------------------------------------------------

resolve_column <- function(data, candidates, label) {
  hits <- candidates[candidates %in% names(data)]
  if (!length(hits)) {
    stop(
      "Could not identify ", label, ". Expected one of: ",
      paste(candidates, collapse = ", ")
    )
  }
  hits[[1L]]
}

numeric_locale <- function(x) {
  as.numeric(gsub(",", ".", trimws(as.character(x)), fixed = TRUE))
}

assert_file <- function(path, label) {
  if (!file.exists(path)) stop(label, " does not exist: ", path)
}

hyras_file <- function(variable, year) {
  spec <- VARIABLE_SPECS[[variable]]
  file.path(
    HYRAS_DIR,
    spec$directory,
    sprintf(spec$filename, year)
  )
}

# Extract one variable-year block. terra::extract(method = "simple") uses the
# containing/nearest raster cell and reproduces the extraction rule used by the
# existing HYRAS workflow.
extract_variable_year <- function(variable, year, locations) {
  path <- hyras_file(variable, year)
  assert_file(path, paste(variable, year, "HYRAS file"))

  message("Extracting ", variable, " for ", year, ": ", basename(path))
  raster <- terra::rast(path)

  dates <- terra::time(raster)
  if (is.null(dates) || !length(dates) || all(is.na(dates))) {
    stop("No valid time metadata in ", path)
  }
  if (!inherits(dates, c("Date", "POSIXct", "POSIXt"))) {
    stop("Unexpected time metadata class in ", path)
  }
  dates <- as.Date(dates)

  if (length(dates) != terra::nlyr(raster)) {
    stop(
      "Time/layer mismatch in ", path, ": ",
      length(dates), " dates but ", terra::nlyr(raster), " layers."
    )
  }
  if (anyDuplicated(dates)) stop("Duplicate dates in ", path)
  if (any(as.integer(format(dates, "%Y")) != year)) {
    stop("The time dimension contains dates outside ", year, " in ", path)
  }

  points <- terra::vect(
    locations,
    geom = c("lon", "lat"),
    crs = "EPSG:4326"
  )
  points <- terra::project(points, terra::crs(raster))

  values <- terra::extract(raster, points, method = "simple")
  values <- as.data.table(values)

  if (ncol(values) - 1L != length(dates)) {
    stop(
      "Extraction/layer mismatch in ", path, ": ",
      ncol(values) - 1L, " value columns but ", length(dates), " dates."
    )
  }

  value_columns <- names(values)[-1L]
  setnames(values, value_columns, paste0("day_", seq_along(value_columns)))

  long <- melt(
    values,
    id.vars = "ID",
    variable.name = "day_index",
    value.name = variable,
    variable.factor = FALSE
  )
  long[, day_index := as.integer(sub("day_", "", day_index, fixed = TRUE))]
  long[, Date := dates[day_index]]

  site_key <- copy(locations[, .(ID = point_id, Location, lat, lon)])
  long <- merge(long, site_key, by = "ID", all.x = TRUE, sort = FALSE)
  long[, source_file := basename(path)]

  result <- long[, .(
    Location,
    Date,
    value = get(variable),
    lat,
    lon,
    source_file
  )]
  setnames(result, "value", variable)
  setorder(result, Location, Date)
  result[]
}

# ------------------------------------------------------------------------------
# 3. Read and validate the six BRIWECS locations
# ------------------------------------------------------------------------------

assert_file(LOCATION_FILE, "BRIWECS location file")
if (!dir.exists(HYRAS_DIR)) stop("HYRAS directory does not exist: ", HYRAS_DIR)

location_raw <- fread(
  LOCATION_FILE,
  na.strings = c("", "NA", "NaN"),
  encoding = "UTF-8"
)

location_col <- resolve_column(
  location_raw,
  c("Location", "location", "Standort", "standort"),
  "location name"
)
lat_col <- resolve_column(
  location_raw,
  c("Lat", "lat", "latitude", "Latitude"),
  "latitude"
)
lon_col <- resolve_column(
  location_raw,
  c("Long", "Lon", "lon", "longitude", "Longitude"),
  "longitude"
)

locations <- unique(location_raw[, .(
  Location = trimws(as.character(get(location_col))),
  lat = numeric_locale(get(lat_col)),
  lon = numeric_locale(get(lon_col))
)])
setorder(locations, Location)
locations[, point_id := .I]

if (nrow(locations) != EXPECTED_N_LOCATIONS) {
  stop(
    "Expected exactly ", EXPECTED_N_LOCATIONS,
    " BRIWECS locations, but found ", nrow(locations), "."
  )
}
if (locations[, anyNA(Location) | anyNA(lat) | anyNA(lon)]) {
  stop("The BRIWECS coordinate table contains missing values.")
}
if (anyDuplicated(locations$Location)) {
  stop("BRIWECS location names are not unique.")
}
if (locations[, any(lat < 47 | lat > 56 | lon < 5 | lon > 16)]) {
  stop(
    "At least one coordinate is outside the expected German domain. ",
    "Check longitude/latitude order and units."
  )
}

# Fail before extraction if any required annual file is absent.
required_files <- CJ(
  variable = names(VARIABLE_SPECS),
  year = seq.int(START_YEAR, END_YEAR)
)
required_files[, file := mapply(hyras_file, variable, year)]
missing_files <- required_files[!file.exists(file)]
if (nrow(missing_files)) {
  stop(
    "Missing ", nrow(missing_files), " required HYRAS files:\n",
    paste(missing_files$file, collapse = "\n")
  )
}

# ------------------------------------------------------------------------------
# 4. Extract, combine and audit the complete daily table
# ------------------------------------------------------------------------------

blocks <- lapply(names(VARIABLE_SPECS), function(variable) {
  result <- rbindlist(
    lapply(seq.int(START_YEAR, END_YEAR), function(year) {
      extract_variable_year(variable, year, locations)
    }),
    use.names = TRUE
  )
  setnames(result, "source_file", paste0("source_", variable, "_file"))
  result
})
names(blocks) <- names(VARIABLE_SPECS)

weather <- Reduce(
  function(left, right) {
    merge(
      left,
      right,
      by = c("Location", "Date", "lat", "lon"),
      all = TRUE,
      sort = FALSE
    )
  },
  blocks
)

setnames(
  weather,
  c("source_TN_file", "source_TX_file", "source_RR_file"),
  c("source_tasmin_file", "source_tasmax_file", "source_pr_file")
)

setorder(weather, Location, Date)

expected_dates <- seq.Date(
  as.Date(sprintf("%d-01-01", START_YEAR)),
  as.Date(sprintf("%d-12-31", END_YEAR)),
  by = "day"
)
expected_rows <- nrow(locations) * length(expected_dates)

if (nrow(weather) != expected_rows) {
  stop(
    "Expected ", expected_rows, " location-date rows but obtained ",
    nrow(weather), "."
  )
}
if (weather[, anyDuplicated(paste(Location, Date)) > 0L]) {
  stop("Duplicate Location-Date rows were produced.")
}
if (weather[, anyNA(TN) | anyNA(TX) | anyNA(RR)]) {
  missing_summary <- weather[, .(
    n_missing_TN = sum(is.na(TN)),
    n_missing_TX = sum(is.na(TX)),
    n_missing_RR = sum(is.na(RR))
  )]
  stop(
    "Missing HYRAS values remain after extraction: ",
    paste(names(missing_summary), unlist(missing_summary), collapse = ", ")
  )
}
if (weather[, any(TN > TX)]) stop("At least one day has TN > TX.")
if (weather[, any(TN < -50 | TX > 50)]) {
  stop("Extracted temperature lies outside the -50 to 50 degree C audit range.")
}
if (weather[, any(RR < 0 | RR > 500)]) {
  stop("Extracted precipitation lies outside the 0 to 500 mm/day audit range.")
}

weather[, `:=`(
  Temperature.Air.Mean_C = (TN + TX) / 2,
  Temperature.Air.Min_C = TN,
  Temperature.Air.Max_C = TX,
  Precipitation_mm_d = RR,
  # Compatibility column only. The validation script uses ERA5-Land ssrd.
  Global.Radiation_MJ_m2_d = NA_real_,
  weather_source = "HYRAS-DE v6-1"
)]

output <- weather[, .(
  Location,
  Date = format(Date, "%Y-%m-%d"),
  Temperature.Air.Mean_C,
  Temperature.Air.Min_C,
  Temperature.Air.Max_C,
  Precipitation_mm_d,
  Global.Radiation_MJ_m2_d,
  lat,
  lon,
  weather_source,
  source_tasmin_file,
  source_tasmax_file,
  source_pr_file
)]

dir.create(dirname(OUTPUT_FILE), recursive = TRUE, showWarnings = FALSE)
fwrite(output, OUTPUT_FILE, sep = ";", na = "NA")

message(
  "Finished. Wrote ", nrow(output), " daily rows for ",
  uniqueN(output$Location), " locations to:\n", OUTPUT_FILE
)
message(
  "Use with the validation script by setting:\n",
  "PHENOLOGY_BRIWECS_WEATHER_PATTERN=",
  "^DWD_HYRAS_2014_2022_all_locations[.]csv$"
)

# Extract ERA5-Land daily radiation, soil water, evaporation, and potential
# evaporation at all Austrian phenology locations.

library(readr)
library(dplyr)
library(terra)
library(data.table)
library(future.apply)

PROJECT_DIR <- normalizePath(
  Sys.getenv("PHENOLOGY_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
DATA_DIR <- file.path(PROJECT_DIR, "Data")
BASE_DIR <- file.path(DATA_DIR, "Era5")

# -----------------------------------------------------------------------------
# Phenology locations
# -----------------------------------------------------------------------------
phen_df <- read_csv(
  file.path(DATA_DIR, "PEP725", "pep_at_filtered_final_19_06.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    location_id = case_when(
      is.na(location_id) ~ as.character(interaction(lat, lon)),
      TRUE ~ as.character(location_id)
    ),
    location = coalesce(Location, location, Location_name)
  )

# Keep station identity and original WGS84 coordinates. 
coords <- phen_df %>%
  distinct(Stations_id,location_id, location, lon, lat) %>%
  transmute(
    Stations_id,
    location_id,
    location,
    lon_orig = as.numeric(lon),
    lat_orig = as.numeric(lat)
  ) %>%
  filter(!is.na(lon_orig), !is.na(lat_orig)) %>%
  arrange(location_id, lon_orig, lat_orig) %>%
  mutate(id = Stations_id) %>%
  select(id, location_id, location, lon_orig, lat_orig) %>%
  as.data.frame()

if (nrow(coords) == 0L) {
  stop("No valid phenology coordinates are available for extraction.")
}

location_coordinate_versions <- coords %>%
  count(location_id, name = "n_coordinate_versions") %>%
  filter(n_coordinate_versions > 1L)

if (nrow(location_coordinate_versions) > 0L) {
  stop(
    "Some Austrian location_id values have multiple coordinate pairs. ",
    "Resolve these before weather extraction; otherwise downstream joins by ",
    "location_id would be ambiguous. Example IDs: ",
    paste(head(location_coordinate_versions$location_id, 10), collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# ERA5-Land input and output
# -----------------------------------------------------------------------------
base_dir <- BASE_DIR

# Annual daily NetCDF files are produced by
# extract_radiation_era5_with_soilwater_new.py. The four SWVL variables are UTC
# daily means of the hourly soil-water states; ssrd, e, and pev are daily UTC
# accumulations.
input_dir <- file.path(base_dir, "era5land_daily_yearly")
input_files <- sort(list.files(
  input_dir,
  pattern = "^era5land_daily_[0-9]{4}\\.nc$",
  full.names = TRUE
))

variables <- c("ssrd", "swvl1", "swvl2", "swvl3", "swvl4", "e", "pev")

outdir <- file.path(
  base_dir,
  "era5land_yearly",
  "temp_extract_geosphere"
)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

if (length(input_files) == 0L) {
  stop("No annual ERA5-Land input files found in: ", input_dir)
}

# -----------------------------------------------------------------------------
# Extract one named NetCDF variable and write one CSV per calendar year
# -----------------------------------------------------------------------------
process_variable <- function(varname, f, outdir, coords) {
  cat("Processing variable:", varname, "from", basename(f), "\n")

  # Select the named NetCDF subdataset. Reading the entire multi-variable file
  # with rast(f) would combine variables and make the time-to-layer mapping
  # ambiguous.
  r <- tryCatch(
    terra::rast(f, subds = varname),
    error = function(e) {
      stop(
        "Could not read NetCDF variable '", varname, "' from ", basename(f),
        ": ", conditionMessage(e)
      )
    }
  )

  if (terra::nlyr(r) == 0L) {
    stop("Variable '", varname, "' contains no raster layers.")
  }

  # Rebuild points from original coordinates, then project to the raster CRS.
  pts <- terra::vect(
    coords,
    geom = c("lon_orig", "lat_orig"),
    crs = "EPSG:4326"
  )
  pts <- terra::project(pts, terra::crs(r))

  # Extract and validate the daily time axis.
  time_vec <- terra::time(r)

  if (is.null(time_vec) || length(time_vec) == 0L || all(is.na(time_vec))) {
    stop(
      "No valid time metadata found for variable '", varname,
      "' in ", basename(f)
    )
  }

  if (!inherits(time_vec, c("Date", "POSIXct", "POSIXt"))) {
    stop(
      "Unexpected time format for variable '", varname,
      "' in ", basename(f), ": ", paste(class(time_vec), collapse = ", ")
    )
  }

  time_vec <- as.Date(time_vec)

  if (length(time_vec) != terra::nlyr(r)) {
    stop(
      "Time dimension mismatch for variable '", varname, "': nlyr(r) = ",
      terra::nlyr(r), ", length(time_vec) = ", length(time_vec)
    )
  }

  if (anyDuplicated(time_vec)) {
    stop("Duplicate dates found for variable '", varname, "'.")
  }

  # Bilinear interpolation is appropriate for these continuous ERA5-Land
  # variables and avoids unnecessary missing values near coastlines.
  vals <- terra::extract(r, pts, method = "bilinear")

  if ((ncol(vals) - 1L) != length(time_vec)) {
    stop(
      "Layer mismatch for variable '", varname,
      "': extracted columns = ", ncol(vals) - 1L,
      ", time steps = ", length(time_vec)
    )
  }

  vals_dt <- as.data.table(vals)
  setnames(
    vals_dt,
    old = names(vals_dt)[-1L],
    new = paste0("layer_", seq_along(time_vec))
  )

  # One row per station and day.
  long_dt <- melt(
    vals_dt,
    id.vars = "ID",
    variable.name = "layer",
    value.name = varname,
    variable.factor = FALSE
  )

  long_dt[, layer := as.integer(sub("layer_", "", layer, fixed = TRUE))]

  time_dt <- data.table(
    layer = seq_along(time_vec),
    Date = time_vec,
    DOY = as.integer(format(time_vec, "%j")),
    year = as.integer(format(time_vec, "%Y"))
  )

  coords_dt <- as.data.table(coords)
  setnames(coords_dt, "id", "ID")

  long_dt <- merge(
    long_dt,
    time_dt,
    by = "layer",
    all.x = TRUE,
    sort = FALSE
  )
  long_dt <- merge(
    long_dt,
    coords_dt,
    by = "ID",
    all.x = TRUE,
    sort = FALSE
  )

  # Keep the native ERA5-Land variables and add analysis-ready unit
  # conversions without overwriting the source values.
  if (varname == "ssrd") {
    long_dt[, ssrd_MJ_m2 := ssrd / 1e6]
  } else if (varname == "e") {
    # ERA5 uses downward-positive fluxes: evaporation is normally negative.
    long_dt[, evaporation_mm := -e * 1000]
  } else if (varname == "pev") {
    long_dt[, potential_evaporation_mm := -pev * 1000]
  }

  setorder(long_dt, ID, Date)

  # Each NetCDF contains one calendar year. The CSV names retain compatibility
  # with read_weather_year() in the downstream cumulative-phenology scripts.
  years <- sort(unique(long_dt$year))

  if (length(years) != 1L) {
    stop(
      "Expected exactly one calendar year in ", basename(f),
      ", but found: ", paste(years, collapse = ", ")
    )
  }

  for (yr in years) {
    year_dt <- long_dt[year == yr]
    output_file <- file.path(outdir, paste0(varname, "_", yr, ".csv"))

    fwrite(year_dt, file = output_file)

    expected_rows <- nrow(coords) * sum(format(time_vec, "%Y") == yr)
    if (nrow(year_dt) != expected_rows) {
      stop(
        "Unexpected row count in ", basename(output_file), ": observed ",
        nrow(year_dt), ", expected ", expected_rows
      )
    }
  }

  missing_n <- sum(is.na(long_dt[[varname]]))
  cat(
    "Finished", varname, "-", nrow(long_dt), "station-days;",
    missing_n, "missing values\n"
  )

  rm(r, pts, time_vec, vals, vals_dt, long_dt, time_dt, coords_dt)
  gc()

  invisible(NULL)
}

# Process two annual files at a time. Variables within a file are processed
# sequentially to avoid opening seven subdatasets from the same file at once.
plan(multisession, workers = 3)

future_lapply(
  input_files,
  function(year_file) {
    for (varname in variables) {
      process_variable(
        varname = varname,
        f = year_file,
        outdir = outdir,
        coords = coords
      )
    }
    invisible(NULL)
  },
  future.packages = c("terra", "data.table"),
  future.seed = TRUE
)

plan(sequential)

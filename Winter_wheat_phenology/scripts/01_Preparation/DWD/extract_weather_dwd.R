# Extract HYRAS daily weather at DWD phenology stations.
# Run from the repository root or set PHENOLOGY_PROJECT_DIR.
library(readr)
library(dplyr)
library(terra)
library(data.table)
library(future)
library(future.apply)

PROJECT_DIR <- normalizePath(
  Sys.getenv("PHENOLOGY_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
WORK_DIR <- file.path(PROJECT_DIR, "Data", "DWD")
if (!dir.exists(WORK_DIR)) stop("DWD input directory not found: ", WORK_DIR)
setwd(WORK_DIR)

phen_df <- read_csv(
  file.path(WORK_DIR, "DWD_phenology_16_06_final.csv"),
  show_col_types = FALSE
)

# Keep station identity and WGS84 coordinates for the site-level extraction.
coords <- phen_df %>%
  distinct(Stations_id, lon, lat) %>%
  mutate(
    lon_orig = as.numeric(lon),
    lat_orig = as.numeric(lat)
  ) %>%
  select(Stations_id, lon_orig, lat_orig) %>%
  arrange(Stations_id) %>%
  mutate(id = row_number()) %>%
  select(id, Stations_id, lon_orig, lat_orig)

coords <- as.data.frame(coords)

files_tmax <- list.files("tasmax", pattern = "tasmax.*\\.nc$", full.names = TRUE)
files_tmin <- list.files("tasmin", pattern = "tasmin.*\\.nc$", full.names = TRUE)
files_prec <- list.files("pr", pattern = "pr.*\\.nc$", full.names = TRUE)
files_rad  <- list.files("rsds", pattern = "rsds.*\\.nc$", full.names = TRUE)

if (any(lengths(list(
  tasmax = files_tmax,
  tasmin = files_tmin,
  precipitation = files_prec,
  radiation = files_rad
)) == 0L)) {
  stop("At least one HYRAS variable has no NetCDF input files in ", WORK_DIR)
}

# Output folders
dir.create("tasmax/temp_extract", showWarnings = FALSE, recursive = TRUE)
dir.create("tasmin/temp_extract", showWarnings = FALSE, recursive = TRUE)
dir.create("pr/temp_extract", showWarnings = FALSE, recursive = TRUE)
dir.create("rsds/temp_extract", showWarnings = FALSE, recursive = TRUE)

process_year <- function(f, varname, outdir) {
  
  cat("Processing:", basename(f), "\n")
  
  r <- rast(f)
  
  # Rebuild points from original coordinates, then project to raster CRS
  pts <- vect(coords, geom = c("lon_orig", "lat_orig"), crs = "EPSG:4326")
  pts <- project(pts, crs(r))
  
  # Extract raster time metadata
  time_vec <- terra::time(r)
  if (is.null(time_vec) || all(is.na(time_vec))) {
    stop("No valid time metadata found in raster: ", basename(f))
  }
  
  if (!inherits(time_vec, c("Date", "POSIXct", "POSIXt"))) {
    stop("Unexpected time format in raster: ", basename(f))
  }
  
  time_vec <- as.Date(time_vec)
  
  if (length(time_vec) != nlyr(r)) {
    stop(
      "Time dimension mismatch in file ", basename(f),
      ": nlyr(r) = ", nlyr(r),
      ", length(time_vec) = ", length(time_vec)
    )
  }
  
  # Extract the containing grid cell using the original simple method.
  vals <- terra::extract(r, pts, method = "simple")
  
  # First column is point ID; remaining columns are layers/time steps
  if ((ncol(vals) - 1) != length(time_vec)) {
    stop(
      "Layer mismatch in file ", basename(f),
      ": extracted columns = ", ncol(vals) - 1,
      ", time steps = ", length(time_vec)
    )
  }
  
  vals_dt <- as.data.table(vals)
  
  # Standardize layer names
  setnames(vals_dt, old = names(vals_dt)[-1], new = paste0("layer_", seq_along(time_vec)))
  
  # Long format: one row per station-point x day
  long_dt <- melt(
    vals_dt,
    id.vars = "ID",
    variable.name = "layer",
    value.name = varname
  )
  
  long_dt[, layer := as.integer(sub("layer_", "", layer))]
  
  time_dt <- data.table(
    layer = seq_along(time_vec),
    Date  = time_vec,
    DOY   = as.integer(format(time_vec, "%j")),
    year  = as.integer(format(time_vec, "%Y"))
  )
  
  coords_dt <- as.data.table(coords)
  setnames(coords_dt, "id", "ID")
  
  long_dt <- merge(long_dt, time_dt, by = "layer", all.x = TRUE, sort = FALSE)
  long_dt <- merge(long_dt, coords_dt, by = "ID", all.x = TRUE, sort = FALSE)
  
  setorder(long_dt, ID, Date)
  
  # Use year tag from the actual raster time axis
  yrs <- unique(format(time_vec, "%Y"))
  year_tag <- if (length(yrs) == 1) {
    yrs
  } else {
    paste0(min(yrs), "_", max(yrs))
  }
  
  fwrite(
    long_dt,
    file = file.path(outdir, paste0(varname, "_", year_tag, ".csv"))
  )
  
  rm(r, pts, time_vec, vals, vals_dt, long_dt, time_dt, coords_dt)
  gc()
  
  invisible(NULL)
}

plan(multisession, workers = 5)

future_lapply(
  files_tmax,
  process_year,
  varname = "tmax",
  outdir = "tasmax/temp_extract/",
  future.packages = c("terra", "data.table"),
  future.globals = c("coords"),
  future.seed = TRUE
)

future_lapply(
  files_tmin,
  process_year,
  varname = "tmin",
  outdir = "tasmin/temp_extract/",
  future.packages = c("terra", "data.table"),
  future.globals = c("coords"),
  future.seed = TRUE
)

future_lapply(
  files_prec,
  process_year,
  varname = "prec",
  outdir = "pr/temp_extract/",
  future.packages = c("terra", "data.table"),
  future.globals = c("coords"),
  future.seed = TRUE
)

future_lapply(
  files_rad,
  process_year,
  varname = "rad",
  outdir = "rsds/temp_extract/",
  future.packages = c("terra", "data.table"),
  future.globals = c("coords"),
  future.seed = TRUE
)

plan(sequential)

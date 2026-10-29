# Extract SPARTACUS daily weather at Austrian phenology locations.
library(dplyr)
library(readr)
library(terra)
library(data.table)
library(future)
library(future.apply)

PROJECT_DIR <- normalizePath(
  Sys.getenv("PHENOLOGY_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
WORK_DIR <- file.path(PROJECT_DIR, "Data", "GeoSphere")
if (!dir.exists(WORK_DIR)) stop("SPARTACUS input directory not found: ", WORK_DIR)
setwd(WORK_DIR)

phen_df <- read_csv(
  file.path(PROJECT_DIR, "Data", "PEP725", "pep_at_filtered_final_19_06.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    location_id = case_when(
      is.na(location_id) ~ interaction(lat, lon),
      TRUE ~ location_id
    ),
    location = coalesce(Location, location, Location_name)
  )

# Keep one row per location and its coordinates
coords <- phen_df %>%
  distinct(Stations_id,location_id, location, lon, lat) %>%
  mutate(
    lon_orig = as.numeric(lon),
    lat_orig = as.numeric(lat)
  ) %>%
  select(location_id,Stations_id, location, lon_orig, lat_orig) %>%
  arrange(location_id, location) %>%
  mutate(id = Stations_id) %>%
  select(id, location_id, location, lon_orig, lat_orig)

coords <- as.data.frame(coords)

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
  
  # Extract values
  vals <- terra::extract(r, pts, method = "simple")
  
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
  
  # Long format
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
  
  # File tag
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

files_tmax <- list.files("TX", pattern = "\\.nc$", full.names = TRUE)
files_tmin <- list.files("TN", pattern = "\\.nc$", full.names = TRUE)
files_prec <- list.files("RR", pattern = "\\.nc$", full.names = TRUE)
files_rad  <- list.files("SA", pattern = "\\.nc$", full.names = TRUE)

dir.create("TX/temp_extract",     showWarnings = FALSE, recursive = TRUE)
dir.create("TN/temp_extract",     showWarnings = FALSE, recursive = TRUE)
dir.create("RR/temp_extract",   showWarnings = FALSE, recursive = TRUE)
dir.create("SA/temp_extract", showWarnings = FALSE, recursive = TRUE)

plan(multisession, workers = 5)

future_lapply(
  files_tmax,
  process_year,
  varname = "tmax",
  outdir = "TX/temp_extract/",
  future.packages = c("terra", "data.table"),
  future.globals = c("coords"),
  future.seed = TRUE
)

future_lapply(
  files_tmin,
  process_year,
  varname = "tmin",
  outdir = "TN/temp_extract/",
  future.packages = c("terra", "data.table"),
  future.globals = c("coords"),
  future.seed = TRUE
)

future_lapply(
  files_prec,
  process_year,
  varname = "prec",
  outdir = "RR/temp_extract/",
  future.packages = c("terra", "data.table"),
  future.globals = c("coords"),
  future.seed = TRUE
)

future_lapply(
  files_rad,
  process_year,
  varname = "rad",
  outdir = "SA/temp_extract/",
  future.packages = c("terra", "data.table"),
  future.globals = c("coords"),
  future.seed = TRUE
)

plan(sequential)


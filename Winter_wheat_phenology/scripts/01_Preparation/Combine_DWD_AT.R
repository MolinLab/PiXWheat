# =============================================================================
# Combine German and Austrian winter-wheat phenology data
# =============================================================================
#
# Inputs are the leakage-safe observation-level outputs of the two country-
# specific GDD/process scripts. Full-season summaries are deliberately not
# reconstructed here. Spatial attributes are added once per physical site.
# =============================================================================


# ----------------------------------------------------------------------------
# 1. Configuration
# ----------------------------------------------------------------------------

PROJECT_DIR <- Sys.getenv(
  "PHENOLOGY_PROJECT_DIR",
  unset = "/data/analysis3/milekj/UseCase_WW02_2024"
)

GERMANY_FILE <- Sys.getenv(
  "PHENOLOGY_GERMANY_FILE",
  unset = file.path(PROJECT_DIR, "DWD", "DWD_phenology_GDD_era5_test1.csv")
)

AUSTRIA_FILE <- Sys.getenv(
  "PHENOLOGY_AUSTRIA_FILE",
  unset = file.path(PROJECT_DIR, "geosphere", "at_phenology_GDD_era5_test.csv")
)

ENV_ZONE_RASTER <- Sys.getenv(
  "PHENOLOGY_ENV_ZONE_RASTER",
  unset = file.path(
    PROJECT_DIR,
    "eea_r_3035_1_km_env-zones_p_2018_v01_r00",
    "eea_r_3035_1_km_env-zones_p_2018_v01_r00.tif"
  )
)

SPATIAL_CACHE_DIR <- file.path(PROJECT_DIR, "spatial_cache")
NUTS_CACHE_FILE <- file.path(SPATIAL_CACHE_DIR, "nuts2_at_de_2024.rds")
OUTPUT_FILE <- file.path(PROJECT_DIR, "pheno_df_CV.csv")
QC_FILE <- file.path(PROJECT_DIR, "pheno_combined_qc_summary.csv")


# ----------------------------------------------------------------------------
# 2. Packages and helpers
# ----------------------------------------------------------------------------

required_packages <- c("dplyr", "readr", "tibble", "sf", "terra", "giscoR")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install these packages before running the combination script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(sf)
})

assert_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

assert_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(
      label, " is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

standardize_germany <- function(data) {
  required <- c(
    "Stations_id", "Stationsname", "traj_id", "phen_date",
    "year", "season_year", "DOY", "BBCH", "lat_orig", "lon_orig",
    "sow_doy", "daylength", "cumDTU_beta"
  )
  assert_columns(data, required, "German phenology data")

  data %>%
    mutate(
      country = "DE",
      location_name = as.character(Stationsname),
      location_id = if_else(
        !is.na(Stations_id),
        paste0("DE_", as.character(Stations_id)),
        NA_character_
      ),
      Location = location_id,
      #site_key = if_else(
       # !is.na(site_key),
        #paste0("DE_", as.character(site_key)),
        #location_id
      #),
      traj_id = if_else(
        !is.na(traj_id),
        paste0("DE_", as.character(traj_id)),
        NA_character_
      ),
      env_id = traj_id,
      Date = as.Date(phen_date),
      SeasonYear = as.integer(season_year),
      sowing_DOY = as.integer(sow_doy),
      Daylength = as.numeric(daylength),
      lat = as.numeric(lat_orig),
      lon = as.numeric(lon_orig)
    )
}

standardize_austria <- function(data) {
  required <- c(
    "Stations_id", "Stationsname", "traj_id", "phen_date",
    "year", "season_year", "DOY", "BBCH", "lat_orig", "lon_orig",
    "sow_doy", "daylength", "cumDTU_beta"
  )
  assert_columns(data, required, "Austrian phenology data")

  data %>%
    mutate(
      country = "AT",
      location_name = as.character(Stationsname),
      location_id = if_else(
        !is.na(Stations_id),
        paste0("AT_", as.character(Stations_id)),
        NA_character_
      ),
      Location = location_id,
      traj_id = if_else(
        !is.na(traj_id),
        paste0("AT_", as.character(traj_id)),
        NA_character_
      ),
      env_id = traj_id,
      Date = as.Date(phen_date),
      SeasonYear = as.integer(season_year),
      sowing_DOY = as.integer(sow_doy),
      Daylength = as.numeric(daylength),
      lat = as.numeric(lat_orig),
      lon = as.numeric(lon_orig)
    )
}


# ----------------------------------------------------------------------------
# 3. Read and standardize country data
# ----------------------------------------------------------------------------

assert_file_exists(GERMANY_FILE, "German input")
assert_file_exists(AUSTRIA_FILE, "Austrian input")
assert_file_exists(ENV_ZONE_RASTER, "Environmental-zone raster")

germany <- read_csv(GERMANY_FILE, show_col_types = FALSE) %>%
  standardize_germany() %>%  mutate(sowing_DOY = if_else(BBCH == 0, DOY, NA_real_)) %>% group_by(traj_id, year) %>%
  mutate(
    sowing_DOY = suppressWarnings(min(sowing_DOY, na.rm = TRUE)),
    sowing_DOY = if_else(is.infinite(sowing_DOY), NA_real_, sowing_DOY)
  ) %>%
  ungroup()

austria <- read_csv(AUSTRIA_FILE, show_col_types = FALSE) %>%group_by (location, lat_orig, lon_orig) %>%
  mutate(Stations_id = cur_group_id()) %>%
  ungroup() %>% mutate(Stationsname=location) %>%
  standardize_austria()

pheno_combined <- bind_rows(germany, austria) %>%
  filter(
    #!is.na(cumGDD_from_sow_Soltani),
    !is.na(BBCH),
    BBCH <= 90,
    !is.na(Date),
    #!is.na(Location),
    #!is.na(env_id),
    !is.na(sowing_DOY),
    between(sowing_DOY, 1L, 366L),
    !is.na(lat),
    !is.na(lon)
  )

if (!nrow(pheno_combined)) {
  stop("No observations remain after the common filters.", call. = FALSE)
}

duplicate_observations <- pheno_combined %>%
  count(country, env_id, Date, BBCH, name = "n") %>%
  filter(n > 1L)

if (nrow(duplicate_observations)) {
  stop(
    "Duplicate country/env_id/date/BBCH observations found. Inspect ",
    "duplicate_observations before combining.",
    call. = FALSE
  )
}


# ----------------------------------------------------------------------------
# 4. Build one spatial record per physical site
# ----------------------------------------------------------------------------

site_coordinate_check <- pheno_combined %>%
  distinct(Location, location_id, location_name, country, lon, lat) %>%
  count(Location, name = "n_coordinate_pairs") %>%
  filter(n_coordinate_pairs > 1L)

if (nrow(site_coordinate_check)) {
  stop(
    "At least one physical Location maps to multiple coordinate pairs. ",
    "Inspect site_coordinate_check.",
    call. = FALSE
  )
}

sites <- pheno_combined %>%
  distinct(Stations_id,Location, location_id, location_name, country, lon, lat)

site_points <- terra::vect(
  sites,
  geom = c("lon", "lat"),
  crs = "EPSG:4326"
)
environmental_zones <- terra::rast(ENV_ZONE_RASTER)
site_points_projected <- terra::project(site_points, environmental_zones)
zone_values <- terra::extract(terra::as.int(environmental_zones), site_points_projected)

if (ncol(zone_values) < 2L) {
  stop("Environmental-zone extraction returned no raster values.", call. = FALSE)
}

zone_lookup <- tibble::tribble(
  ~zone_code, ~env_zone,
  1L,  "ALN",
  2L,  "BOR",
  3L,  "NEM",
  4L,  "ATN",
  5L,  "ALS",
  6L,  "CON",
  7L,  "ATC",
  8L,  "PAN",
  9L,  "LUS",
  10L, "ANA",
  11L, "MDM",
  12L, "MDN",
  13L, "MDS"
)

site_spatial <- sites %>%
  mutate(zone_code = as.integer(zone_values[[2L]])) %>%
  left_join(zone_lookup, by = "zone_code")

if (anyNA(site_spatial$env_zone)) {
  missing_sites <- site_spatial %>% filter(is.na(env_zone)) %>% pull(Location)
  stop(
    "Environmental zone is missing for: ",
    paste(missing_sites, collapse = ", "),
    call. = FALSE
  )
}


# ----------------------------------------------------------------------------
# 5. Add fixed-year NUTS2 regions
# ----------------------------------------------------------------------------

dir.create(SPATIAL_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

if (file.exists(NUTS_CACHE_FILE)) {
  nuts2 <- readRDS(NUTS_CACHE_FILE)
} else {
  nuts2 <- giscoR::gisco_get_nuts(
    country = c("AT", "DE"),
    nuts_level = 2,
    year = "2024"
  )
  saveRDS(nuts2, NUTS_CACHE_FILE)
}

site_sf <- st_as_sf(
  site_spatial,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)
nuts2 <- st_transform(nuts2, st_crs(site_sf))

site_sf <- st_join(
  site_sf,
  nuts2[, c("NUTS_ID", "NUTS_NAME", "CNTR_CODE")],
  join = st_within,
  left = TRUE
)

missing_region <- which(is.na(site_sf$NUTS_ID))
if (length(missing_region)) {
  nearest <- st_nearest_feature(site_sf[missing_region, ], nuts2)
  site_sf$NUTS_ID[missing_region] <- nuts2$NUTS_ID[nearest]
  site_sf$NUTS_NAME[missing_region] <- nuts2$NUTS_NAME[nearest]
  site_sf$CNTR_CODE[missing_region] <- nuts2$CNTR_CODE[nearest]
}

site_attributes <- site_sf %>%
  st_drop_geometry() %>%
  transmute(
    Location,
    env_zone,
    region = as.character(NUTS_NAME),
    NUTS_ID = as.character(NUTS_ID)
  )

if (anyDuplicated(site_attributes$Location)) {
  stop(
    "Spatial enrichment produced more than one row for a physical site.",
    call. = FALSE
  )
}


# ----------------------------------------------------------------------------
# 6. Final checks and output
# ----------------------------------------------------------------------------

pheno_final <- pheno_combined %>%
  left_join(site_attributes, by = "Location", relationship = "many-to-one")

assert_columns(
  pheno_final,
  c(
    "BBCH", "Date", "SeasonYear", "Location", "env_id", "env_zone",
    "lat", "lon", "sowing_DOY", "Daylength",
    "cumDTU_beta", "swvl1", "swvl2", "swvl3", "swvl4",
    "soil_water_0_100_vwc", "FTSW_0_100_proxy"
  ),
  "Combined model input"
)

if (anyNA(pheno_final$region)) {
  stop("NUTS2 region is missing after spatial enrichment.", call. = FALSE)
}

qc_summary <- pheno_final %>%
  group_by(country) %>%
  summarise(
    n_observations = n(),
    n_sites = n_distinct(Location),
    n_environments = n_distinct(env_id),
    first_date = min(Date),
    last_date = max(Date),
    swvl1_missing_fraction = mean(is.na(swvl1)),
    swvl2_missing_fraction = mean(is.na(swvl2)),
    swvl3_missing_fraction = mean(is.na(swvl3)),
    swvl4_missing_fraction = mean(is.na(swvl4)),
    ftsw_0_100_missing_fraction = mean(is.na(FTSW_0_100_proxy)),
    .groups = "drop"
  )

write_csv(pheno_final, OUTPUT_FILE, na = "NA")
write_csv(qc_summary, QC_FILE, na = "NA")

message("Wrote combined model input: ", OUTPUT_FILE)
message("Wrote QC summary: ", QC_FILE)
print(qc_summary)

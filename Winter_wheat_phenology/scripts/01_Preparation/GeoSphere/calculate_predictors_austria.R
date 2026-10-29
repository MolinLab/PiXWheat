# =============================================================================
# Daily process-informed phenology and weather features for Austrian winter wheat
# =============================================================================
#
# Four response indices are calculated:
#   1. devfun                 = unphased response-function ablation
#   2. devfun_phase           = published Soltani/Sinclair phase logic
#   3. devfun_water           = unphased index with drought modifier
#   4. devfun_phase_water     = phase-aware index with drought modifier
#
# All cumulative and rolling predictors use information available on or before
# the current date. Full-season weather summaries are deliberately excluded.
# The process indices are fixed-parameter representations for ML benchmarking;
# they are not cultivar-calibrated predictions from a complete crop model.
# =============================================================================

library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(future)
library(future.apply)
library(data.table)
library(geosphere)
library(sf)

PROJECT_DIR <- Sys.getenv(
  "PHENOLOGY_PROJECT_DIR",
  unset = "."
)
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = TRUE)
WORK_DIR <- file.path(PROJECT_DIR, "Data", "GeoSphere")
if (!dir.exists(WORK_DIR)) stop("Austrian working directory not found: ", WORK_DIR)
setwd(WORK_DIR)


# =============================================================================
# 1. Configuration
# =============================================================================
PHENOLOGY_FILE      <- file.path(PROJECT_DIR, "Data", "PEP725", "pep_at_filtered_final_19_06.csv")
SOIL_HYDRAULIC_FILE <- "../Soil_hydraulic_ERA5_0_100.csv"
GDD_OUTPUT_DIR      <- "gdd"
FINAL_OUTPUT_FILE   <- "at_phenology_GDD_era5_test.csv"

SOIL_MAX_DISTANCE_M <- 200
SITE_KEY_DIGITS     <- 5
N_WORKERS           <- 5
RECENT_LOOKBACK_DAYS <- 30L

# SSM-Wheat beta temperature-response parameters
TBD <- 0
TOD <- 27.5
TCD <- 40

# Photoperiod response
CPP   <- 21
ppsen <- 0.004
SUN_ANGLE_BELOW_HORIZON <- 6

# Vernalization response
TBV   <- -1
TP1V  <- 0
TP2V  <- 8
TCV   <- 12
VDSAT <- 50
vsen  <- 0.030

# Fixed Austrian SSM-iCrop-type phase parameterization.
# The internally modelled first-node boundary is:
# 4.5 BD sowing-to-emergence +
# 1.2 BD emergence-to-tillering +
# 12.1 BD tillering-to-first-node = 17.8 BD.
# These are not universal SSM-Wheat constants and do not use observed BBCH.
BD_SOW_EMR     <- 4.5
BD_PHOTO_START <- 4.5
BD_PHOTO_END   <- 17.8
BD_VERN_START  <- 4.5
BD_VERN_END    <- 17.8
WSSG           <- 0.30
WSSD           <- 0.40

# Early de-vernalization parameters from Soltani and Sinclair (2012), Eq. 8.5.
DEVERN_CUMVER_LIMIT <- 10
DEVERN_TMAX_LIMIT   <- 30
DEVERN_RATE         <- 0.5

# Native ERA5-Land soil-layer thicknesses. Layer 4 is retained as a deep-soil
# state, but SoilGrids only reaches 200 cm and therefore cannot provide an exact
# field-capacity/wilting-point match for the full 100-289 cm ERA5 layer.
SWVL1_DEPTH_MM <- 70
SWVL2_DEPTH_MM <- 210
SWVL3_DEPTH_MM <- 720
SWVL4_DEPTH_MM <- 1890
SW_DEPTH_0_28_MM  <- SWVL1_DEPTH_MM + SWVL2_DEPTH_MM
SW_DEPTH_0_100_MM <- SW_DEPTH_0_28_MM + SWVL3_DEPTH_MM
SW_DEPTH_0_289_MM <- SW_DEPTH_0_100_MM + SWVL4_DEPTH_MM

# Temperature stress summaries
TLOW  <- 0
THIGH <- 30


# =============================================================================
# 2. Helpers
# =============================================================================

is_leap_year <- function(y) {
  ((y %% 4 == 0) & (y %% 100 != 0)) | (y %% 400 == 0)
}

days_in_year <- function(y) {
  if_else(is_leap_year(y), 366L, 365L)
}

make_site_key <- function(station_id, lat, lon, digits = SITE_KEY_DIGITS) {
  paste(
    station_id,
    sprintf(paste0("%.", digits, "f"), as.numeric(lat)),
    sprintf(paste0("%.", digits, "f"), as.numeric(lon)),
    sep = "__"
  )
}

# Create a stable physical extraction ID. Nearby experiments can share a
# source location_id, so the coordinate pair is part of the identifier;
# traj_id remains the experiment-season identifier.
make_location_id <- function(source_id, lat, lon, digits = SITE_KEY_DIGITS) {
  source_id <- as.character(source_id)
  coordinate_suffix <- paste0(
    "__", sprintf(paste0("%.", digits, "f"), as.numeric(lat)),
    "__", sprintf(paste0("%.", digits, "f"), as.numeric(lon))
  )
  if_else(
    endsWith(source_id, coordinate_suffix),
    source_id,
    paste0(source_id, coordinate_suffix)
  )
}

assert_unique_keys <- function(x, keys, object_name) {
  duplicates <- x %>%
    group_by(across(all_of(keys))) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n > 1L)
  
  if (nrow(duplicates) > 0L) {
    stop(
      object_name, " is not unique for keys ",
      paste(keys, collapse = ", "), ". Example duplicate keys:\n",
      paste(capture.output(print(head(duplicates, 10))), collapse = "\n"),
      call. = FALSE
    )
  }
  
  invisible(x)
}

# Return an all-NA vector if a variable is completely unavailable. Otherwise,
# skip isolated missing values and retain the corresponding observation count.
cumsum_available <- function(x) {
  if (all(is.na(x))) {
    rep(NA_real_, length(x))
  } else {
    cumsum(replace_na(x, 0))
  }
}

# Process states require a complete daily history. Once a required daily value
# is missing, the state remains NA rather than silently treating that day as 0.
cumsum_complete_to_date <- function(x) {
  running_sum <- cumsum(replace_na(x, 0))
  n_missing <- cumsum(is.na(x))
  if_else(n_missing == 0L, running_sum, NA_real_)
}

# Photoperiod including civil twilight, following the Kiesling calculation
# reproduced by Soltani and Sinclair (2012). geosphere::daylength() is retained
# separately. 
calc_civil_photoperiod <- function(
    lat,
    date,
    sun_angle_below_horizon = SUN_ANGLE_BELOW_HORIZON
) {
  doy <- yday(date)
  rad <- pi / 180
  solar_mean_anomaly <- 0.9856 * doy - 3.251
  solar_longitude <- solar_mean_anomaly +
    1.916 * sin(solar_mean_anomaly * rad) +
    0.020 * sin(2 * solar_mean_anomaly * rad) +
    282.565
  declination <- asin(0.39779 * sin(solar_longitude * rad))
  latitude <- as.numeric(lat) * rad
  solar_altitude <- -(sun_angle_below_horizon * rad)

  cos_hour_angle <- (
    sin(solar_altitude) - sin(latitude) * sin(declination)
  ) / (cos(latitude) * cos(declination))
  cos_hour_angle <- pmin(1, pmax(-1, cos_hour_angle))

  24 * acos(cos_hour_angle) / pi
}

# Eq. 8.4 with the early de-vernalization rule from Eq. 8.5 / Box 8.1.
# A missing required day invalidates the state from that day onward.
accumulate_vernalization <- function(
    verday,
    tmax,
    devernalization_limit = DEVERN_CUMVER_LIMIT,
    devernalization_tmax = DEVERN_TMAX_LIMIT,
    devernalization_rate = DEVERN_RATE
) {
  n <- length(verday)
  out <- rep(NA_real_, n)
  state <- 0
  complete_history <- TRUE

  for (i in seq_len(n)) {
    if (!complete_history || is.na(verday[i]) || is.na(tmax[i])) {
      complete_history <- FALSE
      next
    }

    state <- state + verday[i]
    if (
      state < devernalization_limit &&
      tmax[i] > devernalization_tmax
    ) {
      state <- state -
        devernalization_rate * (tmax[i] - devernalization_tmax)
    }
    state <- max(0, state)
    out[i] <- state
  }

  out
}

# Phase gating must depend on the model's own cumulative biological day, not
# observed BBCH. The current day's gate is evaluated from the previous day's
# state, matching the order in Soltani and Sinclair (2012), Box 8.1.
calc_soltani_phase_state <- function(
    tempfun,
    ppfun,
    verfun,
    emergence_bd = BD_SOW_EMR,
    photo_start_bd = BD_PHOTO_START,
    photo_end_bd = BD_PHOTO_END,
    vern_start_bd = BD_VERN_START,
    vern_end_bd = BD_VERN_END
) {
  n <- length(tempfun)
  ppfun_phase <- rep(NA_real_, n)
  verfun_phase <- rep(NA_real_, n)
  devfun_phase <- rep(NA_real_, n)
  photo_active <- rep(NA, n)
  vern_active <- rep(NA, n)
  post_emergence <- rep(NA, n)
  emerged_proxy <- rep(NA, n)

  cbd <- 0
  complete_history <- TRUE

  for (i in seq_len(n)) {
    if (
      !complete_history ||
      is.na(tempfun[i]) || is.na(ppfun[i]) || is.na(verfun[i])
    ) {
      complete_history <- FALSE
      next
    }

    photo_active[i] <- cbd >= photo_start_bd && cbd <= photo_end_bd
    vern_active[i] <- cbd >= vern_start_bd && cbd <= vern_end_bd
    post_emergence[i] <- cbd >= emergence_bd
    ppfun_phase[i] <- if (photo_active[i]) ppfun[i] else 1
    verfun_phase[i] <- if (vern_active[i]) verfun[i] else 1
    devfun_phase[i] <- tempfun[i] * ppfun_phase[i] * verfun_phase[i]
    cbd <- cbd + devfun_phase[i]
    emerged_proxy[i] <- cbd >= emergence_bd
  }

  tibble::tibble(
    ppfun_phase = ppfun_phase,
    verfun_phase = verfun_phase,
    soltani_photo_active = photo_active,
    soltani_vern_active = vern_active,
    soltani_response_active = photo_active | vern_active,
    post_emergence = post_emergence,
    emerged_proxy = emerged_proxy,
    devfun_phase = devfun_phase
  )
}

# The water-modified index needs its own cumulative biological-day state,
# because WSFD changes the day on which the later phase boundary is reached.
calc_soltani_phase_water_state <- function(
    tempfun,
    ppfun,
    verfun,
    wsfd,
    photo_start_bd = BD_PHOTO_START,
    photo_end_bd = BD_PHOTO_END,
    vern_start_bd = BD_VERN_START,
    vern_end_bd = BD_VERN_END
) {
  n <- length(tempfun)
  ppfun_phase_water <- rep(NA_real_, n)
  verfun_phase_water <- rep(NA_real_, n)
  devfun_phase_water <- rep(NA_real_, n)
  photo_active_water <- rep(NA, n)
  vern_active_water <- rep(NA, n)

  cbd_water <- 0
  complete_history <- TRUE

  for (i in seq_len(n)) {
    if (
      !complete_history ||
      anyNA(c(tempfun[i], ppfun[i], verfun[i], wsfd[i]))
    ) {
      complete_history <- FALSE
      next
    }

    photo_active_water[i] <-
      cbd_water >= photo_start_bd && cbd_water <= photo_end_bd
    vern_active_water[i] <-
      cbd_water >= vern_start_bd && cbd_water <= vern_end_bd

    ppfun_phase_water[i] <- if (photo_active_water[i]) ppfun[i] else 1
    verfun_phase_water[i] <- if (vern_active_water[i]) verfun[i] else 1
    devfun_phase_water[i] <-
      tempfun[i] *
      ppfun_phase_water[i] *
      verfun_phase_water[i] *
      wsfd[i]

    cbd_water <- cbd_water + devfun_phase_water[i]
  }

  tibble::tibble(
    ppfun_phase_water = ppfun_phase_water,
    verfun_phase_water = verfun_phase_water,
    soltani_photo_active_water = photo_active_water,
    soltani_vern_active_water = vern_active_water,
    soltani_response_active_water =
      photo_active_water | vern_active_water,
    devfun_phase_water = devfun_phase_water
  )
}

rolling_sum_available <- function(x, width) {
  x_sum <- cumsum(replace_na(x, 0))
  x_sum - lag(x_sum, width, default = 0)
}

rolling_n_available <- function(x, width) {
  x_n <- cumsum(!is.na(x))
  x_n - lag(x_n, width, default = 0)
}

rolling_calendar_complete <- function(date, width) {
  date_integer <- as.integer(date)
  seq_along(date_integer) >= width &
    date_integer - lag(date_integer, width - 1L) == width - 1L
}

rolling_sum_complete <- function(x, date, width) {
  if_else(
    rolling_calendar_complete(date, width) &
      rolling_n_available(x, width) == width,
    rolling_sum_available(x, width),
    NA_real_
  )
}

rolling_mean_complete <- function(x, date, width) {
  if_else(
    rolling_calendar_complete(date, width) &
      rolling_n_available(x, width) == width,
    rolling_sum_available(x, width) / width,
    NA_real_
  )
}

read_weather_year <- function(var, year, folder) {
  file <- file.path(folder, paste0(var, "_", year, ".csv"))
  if (!file.exists(file)) return(NULL)
  
  fread(file, showProgress = FALSE) %>%
    mutate(Date = as.Date(Date))
}

read_weather_pair <- function(var, harvest_year, folder) {
  previous <- read_weather_year(var, harvest_year - 1L, folder)
  current  <- read_weather_year(var, harvest_year, folder)
  
  if (is.null(previous) || is.null(current)) return(NULL)
  bind_rows(previous, current)
}

normalize_weather_keys <- function(x) {
  x %>%
    mutate(
      location_id = make_location_id(location_id, lat_orig, lon_orig),
      year = as.integer(year),
      DOY = as.integer(DOY),
      Date = as.Date(Date)
    )
}


# =============================================================================
# 3. Phenology sites and sowing/harvest windows
# =============================================================================

phen_df <- read_csv(PHENOLOGY_FILE, show_col_types = FALSE) %>%
  filter(coalesce(as.character(traj_class), "") != "too_sparse") %>%
  select(
    any_of("traj_id"),
    location_id,
    location,
    Location,
    Location_name,
    SeasonYear,
    season_date,
    lat,
    lon,
    BBCH_adjusted
    ) %>%
  transmute(
    across(any_of("traj_id")),
    location_id = make_location_id(location_id, lat, lon),
    location=coalesce(Location,location,Location_name),
    season_year = as.integer(SeasonYear),
    year = as.integer(SeasonYear),
    phen_date = as.Date(season_date),
    DOY = yday(phen_date),
    lat_orig = as.numeric(lat),
    lon_orig = as.numeric(lon),
    BBCH = BBCH_adjusted
  ) %>%
  filter(
    !is.na(traj_id),
    !is.na(location_id),
    !is.na(phen_date),
    !is.na(lat_orig),
    !is.na(lon_orig)
  ) %>%
  mutate(site_key = make_site_key(traj_id, lat_orig, lon_orig))

phen_sites <- phen_df %>%
  distinct(site_key, traj_id, lat_orig, lon_orig)

assert_unique_keys(phen_sites, "site_key", "phen_sites")

station_coordinate_versions <- phen_sites %>%
  count(traj_id, name = "n_coordinate_versions") %>%
  filter(n_coordinate_versions > 1L)

if (nrow(station_coordinate_versions) > 0L) {
  message(
    nrow(station_coordinate_versions),
    " station IDs have multiple coordinate versions; site_key keeps them separate."
  )
}

sow_harv <- phen_df %>%
  group_by(
    site_key, traj_id, location_id, location, year, lat_orig, lon_orig
  ) %>%
  summarise(
    sow_doy = {
      observed_sowing <- DOY[BBCH == 0 & !is.na(DOY)]
      if (length(observed_sowing) > 0L) {
        min(observed_sowing)
      } else {
        NA_real_
      }
    },
    harv_doy = {
      bbch89 <- DOY[BBCH == 89 & !is.na(DOY)]
      bbch87 <- DOY[BBCH == 87 & !is.na(DOY)]
      
      if (length(bbch89) > 0L) {
        max(bbch89)
      } else if (length(bbch87) > 0L) {
        max(bbch87) + 14
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  ) %>%
  mutate(
    sow_year = year - 1L,
    harv_year = year,
    harv_missing = is.na(harv_doy),
    harv_doy_use = if_else(
      harv_missing,
      as.numeric(days_in_year(harv_year)),
      pmin(harv_doy, as.numeric(days_in_year(harv_year)))
    ),
    sow_year_days = as.numeric(days_in_year(sow_year))
  ) %>%
  filter(!is.na(sow_doy))

assert_unique_keys(
  sow_harv,
  c("site_key", "harv_year"),
  "sow_harv"
)


# =============================================================================
# 4. Static SoilGrids hydraulic limits: one nearest match per site
# =============================================================================

if (!file.exists(SOIL_HYDRAULIC_FILE)) {
  stop("Soil hydraulic file not found: ", SOIL_HYDRAULIC_FILE, call. = FALSE)
}

required_hydraulic_columns <- c(
  "lat", "lon",
  "theta_fc10_l1", "theta_fc10_l2", "theta_fc10_l3",
  "theta_fc33_l1", "theta_fc33_l2", "theta_fc33_l3",
  "theta_wp_l1", "theta_wp_l2", "theta_wp_l3",
  "theta_fc10_0_28", "theta_fc33_0_28", "theta_wp_0_28",
  "theta_fc10_0_100", "theta_fc33_0_100", "theta_wp_0_100",
  "paw_fc10_0_28_mm", "paw_fc33_0_28_mm",
  "paw_fc10_0_100_mm", "paw_fc33_0_100_mm"
)

soil_df <- read_csv(SOIL_HYDRAULIC_FILE, show_col_types = FALSE)

missing_hydraulic_columns <- setdiff(
  required_hydraulic_columns,
  names(soil_df)
)

if (length(missing_hydraulic_columns) > 0L) {
  stop(
    "The soil hydraulic file is missing: ",
    paste(missing_hydraulic_columns, collapse = ", "),
    call. = FALSE
  )
}

soil_df <- soil_df %>%
  mutate(
    lat = as.numeric(lat),
    lon = as.numeric(lon)
  ) %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  distinct(lat, lon, .keep_all = TRUE) %>%
  mutate(soil_id = row_number())

# Fail early if the SoilGrids water-retention values still have a unit error.
median_fc33 <- median(soil_df$theta_fc33_0_100, na.rm = TRUE)
if (!is.finite(median_fc33) || median_fc33 < 0.05 || median_fc33 > 0.80) {
  stop(
    "Implausible theta_fc33_0_100 median (", signif(median_fc33, 4),
    "). Hydraulic limits must be volumetric fractions in m3 m-3 before FTSW is calculated.",
    call. = FALSE
  )
}

phen_sites_sf <- st_as_sf(
  phen_sites,
  coords = c("lon_orig", "lat_orig"),
  crs = 4326,
  remove = FALSE
)

soil_sf <- st_as_sf(
  soil_df,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

nearest_idx <- st_nearest_feature(phen_sites_sf, soil_sf)

site_soil_match_all <- phen_sites %>%
  mutate(
    soil_id = soil_df$soil_id[nearest_idx],
    soil_dist_m = as.numeric(
      st_distance(
        phen_sites_sf,
        soil_sf[nearest_idx, ],
        by_element = TRUE
      )
    )
  ) %>%
  left_join(
    soil_df %>% rename(soil_lat = lat, soil_lon = lon),
    by = "soil_id",
    relationship = "many-to-one"
  )

site_soil_match <- site_soil_match_all %>%
  filter(soil_dist_m <= SOIL_MAX_DISTANCE_M)

assert_unique_keys(site_soil_match, "site_key", "site_soil_match")

soil_match_summary <- site_soil_match_all %>%
  summarise(
    n_sites = n(),
    n_within_threshold = sum(soil_dist_m <= SOIL_MAX_DISTANCE_M),
    match_rate = mean(soil_dist_m <= SOIL_MAX_DISTANCE_M),
    median_distance_m = median(soil_dist_m),
    max_distance_m = max(soil_dist_m)
  )

print(soil_match_summary)


# =============================================================================
# 5. Calculate daily features for one harvest year
# =============================================================================

calc_gdd_winter_geo <- function(y) {
  
  message("Processing harvest year: ", y)
  
  season_table <- sow_harv %>%
    filter(harv_year == y) %>%
    select(
      site_key, traj_id, location_id, location, lat_orig, lon_orig,
      sow_year, harv_year, sow_doy, harv_doy_use, harv_missing
    )
  
  if (nrow(season_table) == 0L) return(invisible(NULL))
  
  assert_unique_keys(
    season_table,
    c("site_key", "harv_year"),
    paste0("season_table for ", y)
  )
  
  weather_raw <- list(
    tmax  = read_weather_pair("tmax",  y, "TX/temp_extract"),
    tmin  = read_weather_pair("tmin",  y, "TN/temp_extract"),
    prec  = read_weather_pair("prec",  y, "RR/temp_extract"),
    ssrd  = read_weather_pair("ssrd",  y, "../Era5/era5land_yearly/temp_extract_geosphere/"),
    e     = read_weather_pair("e",     y, "../Era5/era5land_yearly/temp_extract_geosphere/"),
    pev   = read_weather_pair("pev",   y, "../Era5/era5land_yearly/temp_extract_geosphere/"),
    swvl1 = read_weather_pair("swvl1", y, "../Era5/era5land_yearly/temp_extract_geosphere/"),
    swvl2 = read_weather_pair("swvl2", y, "../Era5/era5land_yearly/temp_extract_geosphere/"),
    swvl3 = read_weather_pair("swvl3", y, "../Era5/era5land_yearly/temp_extract_geosphere/"),
    swvl4 = read_weather_pair("swvl4", y, "../Era5/era5land_yearly/temp_extract_geosphere/")
  )
  
  missing_weather <- names(weather_raw)[vapply(weather_raw, is.null, logical(1))]
  if (length(missing_weather) > 0L) {
    warning(
      "Skipping harvest year ", y, ". Missing previous/current files for: ",
      paste(missing_weather, collapse = ", ")
    )
    return(invisible(NULL))
  }
  
  daily_keys <- c("location_id", "year", "DOY", "Date")
  
  tmax_df <- normalize_weather_keys(weather_raw$tmax) %>%
    select(
      location_id, location, lat_orig, lon_orig,
      year, DOY, Date, tmax
    )
  
  tmin_df <- normalize_weather_keys(weather_raw$tmin) %>%
    select(all_of(daily_keys), tmin)
  
  prec_df <- normalize_weather_keys(weather_raw$prec) %>%
    select(all_of(daily_keys), prec)
  
  rad_df <- normalize_weather_keys(weather_raw$ssrd) %>%
    select(all_of(daily_keys), ssrd) %>%
    mutate(rad = ssrd / 1e6) %>% distinct()
  
  e_df <- normalize_weather_keys(weather_raw$e) %>%
    select(all_of(daily_keys), any_of(c("e", "evaporation_mm"))) %>% distinct()
  if (!"e" %in% names(e_df)) e_df$e <- NA_real_
  if (!"evaporation_mm" %in% names(e_df)) {
    e_df$evaporation_mm <- -1000 * e_df$e
  } 
  
  pev_df <- normalize_weather_keys(weather_raw$pev) %>%
    select(all_of(daily_keys), any_of(c("pev", "potential_evaporation_mm"))) %>% distinct()
  if (!"pev" %in% names(pev_df)) pev_df$pev <- NA_real_
  if (!"potential_evaporation_mm" %in% names(pev_df)) {
    pev_df$potential_evaporation_mm <- -1000 * pev_df$pev
  }
  
  swvl1_df <- normalize_weather_keys(weather_raw$swvl1) %>%
    select(all_of(daily_keys), swvl1) %>% distinct()
  
  swvl2_df <- normalize_weather_keys(weather_raw$swvl2) %>%
    select(all_of(daily_keys), swvl2) %>% distinct()

  swvl3_df <- normalize_weather_keys(weather_raw$swvl3) %>%
    select(all_of(daily_keys), swvl3) %>% distinct()

  swvl4_df <- normalize_weather_keys(weather_raw$swvl4) %>%
    select(all_of(daily_keys), swvl4) %>% distinct()
  
  weather_tables <- list(
    tmax = tmax_df,
    tmin = tmin_df,
    prec = prec_df,
    rad = rad_df,
    e = e_df,
    pev = pev_df,
    swvl1 = swvl1_df,
    swvl2 = swvl2_df,
    swvl3 = swvl3_df,
    swvl4 = swvl4_df
  )
  
  for (table_name in names(weather_tables)) {
    assert_unique_keys(
      weather_tables[[table_name]],
      daily_keys,
      paste0(table_name, " weather for ", y)
    )
  }
  
  df <- tmax_df %>%
    left_join(tmin_df,  by = daily_keys, relationship = "one-to-one") %>%
    left_join(prec_df,  by = daily_keys, relationship = "one-to-one") %>%
    left_join(rad_df,   by = daily_keys, relationship = "one-to-one") %>%
    left_join(e_df,     by = daily_keys, relationship = "one-to-one") %>%
    left_join(pev_df,   by = daily_keys, relationship = "one-to-one") %>%
    left_join(swvl1_df, by = daily_keys, relationship = "one-to-one") %>%
    left_join(swvl2_df, by = daily_keys, relationship = "one-to-one") %>%
    left_join(swvl3_df, by = daily_keys, relationship = "one-to-one") %>%
    left_join(swvl4_df, by = daily_keys, relationship = "one-to-one")
  
  assert_unique_keys(df, daily_keys, paste0("weather-joined data for ", y))
  
  # Restrict the daily trajectory to sowing--harvest/year-end. Observed harvest
  # only defines the extraction window; it is never used as a predictor.
  df <- df %>%
    inner_join(
      season_table %>%
        select(
          site_key, traj_id, location_id, sow_year, harv_year,
          sow_doy, harv_doy_use, harv_missing
        ),
      by = "location_id",
      relationship = "many-to-many"
    ) %>%
    mutate(
      season_year = harv_year,
      sowing_date = as.Date(paste0(sow_year, "-01-01")) +
        as.integer(sow_doy) - 1L,
      season_end_date = as.Date(paste0(harv_year, "-01-01")) +
        as.integer(harv_doy_use) - 1L
    ) %>%
    filter(
      Date >= sowing_date - (RECENT_LOOKBACK_DAYS - 1L),
      Date <= season_end_date
    ) %>%
    mutate(
      days_since_sowing = as.integer(Date - sowing_date),
      tmean = (tmax + tmin) / 2,
      evaporation_mm = coalesce(evaporation_mm, -1000 * e),
      potential_evaporation_mm = coalesce(
        potential_evaporation_mm,
        -1000 * pev
      ),
      soil_water_0_28_vwc = (
        SWVL1_DEPTH_MM * swvl1 + SWVL2_DEPTH_MM * swvl2
      ) / SW_DEPTH_0_28_MM,
      soil_water_0_28_mm =
        SWVL1_DEPTH_MM * swvl1 + SWVL2_DEPTH_MM * swvl2,
      soil_water_0_100_vwc = (
        SWVL1_DEPTH_MM * swvl1 +
          SWVL2_DEPTH_MM * swvl2 +
          SWVL3_DEPTH_MM * swvl3
      ) / SW_DEPTH_0_100_MM,
      soil_water_0_100_mm =
        SWVL1_DEPTH_MM * swvl1 +
        SWVL2_DEPTH_MM * swvl2 +
        SWVL3_DEPTH_MM * swvl3,
      soil_water_0_289_vwc = (
        SWVL1_DEPTH_MM * swvl1 +
          SWVL2_DEPTH_MM * swvl2 +
          SWVL3_DEPTH_MM * swvl3 +
          SWVL4_DEPTH_MM * swvl4
      ) / SW_DEPTH_0_289_MM,
      soil_water_0_289_mm =
        SWVL1_DEPTH_MM * swvl1 +
        SWVL2_DEPTH_MM * swvl2 +
        SWVL3_DEPTH_MM * swvl3 +
        SWVL4_DEPTH_MM * swvl4,
      swvl1_minus_swvl2 = swvl1 - swvl2,
      swvl2_minus_swvl3 = swvl2 - swvl3,
      swvl3_minus_swvl4 = swvl3 - swvl4,
      et_gap_mm = potential_evaporation_mm - evaporation_mm,
      et_fraction = if_else(
        !is.na(potential_evaporation_mm) & potential_evaporation_mm > 0.05,
        evaporation_mm / potential_evaporation_mm,
        NA_real_
      ),
      cold_day = if_else(is.na(tmin), NA_integer_, as.integer(tmin < TLOW)),
      hot_day = if_else(is.na(tmax), NA_integer_, as.integer(tmax > THIGH)),
      CDD = case_when(
        is.na(tmin) ~ NA_real_,
        tmin <= TLOW ~ TLOW - tmin,
        TRUE ~ 0
      ),
      HDD = case_when(
        is.na(tmax) ~ NA_real_,
        tmax >= THIGH ~ tmax - THIGH,
        TRUE ~ 0
      )
    )
  
  assert_unique_keys(
    df,
    c("site_key", "season_year", "Date"),
    paste0("season-window data for ", y)
  )
  
  # One static soil record per station-coordinate site. This join cannot
  # multiply daily rows because site_soil_match is unique by site_key.
  df <- df %>%
    left_join(
      site_soil_match %>%
        select(-traj_id,-lat_orig, -lon_orig),
      by = "site_key",
      relationship = "many-to-one"
    )
  
  assert_unique_keys(
    df,
    c("site_key", "season_year", "Date"),
    paste0("soil-joined data for ", y)
  )
  
  # Optional QC columns exist in the current hydraulic file. These fallbacks
  # keep the calculation explicit if an older compatible file is supplied.
  if (!"hydraulic_qc_p33_0_28" %in% names(df)) {
    df$hydraulic_qc_p33_0_28 <- TRUE
  }
  if (!"hydraulic_qc_p10_0_28" %in% names(df)) {
    df$hydraulic_qc_p10_0_28 <- TRUE
  }
  if (!"hydraulic_qc_p33_0_100" %in% names(df)) {
    df$hydraulic_qc_p33_0_100 <- TRUE
  }
  if (!"hydraulic_qc_p10_0_100" %in% names(df)) {
    df$hydraulic_qc_p10_0_100 <- TRUE
  }
  if (!"hydraulic_high_uncertainty" %in% names(df)) {
    df$hydraulic_high_uncertainty <- NA
  }
  if (!"hydraulic_max_cv_0_100" %in% names(df)) {
    df$hydraulic_max_cv_0_100 <- NA_real_
  }
  if (!"standort" %in% names(df)) df$standort <- NA_character_

  swvl_columns <- c("swvl1", "swvl2", "swvl3", "swvl4")
  missing_swvl <- setdiff(swvl_columns, names(df))
  if (length(missing_swvl) > 0L) {
    stop(
      "Missing ERA5-Land soil-water columns: ",
      paste(missing_swvl, collapse = ", "),
      call. = FALSE
    )
  }
  swvl_values <- unlist(
    lapply(swvl_columns, function(column) df[[column]]),
    use.names = FALSE
  )
  if (any(swvl_values < 0 | swvl_values > 1, na.rm = TRUE)) {
    stop(
      "Harvest year ", y,
      " contains ERA5-Land soil-water values outside [0, 1] m3 m-3.",
      call. = FALSE
    )
  }

  # Fixed retrospective windows include pre-sowing weather. The 30-day
  # lookback rows are used only here and are removed before process states and
  # cumulative-from-sowing features are calculated.
  df <- df %>%
    arrange(site_key, season_year, Date) %>%
    group_by(site_key, traj_id, lat_orig, lon_orig, season_year) %>%
    mutate(
      window_complete_14 = rolling_calendar_complete(Date, 14),
      window_complete_30 = rolling_calendar_complete(Date, 30),
      sumRR_14 = rolling_sum_complete(prec, Date, 14),
      sumRR_30 = rolling_sum_complete(prec, Date, 30),
      sumRad_14 = rolling_sum_complete(rad, Date, 14),
      sumRad_30 = rolling_sum_complete(rad, Date, 30),
      sumCDD_14 = rolling_sum_complete(CDD, Date, 14),
      sumHDD_14 = rolling_sum_complete(HDD, Date, 14),
      sumCDD_30 = rolling_sum_complete(CDD, Date, 30),
      sumHDD_30 = rolling_sum_complete(HDD, Date, 30),
      sumEvap_14 = rolling_sum_complete(evaporation_mm, Date, 14),
      sumPEvap_14 = rolling_sum_complete(
        potential_evaporation_mm, Date, 14
      ),
      sumEvap_30 = rolling_sum_complete(evaporation_mm, Date, 30),
      sumPEvap_30 = rolling_sum_complete(
        potential_evaporation_mm, Date, 30
      ),
      sumEtGap_14 = rolling_sum_complete(et_gap_mm, Date, 14),
      meanSWVL1_14 = rolling_mean_complete(swvl1, Date, 14),
      meanSWVL2_14 = rolling_mean_complete(swvl2, Date, 14),
      meanSWVL3_14 = rolling_mean_complete(swvl3, Date, 14),
      meanSWVL4_14 = rolling_mean_complete(swvl4, Date, 14),
      meanSoilWater_0_28_14 = rolling_mean_complete(
        soil_water_0_28_vwc, Date, 14
      ),
      meanSoilWater_0_100_14 = rolling_mean_complete(
        soil_water_0_100_vwc, Date, 14
      ),
      nSoilWater_14 = rolling_n_available(soil_water_0_28_vwc, 14),
      nSoilWater_0_100_14 = rolling_n_available(soil_water_0_100_vwc, 14),
      nEvap_14 = rolling_n_available(evaporation_mm, 14),
      nPEvap_14 = rolling_n_available(potential_evaporation_mm, 14),
      etFraction_14 = if_else(
        sumPEvap_14 > 0.05,
        sumEvap_14 / sumPEvap_14,
        NA_real_
      ),
      climaticWaterBalance_14 = sumRR_14 - sumPEvap_14,
      sumEvap_30 = rolling_sum_complete(evaporation_mm, Date, 30),
      sumPEvap_30 = rolling_sum_complete(
        potential_evaporation_mm, Date, 30
      ),
      sumEtGap_30 = rolling_sum_complete(et_gap_mm, Date, 30),
      meanSWVL1_30 = rolling_mean_complete(swvl1, Date, 30),
      meanSWVL2_30 = rolling_mean_complete(swvl2, Date, 30),
      meanSWVL3_30 = rolling_mean_complete(swvl3, Date, 30),
      meanSWVL4_30 = rolling_mean_complete(swvl4, Date, 30),
      meanSoilWater_0_28_30 = rolling_mean_complete(
        soil_water_0_28_vwc, Date, 30
      ),
      meanSoilWater_0_100_30 = rolling_mean_complete(
        soil_water_0_100_vwc, Date, 30
      ),
      nSoilWater_30 = rolling_n_available(soil_water_0_28_vwc, 30),
      nSoilWater_0_100_30 = rolling_n_available(soil_water_0_100_vwc, 30),
      nEvap_30 = rolling_n_available(evaporation_mm, 30),
      nPEvap_30 = rolling_n_available(potential_evaporation_mm, 30),
      etFraction_30 = if_else(
        sumPEvap_30 > 0.05,
        sumEvap_30 / sumPEvap_30,
        NA_real_
      ),
      climaticWaterBalance_30 = sumRR_30 - sumPEvap_30,
      soilWaterChange_7 = soil_water_0_28_vwc - lag(soil_water_0_28_vwc, 7),
      soilWaterChange_14 =
        soil_water_0_28_vwc - lag(soil_water_0_28_vwc, 14),
      soilWaterChange_0_100_30 =
        soil_water_0_100_vwc - lag(soil_water_0_100_vwc, 30),
      soilWaterChange_0_100_14 =
        soil_water_0_100_vwc - lag(soil_water_0_100_vwc, 14)
    ) %>%
    filter(days_since_sowing >= 0L) %>%
    ungroup()

  missing_tmean <- sum(is.na(df$tmean))
  if (missing_tmean > 0L) {
    warning(
      "Harvest year ", y, " contains ", missing_tmean,
      " daily rows with missing mean temperature. Process states become NA from the first missing day."
    )
  }
  
  # ---------------------------------------------------------------------------
  # 5a. Daily biological response functions
  # ---------------------------------------------------------------------------
  df <- df %>%
    arrange(site_key, season_year, Date) %>%
    group_by(site_key,traj_id, lat_orig, lon_orig, season_year) %>%
    mutate(
      daylength = geosphere::daylength(lat_orig, yday(Date)),
      photoperiod_civil = calc_civil_photoperiod(lat_orig, Date),
      TMP = tmean,
      
      #tempfun = case_when(
      #  is.na(TMP) ~ NA_real_,
      #  TMP <= TBD ~ 0,
      #  TMP < TP1D ~ (TMP - TBD) / (TP1D - TBD),
      #  TMP <= TP2D ~ 1,
      #  TMP < TCD ~ (TCD - TMP) / (TCD - TP2D),
      #  TRUE ~ 0
      #),
      
      tempfun = case_when(
        is.na(TMP) ~ NA_real_,
        TMP <= TBD | TMP >= TCD ~ 0,
        TRUE ~
          ((TCD - TMP) / (TCD - TOD)) *
          ((TMP - TBD) / (TOD - TBD)) ^
          ((TOD - TBD) / (TCD - TOD))
      ),
      
      GDD = if_else(is.na(TMP), NA_real_, pmax(TMP - TBD, 0)),
      DTU_beta  = (TOD - TBD) * tempfun,
      
      ppfun_raw = case_when(
        is.na(photoperiod_civil) ~ NA_real_,
        photoperiod_civil < CPP ~
          1 - ppsen * (CPP - photoperiod_civil)^2,
        TRUE ~ 1
      ),
      ppfun = pmin(1, pmax(0, ppfun_raw)),
      
      # Mean air temperature is an explicit proxy for crown temperature because
      # snow depth/crown temperature is unavailable.
      VERDAY = case_when(
        is.na(TMP) ~ NA_real_,
        TMP <= TBV ~ 0,
        TMP < TP1V ~ (TMP - TBV) / (TP1V - TBV),
        TMP <= TP2V ~ 1,
        TMP < TCV ~ (TCV - TMP) / (TCV - TP2V),
        TRUE ~ 0
      )
    ) %>%
    
    # -------------------------------------------------------------------------
  # 5b. Vernalization and published biological-day phase logic
  # -------------------------------------------------------------------------
  mutate(
    cumGDD = cumsum_complete_to_date(GDD),
    cumDTU_beta = cumsum_complete_to_date(DTU_beta),
    cumVER = accumulate_vernalization(VERDAY, tmax),
    
    verfun = case_when(
      is.na(cumVER) ~ NA_real_,
      TRUE ~ pmin(1, pmax(0, 1 - vsen * (VDSAT - cumVER)))
    ),
    vern_saturated = coalesce(cumVER >= VDSAT, FALSE),
    # Unphased product is retained only as an explicit ablation.
    devfun = tempfun * ppfun * verfun
  ) %>%
  group_modify(~ {
    phase_state <- calc_soltani_phase_state(
      tempfun = .x$tempfun,
      ppfun = .x$ppfun,
      verfun = .x$verfun
    )
    bind_cols(.x, phase_state)
  }) %>%
    
    # -------------------------------------------------------------------------
  # 5c. Layer-matched FTSW and water-modified development
  # -------------------------------------------------------------------------
  mutate(
    theta_fc_0_28 = theta_fc33_0_28,
    theta_fc_0_100 = theta_fc33_0_100,
    hydraulic_source_0_28 = if_else(
      coalesce(hydraulic_qc_p33_0_28, FALSE) &
        !is.na(theta_fc33_l1) & !is.na(theta_fc33_l2),
      "SoilGrids_p33_layer_matched_0_28",
      NA_character_
    ),
    hydraulic_source_0_100 = if_else(
      coalesce(hydraulic_qc_p33_0_100, FALSE) &
        !is.na(theta_fc33_l1) & !is.na(theta_fc33_l2) &
        !is.na(theta_fc33_l3),
      "SoilGrids_p33_layer_matched_0_100",
      NA_character_
    ),
    hydraulic_limits_0_28_available =
      coalesce(hydraulic_qc_p33_0_28, FALSE) &
      !is.na(swvl1) & !is.na(swvl2) &
      !is.na(theta_fc33_l1) & !is.na(theta_fc33_l2) &
      !is.na(theta_wp_l1) & !is.na(theta_wp_l2) &
      theta_fc33_l1 > theta_wp_l1 &
      theta_fc33_l2 > theta_wp_l2,
    hydraulic_limits_p10_0_28_available =
      coalesce(hydraulic_qc_p10_0_28, FALSE) &
      !is.na(swvl1) & !is.na(swvl2) &
      !is.na(theta_fc10_l1) & !is.na(theta_fc10_l2) &
      !is.na(theta_wp_l1) & !is.na(theta_wp_l2) &
      theta_fc10_l1 > theta_wp_l1 &
      theta_fc10_l2 > theta_wp_l2,
    hydraulic_limits_0_100_available =
      coalesce(hydraulic_qc_p33_0_100, FALSE) &
      !is.na(swvl1) & !is.na(swvl2) & !is.na(swvl3) &
      !is.na(theta_fc33_l1) & !is.na(theta_fc33_l2) &
      !is.na(theta_fc33_l3) & !is.na(theta_wp_l1) &
      !is.na(theta_wp_l2) & !is.na(theta_wp_l3) &
      theta_fc33_l1 > theta_wp_l1 &
      theta_fc33_l2 > theta_wp_l2 &
      theta_fc33_l3 > theta_wp_l3,
    hydraulic_limits_p10_0_100_available =
      coalesce(hydraulic_qc_p10_0_100, FALSE) &
      !is.na(swvl1) & !is.na(swvl2) & !is.na(swvl3) &
      !is.na(theta_fc10_l1) & !is.na(theta_fc10_l2) &
      !is.na(theta_fc10_l3) & !is.na(theta_wp_l1) &
      !is.na(theta_wp_l2) & !is.na(theta_wp_l3) &
      theta_fc10_l1 > theta_wp_l1 &
      theta_fc10_l2 > theta_wp_l2 &
      theta_fc10_l3 > theta_wp_l3,
    available_soil_water_0_28_mm =
      SWVL1_DEPTH_MM * pmax(swvl1 - theta_wp_l1,0) +
      SWVL2_DEPTH_MM * pmax(swvl2 - theta_wp_l2,0),
    total_available_water_p33_0_28_mm =
      SWVL1_DEPTH_MM * (theta_fc33_l1 - theta_wp_l1) +
      SWVL2_DEPTH_MM * (theta_fc33_l2 - theta_wp_l2),
    total_available_water_p10_0_28_mm =
      SWVL1_DEPTH_MM * (theta_fc10_l1 - theta_wp_l1) +
      SWVL2_DEPTH_MM * (theta_fc10_l2 - theta_wp_l2),
    
    FTSW_0_28_raw = if_else(
      hydraulic_limits_0_28_available &
        total_available_water_p33_0_28_mm > 0,
      available_soil_water_0_28_mm / total_available_water_p33_0_28_mm,
      NA_real_
    ),
    FTSW_0_28_proxy = pmin(1, pmax(0, FTSW_0_28_raw)),
    FTSW_0_28_p10_raw = if_else(
      hydraulic_limits_p10_0_28_available &
        total_available_water_p10_0_28_mm > 0,
      available_soil_water_0_28_mm / total_available_water_p10_0_28_mm,
      NA_real_
    ),
    FTSW_0_28_p10_proxy = pmin(1, pmax(0, FTSW_0_28_p10_raw)),

    available_soil_water_0_100_mm =
      available_soil_water_0_28_mm +
      SWVL3_DEPTH_MM * pmax(swvl3 - theta_wp_l3,0),
    total_available_water_p33_0_100_mm =
      SWVL1_DEPTH_MM * (theta_fc33_l1 - theta_wp_l1) +
      SWVL2_DEPTH_MM * (theta_fc33_l2 - theta_wp_l2) +
      SWVL3_DEPTH_MM * (theta_fc33_l3 - theta_wp_l3),
    total_available_water_p10_0_100_mm =
      SWVL1_DEPTH_MM * (theta_fc10_l1 - theta_wp_l1) +
      SWVL2_DEPTH_MM * (theta_fc10_l2 - theta_wp_l2) +
      SWVL3_DEPTH_MM * (theta_fc10_l3 - theta_wp_l3),
    FTSW_0_100_raw = if_else(
      hydraulic_limits_0_100_available &
        total_available_water_p33_0_100_mm > 0,
      available_soil_water_0_100_mm /
        total_available_water_p33_0_100_mm,
      NA_real_
    ),
    FTSW_0_100_proxy = pmin(1, pmax(0, FTSW_0_100_raw)),
    FTSW_0_100_p10_raw = if_else(
      hydraulic_limits_p10_0_100_available &
        total_available_water_p10_0_100_mm > 0,
      available_soil_water_0_100_mm /
        total_available_water_p10_0_100_mm,
      NA_real_
    ),
    FTSW_0_100_p10_proxy = pmin(1, pmax(0, FTSW_0_100_p10_raw)),
    
    # Eq. 15.3: growth is unaffected above WSSG and decreases linearly
    # below WSSG. Both depths are retained as explicitly labelled proxies.
    WSFG_0_28_proxy = case_when(
      !hydraulic_limits_0_28_available ~ NA_real_,
      is.na(FTSW_0_28_proxy) ~ NA_real_,
      FTSW_0_28_proxy >= WSSG ~ 1,
      TRUE ~ FTSW_0_28_proxy / WSSG
    ),
    WSFG_0_100_proxy = case_when(
      !hydraulic_limits_0_100_available ~ NA_real_,
      is.na(FTSW_0_100_proxy) ~ NA_real_,
      FTSW_0_100_proxy >= WSSG ~ 1,
      TRUE ~ FTSW_0_100_proxy / WSSG
    ),

    # Eq. 15.8: WSFD = (1 - WSFG) * WSSD + 1. Following Eq. 15.9,
    # the development modifier is neutral until the internally simulated
    # biological-day state has passed emergence.
    WSFD_0_28_proxy = case_when(
      is.na(post_emergence) ~ NA_real_,
      !post_emergence ~ 1,
      is.na(WSFG_0_28_proxy) ~ NA_real_,
      TRUE ~ (1 - WSFG_0_28_proxy) * WSSD + 1
    ),
    WSFD_0_100_proxy = case_when(
      is.na(post_emergence) ~ NA_real_,
      !post_emergence ~ 1,
      is.na(WSFG_0_100_proxy) ~ NA_real_,
      TRUE ~ (1 - WSFG_0_100_proxy) * WSSD + 1
    ),

    # Preserve the established process definition: the development index
    # uses the shallow 0-28 cm proxy. The 0-100 cm modifier is exported as a
    # sensitivity feature but does not silently change this index.
    devfun_water_0_28 = devfun * WSFD_0_28_proxy,
    # changed to 0-100 cm as water limited function performed worse in benchmark using shallow proxy
    devfun_water_0_100 =devfun * WSFD_0_100_proxy
  ) %>%
    group_modify(~ {
      water_phase_state_0_28 <-
        calc_soltani_phase_water_state(
          tempfun = .x$tempfun,
          ppfun  = .x$ppfun,
          verfun = .x$verfun,
          wsfd   = .x$WSFD_0_28_proxy
        ) %>%
        rename_with(~ paste0(.x, "_0_28"))
      
      water_phase_state_0_100 <-
        calc_soltani_phase_water_state(
          tempfun = .x$tempfun,
          ppfun  = .x$ppfun,
          verfun = .x$verfun,
          wsfd   = .x$WSFD_0_100_proxy
        ) %>%
        rename_with(~ paste0(.x, "_0_100"))
      
      bind_cols(
        .x,
        water_phase_state_0_28,
        water_phase_state_0_100
      )
    }) %>%
    # -------------------------------------------------------------------------
  # 5d. Cumulative states from sowing
  # -------------------------------------------------------------------------
  mutate(
    cumdevfun = cumsum_complete_to_date(devfun),
    cumdevfun_phase = cumsum_complete_to_date(devfun_phase),
    cumdevfun_water_0_28 = cumsum_complete_to_date(devfun_water_0_28),
    cumdevfun_water_0_100 = cumsum_complete_to_date(devfun_water_0_100),
    cumdevfun_phase_water_0_28 = cumsum_complete_to_date(devfun_phase_water_0_28),
    cumdevfun_phase_water_0_100 = cumsum_complete_to_date(devfun_phase_water_0_100),
    
    cumTmax = cumsum_available(tmax),
    cumTmin = cumsum_available(tmin),
    cumTemp = cumsum_available(tmean),
    cumtempfun = cumsum_available(tempfun),
    cumppfun = cumsum_available(ppfun),
    cumppfun_phase = cumsum_available(ppfun_phase),
    cumvernfun = cumsum_available(verfun),
    cumvernfun_phase = cumsum_available(verfun_phase),
    
    cumRR = cumsum_available(prec),
    cumRad = cumsum_available(rad),
    cumCDD = cumsum_available(CDD),
    cumHDD = cumsum_available(HDD),
    cumColdDays = cumsum_available(cold_day),
    cumHotDays = cumsum_available(hot_day),
    
    n_tmax = cumsum(!is.na(tmax)),
    n_tmin = cumsum(!is.na(tmin)),
    n_tmean = cumsum(!is.na(tmean)),
    n_tempfun = cumsum(!is.na(tempfun)),
    n_pp = cumsum(!is.na(ppfun)),
    n_pp_phase = cumsum(!is.na(ppfun_phase)),
    n_vern = cumsum(!is.na(verfun)),
    n_vern_phase = cumsum(!is.na(verfun_phase)),
    n_dev = cumsum(!is.na(devfun)),
    n_dev_phase = cumsum(!is.na(devfun_phase)),
    n_dev_water_0_28 = cumsum(!is.na(devfun_water_0_28)),
    n_dev_water_0_100 = cumsum(!is.na(devfun_water_0_100)),
    n_dev_phase_water_0_28 = cumsum(!is.na(devfun_phase_water_0_28)),
    n_dev_phase_water_0_100 = cumsum(!is.na(devfun_phase_water_0_100))
  ) %>%
    
    # -------------------------------------------------------------------------
  # 5e. To-date means: retained for explicit ablation experiments
  # -------------------------------------------------------------------------
  mutate(
    meanTmax_to_date = if_else(n_tmax > 0, cumTmax / n_tmax, NA_real_),
    meanTmin_to_date = if_else(n_tmin > 0, cumTmin / n_tmin, NA_real_),
    meanTemp_to_date = if_else(n_tmean > 0, cumTemp / n_tmean, NA_real_),
    meantempfun_to_date = if_else(
      n_tempfun > 0,
      cumtempfun / n_tempfun,
      NA_real_
    ),
    meanppfun_to_date = if_else(n_pp > 0, cumppfun / n_pp, NA_real_),
    meanppfun_phase_to_date = if_else(
      n_pp_phase > 0,
      cumppfun_phase / n_pp_phase,
      NA_real_
    ),
    meanvernfun_to_date = if_else(
      n_vern > 0,
      cumvernfun / n_vern,
      NA_real_
    ),
    meanvernfun_phase_to_date = if_else(
      n_vern_phase > 0,
      cumvernfun_phase / n_vern_phase,
      NA_real_
    ),
    mean_devfun_to_date = if_else(
      n_dev > 0,
      cumdevfun / n_dev,
      NA_real_
    ),
    mean_devfun_phase_to_date = if_else(
      n_dev_phase > 0,
      cumdevfun_phase / n_dev_phase,
      NA_real_
    ),
    mean_devfun_water_0_28_to_date = if_else(
      n_dev_water_0_28 > 0,
      cumdevfun_water_0_28 / n_dev_water_0_28,
      NA_real_
    ),
    
    mean_devfun_water_0_100_to_date = if_else(
      n_dev_water_0_100 > 0,
      cumdevfun_water_0_100 / n_dev_water_0_100,
      NA_real_
    ),
    
    mean_devfun_phase_water_0_28_to_date = if_else(
      n_dev_phase_water_0_28 > 0,
      cumdevfun_phase_water_0_28 / n_dev_phase_water_0_28,
      NA_real_
    ),
    
    mean_devfun_phase_water_0_100_to_date = if_else(
      n_dev_phase_water_0_100 > 0,
      cumdevfun_phase_water_0_100 / n_dev_phase_water_0_100,
      NA_real_
    )
  ) %>%
    ungroup()
  
  # Detect the earlier SoilGrids/ERA5 unit mismatch before writing bad FTSW.
  n_extreme_ftsw <- sum(
    abs(df$FTSW_0_28_raw) > 10 | abs(df$FTSW_0_100_raw) > 10,
    na.rm = TRUE
  )
  if (n_extreme_ftsw > 0L) {
    stop(
      "Harvest year ", y, " contains ", n_extreme_ftsw,
      " rows with abs(FTSW) > 10 in the 0-28 or 0-100 cm profile. ",
      "Check soil-water units.",
      call. = FALSE
    )
  }
  
  assert_unique_keys(
    df,
    c("site_key", "season_year", "Date"),
    paste0("final daily output for ", y)
  )
  
  output_columns <- c(
    # Identifiers and time
    "site_key", "traj_id", "location_id", "location", "lat_orig",
    "lon_orig", "season_year",
    "Date", "DOY", "sow_doy", "sowing_date", "days_since_sowing",
    
    # Temperature, photoperiod, vernalization, and phase states
    "TMP", "daylength", "photoperiod_civil", "GDD", "DTU_beta",
    "cumGDD", "cumDTU_beta", "tempfun", "ppfun_raw",
    "ppfun", "ppfun_phase", "VERDAY", "cumVER", "vern_saturated",
    "verfun", "verfun_phase", "soltani_photo_active",
    "soltani_vern_active", "soltani_response_active",
    "ppfun_phase_water_0_28",
    "verfun_phase_water_0_28",
    "soltani_photo_active_water_0_28",
    "soltani_vern_active_water_0_28",
    "soltani_response_active_water_0_28",
    "ppfun_phase_water_0_100",
    "verfun_phase_water_0_100",
    "soltani_photo_active_water_0_100",
    "soltani_vern_active_water_0_100",
    "soltani_response_active_water_0_100",
    "post_emergence", "emerged_proxy",
    
    # All four daily and cumulative development variants
    "devfun", "devfun_phase", "devfun_water_0_28",
    "devfun_water_0_100", "devfun_phase_water_0_28",
    "devfun_phase_water_0_100",
    "cumdevfun", "cumdevfun_phase", "cumdevfun_water_0_28",
    "cumdevfun_water_0_100",
    "cumdevfun_phase_water_0_28",
    "cumdevfun_phase_water_0_100",
    
    # Cumulative and recent weather
    "cumTmax", "cumTmin", "cumTemp", "cumtempfun", "cumppfun",
    "cumppfun_phase", "cumvernfun", "cumvernfun_phase", "cumColdDays",
    "cumHotDays", "cumCDD", "cumHDD", "cumRad", "cumRR",
    "window_complete_14", "window_complete_30", "sumRR_14",
    "sumRR_30", "sumRad_14", "sumRad_30", "sumCDD_14", "sumHDD_14", "sumCDD_30", "sumHDD_30",
    
    # Direct water features
    "e", "pev", "evaporation_mm", "potential_evaporation_mm", "swvl1",
    "swvl2", "swvl3", "swvl4", "soil_water_0_28_vwc",
    "soil_water_0_28_mm", "soil_water_0_100_vwc", "soil_water_0_100_mm",
    "soil_water_0_289_vwc", "soil_water_0_289_mm",
    "swvl1_minus_swvl2", "swvl2_minus_swvl3", "swvl3_minus_swvl4",
    "et_gap_mm", "et_fraction", "sumEvap_14", "sumPEvap_14",
    "sumEtGap_14", "meanSWVL1_14", "meanSWVL2_14", "meanSWVL3_14",
    "meanSWVL4_14", "meanSoilWater_0_28_14",
    "meanSoilWater_0_100_14", "nSoilWater_14",
    "nSoilWater_0_100_14", "nEvap_14", "nPEvap_14", "etFraction_14",
    "climaticWaterBalance_14", "sumEvap_30", "sumPEvap_30",
    "sumEtGap_30", "meanSWVL1_30", "meanSWVL2_30", "meanSWVL3_30",
    "meanSWVL4_30", "meanSoilWater_0_28_30",
    "meanSoilWater_0_100_30", "nSoilWater_30",
    "nSoilWater_0_100_30", "nEvap_30", "nPEvap_30", "etFraction_30",
    "climaticWaterBalance_30", "soilWaterChange_7",
    "soilWaterChange_14", "soilWaterChange_0_100_30",
    "soilWaterChange_0_100_14",
    
    # Static hydraulic limits and process-water state
    "standort", "soil_dist_m", "hydraulic_source_0_28",
    "hydraulic_source_0_100", "hydraulic_qc_p33_0_28",
    "hydraulic_qc_p10_0_28", "hydraulic_qc_p33_0_100",
    "hydraulic_qc_p10_0_100", "hydraulic_high_uncertainty",
    "hydraulic_max_cv_0_100", "theta_fc33_l1", "theta_fc33_l2",
    "theta_fc33_l3", "theta_fc10_l1", "theta_fc10_l2",
    "theta_fc10_l3", "theta_wp_l1", "theta_wp_l2", "theta_wp_l3",
    "theta_fc_0_28", "theta_fc33_0_28", "theta_fc10_0_28",
    "theta_wp_0_28", "theta_fc_0_100", "theta_fc33_0_100",
    "theta_fc10_0_100", "theta_wp_0_100", "paw_fc33_0_28_mm",
    "paw_fc10_0_28_mm", "paw_fc33_0_100_mm", "paw_fc10_0_100_mm",
    "hydraulic_limits_0_28_available",
    "hydraulic_limits_p10_0_28_available",
    "hydraulic_limits_0_100_available",
    "hydraulic_limits_p10_0_100_available",
    "available_soil_water_0_28_mm", "total_available_water_p33_0_28_mm",
    "total_available_water_p10_0_28_mm", "FTSW_0_28_raw",
    "FTSW_0_28_proxy", "FTSW_0_28_p10_raw", "FTSW_0_28_p10_proxy",
    "available_soil_water_0_100_mm", "total_available_water_p33_0_100_mm",
    "total_available_water_p10_0_100_mm", "FTSW_0_100_raw",
    "FTSW_0_100_proxy", "FTSW_0_100_p10_raw", "FTSW_0_100_p10_proxy",
    "WSFG_0_28_proxy", "WSFG_0_100_proxy", "WSFD_0_28_proxy",
    "WSFD_0_100_proxy",
    
    # To-date mean features and their availability counts
    "meanTmax_to_date", "meanTmin_to_date", "meanTemp_to_date",
    "meantempfun_to_date", "meanppfun_to_date",
    "meanppfun_phase_to_date", "meanvernfun_to_date",
    "meanvernfun_phase_to_date", "mean_devfun_to_date",
    "mean_devfun_phase_to_date","mean_devfun_water_0_28_to_date",
    "mean_devfun_water_0_100_to_date",
    "mean_devfun_phase_water_0_28_to_date",
    "mean_devfun_phase_water_0_100_to_date", "n_tmax", "n_tmin", "n_tmean",
    "n_tempfun", "n_pp", "n_pp_phase", "n_vern", "n_vern_phase",
    "n_dev", "n_dev_phase", "n_dev_water_0_28",
    "n_dev_water_0_100",
    "n_dev_phase_water_0_28",
    "n_dev_phase_water_0_100"
  )
  
  fwrite(
    df %>% select(all_of(output_columns)),
    file.path(GDD_OUTPUT_DIR, paste0("gdd_harv_", y, ".csv")),
    na = "NA"
  )
  
  rm(list = intersect(c("weather_raw", "weather_tables", "df"), ls()))
  gc()
  invisible(NULL)
}


# =============================================================================
# 6. Run yearly calculations
# =============================================================================

dir.create(GDD_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

years <- sort(unique(sow_harv$harv_year))

plan(multisession, workers = N_WORKERS)
future_lapply(years, calc_gdd_winter_geo,future.seed = T)
plan(sequential)


# =============================================================================
# 7. Fast, single-pass merge back to irregular phenology observations
# =============================================================================

phen_df <- phen_df %>%
  mutate(phen_date = as.Date(phen_date))

phen_by_year <- split(phen_df, phen_df$year)
phen_list <- setNames(vector("list", length(years)), as.character(years))

for (y in years) {
  
  message("Merging year: ", y)
  
  gdd_file <- file.path(GDD_OUTPUT_DIR, paste0("gdd_harv_", y, ".csv"))
  if (!file.exists(gdd_file)) {
    warning("Skipping year ", y, ": GDD file is missing.")
    next
  }
  
  phen_y <- phen_by_year[[as.character(y)]]
  if (is.null(phen_y) || nrow(phen_y) == 0L) next
  
  gdd <- fread(gdd_file, showProgress = FALSE) %>%
    tibble::as_tibble() %>%
    mutate(
      Date = as.Date(Date),
      season_year = as.integer(season_year)
    )
  
  assert_unique_keys(
    gdd,
    c("site_key", "traj_id", "season_year", "Date"),
    paste0("GDD output for ", y)
  )
  
  phen_list[[as.character(y)]] <- phen_y %>%
    mutate(season_year = as.integer(season_year)) %>%
    left_join(
      gdd,
      by = c(
        "site_key",
        "traj_id",
        "season_year",
        "phen_date" = "Date"
      ),
      suffix = c(".pheno", ".gdd"),
      relationship = "many-to-one"
    ) %>%
    mutate(
      location_id = coalesce(
        as.character(location_id.pheno),
        as.character(location_id.gdd)
      ),
      location = coalesce(location.pheno, location.gdd),
      lat_orig = coalesce(lat_orig.pheno, lat_orig.gdd),
      lon_orig = coalesce(lon_orig.pheno, lon_orig.gdd),
      DOY = coalesce(DOY.pheno, DOY.gdd)
    ) %>%
    select(-ends_with(".pheno"), -ends_with(".gdd"))
}

phen_df_with_GDD <- bind_rows(phen_list) %>%
  mutate(
    across(where(is.list), unlist),
    lat_orig = as.numeric(lat_orig),
    lon_orig = as.numeric(lon_orig),
    year = as.integer(year),
    season_year = as.integer(season_year)
  )

phen_df_with_GDD2 <- phen_df_with_GDD %>%
  filter(!is.na(year), !is.na(phen_date))


# =============================================================================
# 8. QA checks and final output
# =============================================================================

cumulative_dev_columns <- c(
  "cumDTU_beta",
  "cumdevfun",
  "cumdevfun_phase",
  "cumdevfun_water_0_28",
  "cumdevfun_water_0_100",
  "cumdevfun_phase_water_0_28",
  "cumdevfun_phase_water_0_100"
)

bad_monotonic <- phen_df_with_GDD2 %>%
  arrange(traj_id, season_year, phen_date) %>%
  group_by(traj_id, season_year) %>%
  summarise(
    across(
      all_of(cumulative_dev_columns),
      ~ {
        values <- .x[!is.na(.x)]
        length(values) < 2L || all(diff(values) >= -1e-10)
      },
      .names = "monotonic_{.col}"
    ),
    .groups = "drop"
  ) %>%
  filter(if_any(starts_with("monotonic_"), ~ !.x))

if (nrow(bad_monotonic) > 0L) {
  warning(
    "Some site-season cumulative states are not monotonic. Inspect bad_monotonic."
  )
}

bounded_01_columns <- c(
  "tempfun", "ppfun", "ppfun_phase", "ppfun_phase_water_0_28","ppfun_phase_water_0_100", "VERDAY",
  "verfun", "verfun_phase", "verfun_phase_water_0_28", "verfun_phase_water_0_100",
  "FTSW_0_28_proxy", "FTSW_0_100_proxy",
  "WSFG_0_28_proxy", "WSFG_0_100_proxy"
)

for (column in intersect(bounded_01_columns, names(phen_df_with_GDD2))) {
  values <- phen_df_with_GDD2[[column]]
  if (any(values < -1e-10 | values > 1 + 1e-10, na.rm = TRUE)) {
    stop(column, " contains values outside [0, 1].", call. = FALSE)
  }
}

for (
  column in intersect(
    c("WSFD_0_28_proxy", "WSFD_0_100_proxy"),
    names(phen_df_with_GDD2)
  )
) {
  values <- phen_df_with_GDD2[[column]]
  if (
    any(
      values < 1 - 1e-10 | values > 1 + WSSD + 1e-10,
      na.rm = TRUE
    )
  ) {
    stop(column, " contains values outside [1, 1 + WSSD].")
  }
}

median_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
}

mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

process_stage_summary <- phen_df_with_GDD2 %>%
  filter(!is.na(cumDTU_beta)) %>%
  group_by(BBCH) %>%
  summarise(
    n = n(),
    median_cumDTU_beta = median_or_na(cumDTU_beta),
    median_cumdevfun_phase = median_or_na(cumdevfun_phase),
    median_ppfun_phase = median_or_na(ppfun_phase),
    median_verfun_phase = median_or_na(verfun_phase),
    proportion_vernalized = mean_or_na(vern_saturated),
    median_FTSW_0_28 = median_or_na(FTSW_0_28_proxy),
    median_FTSW_0_100 = median_or_na(FTSW_0_100_proxy),
    median_WSFD_0_28 = median_or_na(WSFD_0_28_proxy),
    median_WSFD_0_100 = median_or_na(WSFD_0_100_proxy),
    proportion_complete_window_14 = mean_or_na(window_complete_14),
    proportion_complete_window_30 = mean_or_na(window_complete_30),
    .groups = "drop"
  ) %>%
  arrange(BBCH)

print(process_stage_summary, n = nrow(process_stage_summary))
write_csv(
  process_stage_summary,
  sub("\\.csv$", "_process_qc.csv", FINAL_OUTPUT_FILE)
)
write_csv(phen_df_with_GDD2, FINAL_OUTPUT_FILE)

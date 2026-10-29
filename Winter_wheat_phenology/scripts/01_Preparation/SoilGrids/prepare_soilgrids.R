# Prepare static SoilGrids predictors and ERA5-Land layer-matched water limits.
# Run from the repository root or set PHENOLOGY_PROJECT_DIR.
library(janitor)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
PROJECT_DIR <- normalizePath(
  Sys.getenv("PHENOLOGY_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
DATA_DIR <- file.path(PROJECT_DIR, "Data")
SOILGRIDS_FILE <- file.path(
  DATA_DIR, "SoilGrids", "SoilGrids_Final_Fixed_NoNulls_17_06.csv"
)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# SoilGrids hydraulic variables are static soil properties. They must be joined
# once by site, not repeated or accumulated by year.
soil_data <- read_csv(
  SOILGRIDS_FILE,
  show_col_types = FALSE
) %>%
  select(-any_of("system:index")) %>%
  clean_names()

soil_data <- soil_data %>%
  mutate(
    coords = str_extract(geo, "\\[.*\\]"),
    coords = str_remove_all(coords, "\\[|\\]"),
    lon = as.numeric(str_split(coords, ",", simplify = TRUE)[, 1]),
    lat = as.numeric(str_split(coords, ",", simplify = TRUE)[, 2])
  ) %>%
  select(-coords, -geo) %>%
  rename_with(~ str_replace(.x, "_mean_mean$", "_mean"))

# Replace sentinel missing values and standardise depth/statistic suffixes.
c_tbl2 <- soil_data %>%
  mutate(across(where(is.numeric), ~ na_if(.x, -999))) %>%
  rename_with(
    ~ .x %>%
      str_replace_all("-", "_") %>%
      str_replace(
        "_(0_5cm|5_15cm|15_30cm|30_60cm|60_100cm|100_200cm)_mean_(true|mean|cv|std_dev)$",
        "_\\1_\\2"
      ) %>%
      str_replace(
        "_(0_5cm|5_15cm|15_30cm|30_60cm|60_100cm|100_200cm)_mean$",
        "_\\1_mean"
      )
  )

depth_weights <- c(
  "0-5cm" = 5,
  "5-15cm" = 10,
  "15-30cm" = 15,
  "30-60cm" = 30,
  "60-100cm" = 40,
  "100-200cm" = 100
)

# The current input table stores water-retention values around 0.003 rather
# than volumetric fractions around 0.3. Multiplication by 100 restores m3 m-3.
# This factor is specific to the current export and is checked below against
# physically plausible converted values; do not transfer it to another export
# without inspecting the raw value distribution first.
SOILGRIDS_WATER_SCALE <- 100

wmean_na <- function(x, w) {
  ok <- !is.na(x) & !is.na(w)
  if (!any(ok)) return(NA_real_)
  weighted.mean(x[ok], w[ok])
}

# Pivot the depth-specific property columns; keep site coordinates as keys.
property_cols <- names(c_tbl2)[
  str_detect(
    names(c_tbl2),
    "_(0_5cm|5_15cm|15_30cm|30_60cm|60_100cm|100_200cm)_(mean|cv)$"
  )
]

if (length(property_cols) == 0L) {
  stop("No SoilGrids depth-specific mean/CV columns were detected.")
}

soil_long <- c_tbl2 %>%
  select(standort, lon, lat, all_of(property_cols)) %>%
  pivot_longer(
    cols = all_of(property_cols),
    names_to = c("trait", "depth", "stat"),
    names_pattern = paste0(
      "^(.*)_(0_5cm|5_15cm|15_30cm|30_60cm|60_100cm|100_200cm)",
      "_(mean|cv)$"
    ),
    values_to = "value"
  ) %>%
  mutate(
    weight_cm = unname(depth_weights[str_replace_all(depth, "_", "-")])
  )

# -----------------------------------------------------------------------------
# Existing 0-200 cm static soil-feature table (kept for compatibility)
# -----------------------------------------------------------------------------
soil_0_200 <- soil_long %>%
  group_by(standort, lon, lat, trait, stat) %>%
  summarise(
    value_0_200 = wmean_na(value, weight_cm),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value_0_200,
    names_prefix = "w_"
  ) %>%
  rename(
    mean_0_200 = w_mean,
    cv_0_200 = w_cv
  )

soil_cv_max <- soil_long %>%
  filter(stat == "cv") %>%
  group_by(standort, trait) %>%
  summarise(
    max_cv = if (all(is.na(value))) NA_real_ else max(value, na.rm = TRUE),
    .groups = "drop"
  )

soil_0_200_full <- soil_0_200 %>%
  left_join(soil_cv_max, by = c("standort", "trait"))

sites_to_remove_0_200 <- soil_0_200_full %>%
  filter(
    !is.na(cv_0_200),
    !is.na(max_cv),
    cv_0_200 >= 0.15 | max_cv >= 0.20
  ) %>%
  distinct(standort)

soil_0_200_full_wide <- soil_0_200_full %>%
  filter(
    !standort %in% sites_to_remove_0_200$standort,
    !trait %in% c("sand", "p33")
  ) %>%
  select(standort, lon, lat, trait, mean_0_200) %>%
  pivot_wider(
    id_cols = c(standort, lon, lat),
    names_from = trait,
    values_from = mean_0_200
  ) %>%
  mutate(
    # Retained in the original SoilGrids storage units for compatibility.
    swrc_0_200 = p10 - p1500
  )

write_csv(soil_0_200_full_wide, file.path(DATA_DIR, "Soil_data_16_07.csv"))

# -----------------------------------------------------------------------------
# Layer-matched hydraulic limits for ERA5-Land SWVL1-SWVL3 (0-100 cm)
# -----------------------------------------------------------------------------
# SoilGrids water-retention values are converted to the m3 m-3 units used by
# ERA5-Land soil water using the explicit scale defined above.
required_hydraulic_traits <- c("p10", "p33", "p1500")
required_hydraulic_depths <- c(
  "0_5cm", "5_15cm", "15_30cm", "30_60cm", "60_100cm"
)

missing_traits <- setdiff(required_hydraulic_traits, unique(soil_long$trait))
missing_depths <- setdiff(required_hydraulic_depths, unique(soil_long$depth))

if (length(missing_traits) > 0L || length(missing_depths) > 0L) {
  stop(
    "Cannot build ERA5-matched hydraulic limits. Missing traits: ",
    paste(missing_traits, collapse = ", "),
    "; missing depths: ", paste(missing_depths, collapse = ", ")
  )
}

# Exact overlap between the SoilGrids depth intervals and ERA5-Land layers.
# ERA5 L1 (0-7 cm) = 5 cm of SG 0-5 + 2 cm of SG 5-15.
# ERA5 L2 (7-28 cm) = 8 cm of SG 5-15 + 13 cm of SG 15-30.
# ERA5 L3 (28-100 cm) = 2 cm of SG 15-30 + 30 cm of SG 30-60 +
# 40 cm of SG 60-100.
era5_soilgrids_overlap <- tribble(
  ~era5_layer, ~depth,    ~overlap_cm,
  "l1",         "0_5cm",   5,
  "l1",         "5_15cm",  2,
  "l2",         "5_15cm",  8,
  "l2",         "15_30cm", 13,
  "l3",         "15_30cm", 2,
  "l3",         "30_60cm", 30,
  "l3",         "60_100cm", 40
)

hydraulic_layer_long <- soil_long %>%
  filter(
    stat == "mean",
    trait %in% required_hydraulic_traits,
    depth %in% required_hydraulic_depths
  ) %>%
  inner_join(era5_soilgrids_overlap, by = "depth") %>%
  mutate(theta_m3_m3 = value * SOILGRIDS_WATER_SCALE) %>%
  group_by(standort, lon, lat, trait, era5_layer) %>%
  summarise(
    theta_m3_m3 = wmean_na(theta_m3_m3, overlap_cm),
    .groups = "drop"
  ) %>%
  mutate(
    hydraulic_parameter = recode(
      trait,
      p10 = "theta_fc10",
      p33 = "theta_fc33",
      p1500 = "theta_wp"
    )
  )

hydraulic_uncertainty <- soil_long %>%
  filter(
    stat == "cv",
    trait %in% required_hydraulic_traits,
    depth %in% required_hydraulic_depths
  ) %>%
  group_by(standort, lon, lat) %>%
  summarise(
    hydraulic_max_cv_0_100 = if (all(is.na(value))) {
      NA_real_
    } else {
      max(value, na.rm = TRUE)
    },
    .groups = "drop"
  )

soil_hydraulic_era5_0_100 <- hydraulic_layer_long %>%
  select(
    standort, lon, lat, hydraulic_parameter, era5_layer, theta_m3_m3
  ) %>%
  pivot_wider(
    names_from = c(hydraulic_parameter, era5_layer),
    names_glue = "{hydraulic_parameter}_{era5_layer}",
    values_from = theta_m3_m3
  ) %>%
  left_join(hydraulic_uncertainty, by = c("standort", "lon", "lat")) %>%
  mutate(
    lat_key = round(lat, 5),
    lon_key = round(lon, 5),
    # Thickness-weighted 0-28 cm equivalents.
    theta_fc10_0_28 = (7 * theta_fc10_l1 + 21 * theta_fc10_l2) / 28,
    theta_fc33_0_28 = (7 * theta_fc33_l1 + 21 * theta_fc33_l2) / 28,
    theta_wp_0_28 = (7 * theta_wp_l1 + 21 * theta_wp_l2) / 28,
    # Thickness-weighted 0-100 cm equivalents.
    theta_fc10_0_100 = (
      7 * theta_fc10_l1 + 21 * theta_fc10_l2 + 72 * theta_fc10_l3
    ) / 100,
    theta_fc33_0_100 = (
      7 * theta_fc33_l1 + 21 * theta_fc33_l2 + 72 * theta_fc33_l3
    ) / 100,
    theta_wp_0_100 = (
      7 * theta_wp_l1 + 21 * theta_wp_l2 + 72 * theta_wp_l3
    ) / 100,
    paw_fc10_0_28_mm =
      70 * (theta_fc10_l1 - theta_wp_l1) +
      210 * (theta_fc10_l2 - theta_wp_l2),
    paw_fc33_0_28_mm =
      70 * (theta_fc33_l1 - theta_wp_l1) +
      210 * (theta_fc33_l2 - theta_wp_l2),
    paw_fc10_0_100_mm =
      70 * (theta_fc10_l1 - theta_wp_l1) +
      210 * (theta_fc10_l2 - theta_wp_l2) +
      720 * (theta_fc10_l3 - theta_wp_l3),
    paw_fc33_0_100_mm =
      70 * (theta_fc33_l1 - theta_wp_l1) +
      210 * (theta_fc33_l2 - theta_wp_l2) +
      720 * (theta_fc33_l3 - theta_wp_l3),
    hydraulic_qc_p10_0_28 =
      theta_fc10_l1 > theta_wp_l1 & theta_fc10_l2 > theta_wp_l2,
    hydraulic_qc_p33_0_28 =
      theta_fc33_l1 > theta_wp_l1 & theta_fc33_l2 > theta_wp_l2,
    hydraulic_qc_p10_0_100 =
      hydraulic_qc_p10_0_28 & theta_fc10_l3 > theta_wp_l3,
    hydraulic_qc_p33_0_100 =
      hydraulic_qc_p33_0_28 & theta_fc33_l3 > theta_wp_l3,
    hydraulic_high_uncertainty =
      !is.na(hydraulic_max_cv_0_100) & hydraulic_max_cv_0_100 >= 0.20
  )

if (anyDuplicated(soil_hydraulic_era5_0_100[c("lat_key", "lon_key")])) {
  stop(
    "The hydraulic table contains duplicate rounded coordinate keys. ",
    "Join by an explicit site identifier instead of coordinates."
  )
}



theta_cols <- names(soil_hydraulic_era5_0_100)[
  str_detect(names(soil_hydraulic_era5_0_100), "^theta_")
]

theta_values <- unlist(
  soil_hydraulic_era5_0_100[theta_cols],
  use.names = FALSE
)

write_csv(
  soil_hydraulic_era5_0_100,
  file.path(DATA_DIR, "Soil_hydraulic_ERA5_0_100.csv")
)

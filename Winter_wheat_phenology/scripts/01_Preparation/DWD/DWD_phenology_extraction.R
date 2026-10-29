# Harmonise DWD annual crop-reporter observations and record QC decisions.
#
# Run from the repository root, or set PHENOLOGY_PROJECT_DIR to that directory.
library(data.table)
library(dplyr)
library(lubridate)
library(readr)
library(tibble)

PROJECT_DIR <- normalizePath(
  Sys.getenv("PHENOLOGY_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
OUTPUT_DIR <- file.path(PROJECT_DIR, "Data", "DWD")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

safe_median_date <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) as.Date(NA) else median(x)
}

# DWD phenology and station metadata
hist_url <- "https://opendata.dwd.de/climate_environment/CDC/observations_germany/phenology/annual_reporters/crops/historical/PH_Jahresmelder_Landwirtschaft_Kulturpflanze_Winterweizen_1925_2024_hist.txt"
recent_url <- "https://opendata.dwd.de/climate_environment/CDC/observations_germany/phenology/annual_reporters/crops/recent/PH_Jahresmelder_Landwirtschaft_Kulturpflanze_Winterweizen_akt.txt"
station <- fread("https://opendata.dwd.de/climate_environment/CDC/help/PH_Beschreibung_Phaenologie_Stationen_Jahresmelder.txt", sep=";")

hist   <- fread(hist_url,   sep = ";", na.strings = c("-999", "-9999"))
recent <- fread(recent_url,  sep = ";", na.strings = c("-999", "-9999"))

pheno_ger <- rbindlist(list(hist, recent), fill = TRUE)

pheno_model <- merge(pheno_ger, station, by = "Stations_id", all.x = TRUE) %>%
  filter( # Winter wheat and the highest available quality category
    Objekt_id == 202,
    Qualitaetsniveau == 10,
    Eintrittsdatum_QB == 1,
    Referenzjahr >= 1980
  ) %>%
  mutate(
    Date = as.Date(sprintf("%08d", as.integer(Eintrittsdatum)), format = "%Y%m%d"),
    DOY  = yday(Date)
  ) %>%
  select(
    Stations_id, Referenzjahr, Phase_id, Date, DOY,
    geograph.Breite, geograph.Laenge, Stationsname
  ) %>%
  rename(lat = geograph.Breite, lon = geograph.Laenge)

bbch_lookup <- tribble(
  ~Phase_id, ~BBCH,
  10, 0,   # Bestellung Beginn
  12, 10,  # Auflaufen Beginn
  15, 31,  # Schossen Beginn
  18, 51,  # Ährenschieben Beginn
  19, 75,  # Milchreife Beginn
  21, 87,  # Gelbreife Beginn
  22, 89,  # Vollreife Beginn
  23, 99,  # Ernte von Hand
  24, 99   # Ernte
)

pheno_model <- pheno_model %>%
  left_join(bbch_lookup, by = "Phase_id") %>%
  filter(!is.na(BBCH))

# Diagnose the interval between yellow ripeness and harvest.
harvest_window <- pheno_model %>%
  filter(BBCH %in% c(87, 89, 99)) %>%
  group_by(Stations_id, Referenzjahr) %>%
  summarise(
    date_87 = safe_median_date(Date[BBCH == 87]),
    date_89 = safe_median_date(Date[BBCH == 89]),
    date_99 = safe_median_date(Date[BBCH == 99]),
    diff_days_87_99 = as.numeric(date_99 - date_87),
    diff_days_87_89 = as.numeric(date_89 - date_87),
    .groups = "drop"
  )

message("Median 87->99 days: ", median(harvest_window$diff_days_87_99, na.rm = TRUE))
message("Median 87->89 days: ", median(harvest_window$diff_days_87_89, na.rm = TRUE))

# DWD reference years denote harvest years. Correct the calendar year of
# autumn sowing and emergence observations before constructing trajectories.

pheno_model <- pheno_model %>%
  mutate(
    year_adj = case_when(
      BBCH == 0 ~ Referenzjahr - 1L,
      BBCH == 10 & DOY > 200 ~ Referenzjahr - 1L,
      TRUE ~ Referenzjahr
    ),
    Date_adj = as.Date(sprintf("%04d-%03d", year_adj, DOY), format = "%Y-%j")
  )
# Quality-control helpers
safe_min_date <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) as.Date(NA) else min(x)
}

fix_reversal_block <- function(df, max_gap_days = 100) {
  df <- df %>% arrange(season_date)
  
  df$BBCH_corrected <- df$BBCH
  df$reversal_action <- "kept"
  
  if (nrow(df) < 2) return(df)
  
  for (i in 2:nrow(df)) {
    dd <- as.numeric(df$season_date[i] - df$season_date[i - 1])
    
    if (is.na(dd) || dd <= 0) next
    
    prev_bbch <- df$BBCH_corrected[i - 1]
    cur_bbch  <- df$BBCH_corrected[i]
    
    if (is.na(prev_bbch) || is.na(cur_bbch)) next
    
    # local reversal
    if (cur_bbch < prev_bbch) {
      if (dd <= max_gap_days) {
        tmp <- df$BBCH_corrected[i - 1]
        df$BBCH_corrected[i - 1] <- df$BBCH_corrected[i]
        df$BBCH_corrected[i] <- tmp
        df$reversal_action[c(i - 1, i)] <- "swapped_reversal"
      } else {
        df$reversal_action[i] <- "manual_review_long_gap_reversal"
      }
    }
  }
  
  df
}
# Stable trajectory identifier and input to quality control
dat_raw <- pheno_model %>%
  distinct() %>%
  mutate(
    traj_id = interaction(Stations_id, Referenzjahr, drop = TRUE),
    season_date = Date_adj,
    SeasonYear = as.integer(Referenzjahr)
  ) %>%
  filter(BBCH < 90)   # harvest excluded

# Summarise observation support and remove trajectories without sowing.
traj_summary <- dat_raw %>%
  group_by(traj_id) %>%
  summarise(
    n_obs = n(),
    n_bbch = n_distinct(BBCH),
    bbch_values = paste(sort(unique(BBCH)), collapse = ", "),
    has_sowing = any(BBCH == 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    traj_class = case_when(
      n_bbch %in% c(1,2) ~ "too_sparse",
      n_bbch %in% c(3,4) ~ "sparse_keep_flagged",
      TRUE ~ "adequate"
    )
  )

traj_keep <- traj_summary %>%
  filter(has_sowing) %>%
  pull(traj_id)

dat_raw2 <- dat_raw %>%
  filter(traj_id %in% traj_keep) %>%
  left_join(traj_summary, by = "traj_id")

# Collapse same-day duplicates; retain BBCH 89 for same-day 87/89 pairs.
daily <- dat_raw2 %>%
  filter(!is.na(season_date), !is.na(BBCH)) %>%
  group_by(traj_id, Stations_id, Referenzjahr, SeasonYear, season_date) %>%
  summarise(
    n_day_obs = n(),
    n_day_bbch = n_distinct(BBCH),
    day_bbch_values = paste(sort(unique(BBCH)), collapse = ", "),
    BBCH = case_when(
      n_day_bbch == 1L ~ first(BBCH),
      setequal(unique(BBCH), c(87L, 89L)) ~ 89L,
      TRUE ~ NA_integer_
    ),
    lat = first(lat),
    lon = first(lon),
    Stationsname = first(Stationsname),
    traj_class = first(traj_class),
    has_sowing = first(has_sowing),
    n_obs = first(n_obs),
    n_bbch = first(n_bbch),
    bbch_values = first(bbch_values),
    .groups = "drop"
  )

same_day_conflicts <- daily %>%
  filter(is.na(BBCH) & n_day_bbch > 1L)

daily_clean <- daily %>%
  filter(!is.na(BBCH)) %>%
  arrange(traj_id, season_date)

# Swap local stage reversals up to 80 days apart; flag longer gaps.
traj_corrected <- daily_clean %>%
  group_by(traj_id) %>%
  group_modify(~ fix_reversal_block(.x, max_gap_days = 80)) %>%
  ungroup() 

# Check corrected trajectories.
traj_flags2 <- traj_corrected %>%
  group_by(traj_id) %>%
  arrange(season_date, .by_group = TRUE) %>%
  mutate(
    prev_date = lag(season_date),
    prev_BBCH = lag(BBCH_corrected),
    days_diff = as.numeric(season_date - prev_date),
    bbch_diff = BBCH_corrected - prev_BBCH,
    rate_per_day = if_else(is.na(days_diff) | days_diff == 0, NA_real_, bbch_diff / days_diff),
    flag_reversal = !is.na(bbch_diff) & bbch_diff < 0,
    flag_long_gap_reversal = flag_reversal & days_diff > 100,
    flag_large_jump_short_gap = !is.na(days_diff) & days_diff <= 3 & bbch_diff >= 20,
    flag_unusually_fast = !is.na(rate_per_day) & rate_per_day >= 10
  ) %>%
  ungroup()


# Preserve QC flags, including those for retained sparse trajectories.
pheno_qc <- traj_flags2 %>%
  mutate(
    qc_flag = case_when(
      is.na(season_date) | is.na(BBCH_corrected) ~ "missing_core_field",
      BBCH_corrected < 0 | BBCH_corrected > 99 ~ "invalid_bbch_range",
      flag_long_gap_reversal ~ "manual_review_long_gap_reversal",
      flag_large_jump_short_gap ~ "large_jump_short_gap",
      flag_unusually_fast ~ "unusually_fast",
      traj_class == "too_sparse" ~ "too_sparse",
      traj_class == "sparse_keep_flagged" ~ "sparse_keep_flagged",
      TRUE ~ "ok"
    )
  )

bad_ids <- pheno_qc %>%
  filter(qc_flag %in% c("too_sparse","manual_review_long_gap_reversal","invalid_bbch_range","missing_core_field"))


# Model-ready subset: remove the clearly problematic trajectories
pheno_model_ready <- pheno_qc %>%
  filter(!traj_id %in% bad_ids$traj_id)

# Sowing DOY for each trajectory
pheno_model_ready <- pheno_model_ready %>%
  group_by(traj_id) %>%
  mutate(
    sowing_DOY = {
      x <- if_else(BBCH_corrected == 0, yday(season_date), NA_integer_)
      if (all(is.na(x))) NA_integer_ else min(x, na.rm = TRUE)
    }
  ) %>%
  ungroup()

traj_summary_qc <- pheno_model_ready %>%
  group_by(traj_id) %>%
  summarise(
    n_obs_qc = n(),
    n_bbch_qc = n_distinct(BBCH_corrected),
    bbch_values_qc = paste(sort(unique(BBCH_corrected)), collapse = ", "),
    has_sowing_qc = any(BBCH_corrected == 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    traj_class_qc = case_when(
      n_bbch_qc %in% c(1, 2) ~ "too_sparse",
      n_bbch_qc %in% c(3, 4) ~ "sparse_keep_flagged",
      TRUE ~ "adequate"
    )
  )
bad_ids_after_qc <- traj_summary_qc %>%
  filter(traj_class_qc %in% c("too_sparse"))

pheno_model_ready <- pheno_model_ready %>% filter(!traj_id %in% bad_ids_after_qc$traj_id)
# Review tables and model-ready output
long_gap_review <- pheno_qc %>%
  filter(qc_flag == "manual_review_long_gap_reversal")

sparse_review <- traj_summary %>%
  filter(traj_class  %in% c("too_sparse", "sparse_keep_flagged"))

sparse_review2 <- traj_summary_qc %>%
  filter(traj_class_qc  %in% c("too_sparse", "sparse_keep_flagged"))
write_csv(pheno_model_ready, file.path(OUTPUT_DIR, "DWD_phenology_16_06_final.csv"))
write_csv(pheno_qc, file.path(OUTPUT_DIR, "DWD_phenology_qc_all.csv"))
write_csv(traj_summary, file.path(OUTPUT_DIR, "DWD_trajectory_summary.csv"))
write_csv(same_day_conflicts, file.path(OUTPUT_DIR, "DWD_same_day_conflicts.csv"))
write_csv(long_gap_review, file.path(OUTPUT_DIR, "DWD_long_gap_reversal_review.csv"))
write_csv(sparse_review, file.path(OUTPUT_DIR, "DWD_sparse_trajectories_review.csv"))
write_csv(sparse_review2, file.path(OUTPUT_DIR, "DWD_sparse_trajectories_after_qc_review.csv"))

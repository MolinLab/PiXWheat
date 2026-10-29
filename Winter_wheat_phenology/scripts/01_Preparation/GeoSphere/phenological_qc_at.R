# Harmonise Austrian PEP725, AGES and WheatVIZ phenology observations.
library(data.table)
library(dplyr)
library(lubridate)
library(readr)
library(stringr)
library(tibble)

PROJECT_DIR <- normalizePath(
  Sys.getenv("PHENOLOGY_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
PEP_DIR <- file.path(PROJECT_DIR, "Data", "PEP725")
GEOSPHERE_DIR <- file.path(PROJECT_DIR, "Data", "GeoSphere")
dir.create(PEP_DIR, recursive = TRUE, showWarnings = FALSE)

pep_geo <- read_csv(
  file.path(PEP_DIR, "pep725_with_geo_labels.csv"),
  show_col_types = FALSE
)

date_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) as.Date(NA) else min(x)
}

pep_audit <- pep_geo %>%
  filter(
    country == "Österreich",
    cult_season == 2,
    phase_id < 90,
    !(phase_id %in% c(0, 10) & day < 200)
  ) %>%
  mutate(
    raw_date = as.Date(date),
    raw_year = lubridate::year(raw_date),
    raw_doy  = lubridate::yday(raw_date),
    
    provisional_SeasonYear = case_when(
      phase_id == 0L ~ raw_year + 1L,
      phase_id == 10L & raw_doy > 200L ~ raw_year + 1L,
      TRUE ~ raw_year
    )
  ) %>%
  mutate(
    phase_id = case_when(
      location_id %in% c("48.7333_15.25", "47.1_15.7") &
        date %in% c("1990-10-30", "2011-10-14", "2011-11-08") ~ NA_real_,
      TRUE ~ phase_id
    )
  ) %>%
  filter(
    !location %in% c("Grabern", "Knittelfeld", "Kriechbaum", "Lutzmannsburg"),
    !is.na(phase_id)
  )

trajectory_audit <- pep_audit %>%
  group_by(location, provisional_SeasonYear) %>%
  summarise(
    n_obs = n(),
    n_bbch = n_distinct(phase_id),
    bbch_values = paste(sort(unique(phase_id)), collapse = ", "),
    early_only = all(phase_id %in% c(0L, 10L)),
    has_sowing = any(phase_id == 0L),
    has_emergence = any(phase_id == 10L),
    has_later = any(phase_id >= 31L),
    first_later_date = date_or_na(raw_date[phase_id >= 31L]),
    .groups = "drop"
  )

early_orphans <- trajectory_audit %>%
  filter(early_only, has_sowing) %>%
  transmute(
    location,
    orphan_SeasonYear = provisional_SeasonYear,
    target_SeasonYear = provisional_SeasonYear - 1L,
    orphan_stages = bbch_values
  )

previous_trajectories <- trajectory_audit %>%
  filter(has_later, !has_sowing) %>%
  transmute(
    location,
    target_SeasonYear = provisional_SeasonYear,
    later_stages = bbch_values,
    first_later_date
  )

year_offset_candidates <- early_orphans %>%
  inner_join(
    previous_trajectories,
    by = c("location", "target_SeasonYear")
  )

candidate_keys <- year_offset_candidates %>%
  distinct(
    location,
    orphan_SeasonYear,
    target_SeasonYear
  )

candidate_early_dates <- pep_audit %>%
  inner_join(
    candidate_keys,
    by = c(
      "location",
      "provisional_SeasonYear" = "orphan_SeasonYear"
    )
  ) %>%
  group_by(
    location,
    orphan_SeasonYear = provisional_SeasonYear,
    target_SeasonYear
  ) %>%
  summarise(
    n_sowing_obs = sum(phase_id == 0L),
    n_emergence_obs = sum(phase_id == 10L),
    
    n_sowing_dates = n_distinct(
      raw_date[phase_id == 0L],
      na.rm = TRUE
    ),
    
    n_emergence_dates = n_distinct(
      raw_date[phase_id == 10L],
      na.rm = TRUE
    ),
    
    raw_sowing_date = date_or_na(
      raw_date[phase_id == 0L]
    ),
    
    raw_emergence_date = date_or_na(
      raw_date[phase_id == 10L]
    ),
    
    .groups = "drop"
  )

target_keys <- candidate_keys %>%
  distinct(location, target_SeasonYear)

candidate_target_dates <- pep_audit %>%
  inner_join(
    target_keys,
    by = c(
      "location",
      "provisional_SeasonYear" = "target_SeasonYear"
    )
  ) %>%
  group_by(
    location,
    target_SeasonYear = provisional_SeasonYear
  ) %>%
  summarise(
    target_has_sowing = any(phase_id == 0L),
    target_has_emergence = any(phase_id == 10L),
    
    n_target_later_stages = n_distinct(
      phase_id[phase_id >= 31L & phase_id < 90L]
    ),
    
    target_later_stages = paste(
      sort(unique(phase_id[phase_id >= 31L & phase_id < 90L])),
      collapse = ", "
    ),
    
    first_later_date = date_or_na(
      raw_date[phase_id >= 31L & phase_id < 90L]
    ),
    
    last_later_date = {
      x <- raw_date[phase_id >= 31L & phase_id < 90L]
      x <- x[!is.na(x)]
      if (!length(x)) as.Date(NA) else max(x)
    },
    
    .groups = "drop"
  )

candidate_checks <- candidate_early_dates %>%
  left_join(
    candidate_target_dates,
    by = c("location", "target_SeasonYear")
  ) %>%
  mutate(
    proposed_sowing_date =
      raw_sowing_date %m-% years(1),
    
    proposed_emergence_date =
      raw_emergence_date %m-% years(1),
    
    # Were the suspicious dates recorded in the harvest year?
    raw_year_matches_target =
      lubridate::year(raw_sowing_date) == target_SeasonYear &
      lubridate::year(raw_emergence_date) == target_SeasonYear,
    
    # After correction they should belong to the preceding autumn.
    corrected_year_matches =
      lubridate::year(proposed_sowing_date) ==
      target_SeasonYear - 1L &
      lubridate::year(proposed_emergence_date) ==
      target_SeasonYear - 1L,
    
    sowing_to_emergence_days = as.numeric(
      proposed_emergence_date - proposed_sowing_date
    ),
    
    sowing_to_first_later_days = as.numeric(
      first_later_date - proposed_sowing_date
    ),
    
    emergence_to_first_later_days = as.numeric(
      first_later_date - proposed_emergence_date
    ),
    
    unique_early_dates =
      n_sowing_dates == 1L &
      n_emergence_dates == 1L,
    
    no_target_early_stages =
      !target_has_sowing &
      !target_has_emergence,
    
    chronology_correct =
      proposed_sowing_date <= proposed_emergence_date &
      proposed_emergence_date < first_later_date,
    
    sowing_emergence_duration_plausible =
      dplyr::between(
        sowing_to_emergence_days,
        0,
        90
      ),
    
    sowing_later_duration_plausible =
      dplyr::between(
        sowing_to_first_later_days,
        60,
        330
      ),
    
    emergence_later_duration_plausible =
      dplyr::between(
        emergence_to_first_later_days,
        30,
        300
      ),
    
    strong_later_stage_anchor =
      n_target_later_stages >= 2L
  )

candidate_checks <- candidate_checks %>%
  mutate(
    automatic_confirmation = coalesce(
      unique_early_dates &
        raw_year_matches_target &
        corrected_year_matches &
        no_target_early_stages &
        chronology_correct &
        sowing_emergence_duration_plausible &
        sowing_later_duration_plausible &
        emergence_later_duration_plausible &
        strong_later_stage_anchor,
      FALSE
    ),
    
    automatic_decision = case_when(
      automatic_confirmation ~
        "confirmed",
      
      target_has_sowing | target_has_emergence ~
        "rejected_target_already_has_early_stage",
      
      !is.na(raw_year_matches_target) &
        !raw_year_matches_target ~
        "rejected_raw_year_does_not_match",
      
      !is.na(chronology_correct) &
        !chronology_correct ~
        "rejected_implausible_chronology",
      
      !is.na(sowing_emergence_duration_plausible) &
        !sowing_emergence_duration_plausible ~
        "manual_review_sowing_emergence_duration",
      
      !is.na(sowing_later_duration_plausible) &
        !sowing_later_duration_plausible ~
        "manual_review_sowing_later_duration",
      
      !unique_early_dates ~
        "manual_review_missing_or_duplicate_early_date",
      
      n_target_later_stages < 2L ~
        "manual_review_weak_later_stage_anchor",
      
      TRUE ~
        "manual_review_other"
    )
  )


candidate_checks <- candidate_checks %>%
  mutate(
    manual_decision = case_when(
      str_starts(automatic_decision, "manual_review") ~ "confirmed",
      TRUE ~ automatic_decision
    )
  )

candidate_decisions <- candidate_checks %>% mutate(
  confirmed_year_offset = case_when(
    manual_decision == "confirmed" ~ TRUE,
    manual_decision != "confirmed"  ~ FALSE
  ),
  final_decision_source = case_when(
    manual_decision %in% c("confirmed", "reject") ~ "manual",
    automatic_confirmation ~ "automatic",
    TRUE ~ "rejected"
  )
)
repair_keys <- candidate_decisions %>%
  select(
    location,
    orphan_SeasonYear,
    confirmed_year_offset,
    automatic_decision,
    final_decision_source
  )

pep_corrected <- pep_audit %>%
  left_join(
    repair_keys,
    by = c(
      "location",
      "provisional_SeasonYear" = "orphan_SeasonYear"
    )
  ) %>%
  mutate(
    confirmed_year_offset =
      coalesce(confirmed_year_offset, FALSE),
    
    Date_original = raw_date,
    
    Date_corrected = case_when(
      confirmed_year_offset &
        phase_id %in% c(0L, 10L) ~
        raw_date %m-% years(1),
      
      TRUE ~ raw_date
    ),
    
    date_shift_years = if_else(
      confirmed_year_offset &
        phase_id %in% c(0L, 10L),
      -1L,
      0L
    ),
    
    date_repair_status = case_when(
      confirmed_year_offset ~ "shifted_back_one_year",
      !is.na(automatic_decision) ~ automatic_decision,
      TRUE ~ "not_candidate"
    ),
    
    # Recalculate these fields from the corrected date.
    year = lubridate::year(Date_corrected),
    day = lubridate::yday(Date_corrected),
    date = Date_corrected,
    
    SeasonYear = case_when(
      phase_id == 0L ~ year + 1L,
      phase_id == 10L & day > 200L ~ year + 1L,
      TRUE ~ year
    )
  )

pep_at <- pep_corrected %>%
  select(-county, -row_id, -license_short, -cult_season, -subspecies,
         -species, -genus, -gss_id, -alt_dem, -alt, -provider_id) %>%
  rename(BBCH = phase_id, Date = Date_corrected) %>%
  mutate(
    traj_id = as.character(interaction(location, SeasonYear, drop = TRUE)),
    source = "PEP"
  ) %>% mutate(DOY=yday(Date))

df_at <- read_csv(
  file.path(GEOSPHERE_DIR, "WVIZ_AGES_pheno.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    lat = case_when(
      Location_name == "Grabenegg, NÖ" & lat == 48.096902 ~ 48.131017,
      TRUE ~ lat
    ),
    lon = case_when(
      Location_name == "Grabenegg, NÖ" & lon == 14.896003 ~ 15.227549,
      TRUE ~ lon
    ),
    Location = coalesce(Location, Location_name),
    DOY = yday(Date)
  )

df_at <- df_at %>%
  rename(traj_id = env_id) %>%
  mutate(
    traj_id = as.character(traj_id),
    source = "AGES"
  )

pheno_model <- bind_rows(pep_at, df_at) %>%
  distinct() 

safe_min_date <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) as.Date(NA) else min(x)
}

fix_reversal_block <- function(df, bbch_col = "BBCH_adjusted", max_gap_days = 80) {
  df <- df %>% arrange(season_date)
  
  df$BBCH_corrected <- df[[bbch_col]]
  df$reversal_action <- "kept"
  
  if (nrow(df) < 2) return(df)
  
  for (i in 2:nrow(df)) {
    dd <- as.numeric(df$season_date[i] - df$season_date[i - 1])
    if (is.na(dd) || dd <= 0) next
    
    prev_bbch <- df$BBCH_corrected[i - 1]
    cur_bbch  <- df$BBCH_corrected[i]
    
    if (is.na(prev_bbch) || is.na(cur_bbch)) next
    
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

# ----------------------------
# Stable trajectory id and raw QC input
# ----------------------------
dat_raw <- pheno_model %>%
  distinct() %>%
  mutate(
    traj_id = as.character(traj_id),
    season_date = as.Date(Date),
    SeasonYear = as.integer(SeasonYear)
  ) %>%
  filter(BBCH < 90)   # harvest excluded

dat_raw <- dat_raw %>%
  filter(!(BBCH %in% c(0, 10) & DOY < 200))
# ----------------------------
# Trajectory summary
# Remove trajectories without sowing
# Flag sparse trajectories, but keep them for transparency
# ----------------------------
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

# ----------------------------
# Manual AT corrections
# ----------------------------
manual_fix <- tribble(
  ~traj_id, ~BBCH, ~BBCH_to, ~fix_type, ~fix_reason,
  
  "Wolfsegg am Hausruck.1990", 60L, 51L, "swap",   "60 and 51 reversed",
  "Wolfsegg am Hausruck.1990", 51L, 60L, "swap",   "60 and 51 reversed",
  
  "Merkendorf.1995",           89L, 60L, "swap",   "89 and 60 reversed",
  "Merkendorf.1995",           60L, 89L, "swap",   "89 and 60 reversed",
  "Merkendorf.1997",           89L, 60L, "swap",   "89 and 60 reversed",
  "Merkendorf.1997",           60L, 89L, "swap",   "89 and 60 reversed",
  
  "Millstatt am See.1995",     51L, 31L, "swap",   "observed 51,60,31; should be 31,51,60",
  "Millstatt am See.1995",     60L, 51L, "swap",   "observed 51,60,31; should be 31,51,60",
  "Millstatt am See.1995",     31L, 60L, "swap",   "observed 51,60,31; should be 31,51,60",
  
  "Friesach.1999",              0L, 10L, "swap",   "0 and 10 reversed",
  "Friesach.1999",             10L, 0L,  "swap",   "0 and 10 reversed",
  
  "Reichersberg.2002",         51L, 60L, "swap",   "51 and 60 reversed",
  "Reichersberg.2002",         60L, 51L, "swap",   "51 and 60 reversed",
  
  "Weitra.2003",               87L, 75L, "swap",   "87 and 75 reversed",
  "Weitra.2003",               75L, 87L, "swap",   "87 and 75 reversed",
  
  "Fürstenfeld.2000",          75L, NA_integer_, "remove", "same-day 75/87 conflict; remove 75"
)

# Apply manual fixes
traj_fixed <- dat_raw2 %>%
  left_join(manual_fix, by = c("traj_id", "BBCH")) %>%
  mutate(
    BBCH_adjusted = case_when(
      fix_type == "remove" ~ NA_integer_,
      !is.na(BBCH_to) ~ as.integer(BBCH_to),
      TRUE ~ as.integer(BBCH)
    ),
    BBCH_changed = if_else(is.na(BBCH_adjusted), TRUE, BBCH_adjusted != BBCH),
    fix_type = coalesce(fix_type, "none"),
    fix_reason = coalesce(fix_reason, "none")
  )

manual_fix_log <- traj_fixed %>%
  filter(fix_type != "none") %>%
  select(traj_id, season_date, BBCH, BBCH_adjusted, fix_type, fix_reason)

manual_removed_review <- traj_fixed %>%
  filter(is.na(BBCH_adjusted)) %>%
  select(traj_id, season_date, BBCH, fix_type, fix_reason)

traj_fixed_clean <- traj_fixed %>%
  filter(!is.na(BBCH_adjusted))

# ----------------------------
# QC on corrected trajectories
# ----------------------------
traj_flags2 <- traj_fixed_clean %>%
  group_by(traj_id) %>%
  arrange(season_date, .by_group = TRUE) %>%
  mutate(
    prev_date = lag(season_date),
    prev_BBCH = lag(BBCH_adjusted),
    days_diff = as.numeric(season_date - prev_date),
    bbch_diff = BBCH_adjusted - prev_BBCH,
    rate_per_day = if_else(is.na(days_diff) | days_diff == 0, NA_real_, bbch_diff / days_diff),
    flag_reversal = !is.na(bbch_diff) & bbch_diff < 0,
    flag_long_gap_reversal = flag_reversal & days_diff > 100,
    flag_large_jump_short_gap = !is.na(days_diff) & days_diff <= 3 & bbch_diff >= 20,
    flag_unusually_fast = !is.na(rate_per_day) & rate_per_day >= 10
  ) %>%
  ungroup()

# ----------------------------
# Final QC flags
# Keep sparse classes flagged; do not auto-delete them here
# ----------------------------
pheno_qc <- traj_flags2 %>%
  mutate(
    qc_flag = case_when(
      is.na(season_date) | is.na(BBCH_adjusted) ~ "missing_core_field",
      BBCH_adjusted < 0 | BBCH_adjusted > 99 ~ "invalid_bbch_range",
      flag_large_jump_short_gap ~ "large_jump_short_gap",
      flag_unusually_fast ~ "unusually_fast",
      traj_class == "too_sparse" ~ "too_sparse",
      traj_class == "sparse_keep_flagged" ~ "sparse_keep_flagged",
      TRUE ~ "ok"
    )
  )

qc_row_summary <- pheno_qc %>%
  count(qc_flag, name = "n_rows")

qc_traj_summary <- pheno_qc %>%
  distinct(traj_id, qc_flag) %>%
  count(qc_flag, name = "n_traj")

# ----------------------------
# Recalculate sparsity AFTER QC removal
# ----------------------------
traj_summary_qc <- pheno_qc %>%
  group_by(traj_id) %>%
  summarise(
    n_obs = n(),
    n_bbch = n_distinct(BBCH_adjusted),
    bbch_values = paste(sort(unique(BBCH_adjusted)), collapse = ", "),
    has_sowing = any(BBCH_adjusted == 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    traj_class = case_when(
      n_obs %in% c(1,2) ~ "too_sparse",
      n_obs %in% c( 3) ~ "sparse_keep_flagged",
      TRUE ~ "adequate"
    )
  )

bad_ids_after_qc <- traj_summary_qc %>%
  filter(traj_class == "too_sparse") %>%
  pull(traj_id)

bad_ids_after_qc <- unique(c(bad_ids_after_qc))

pheno_model_ready <- pheno_qc %>%
  filter(!traj_id %in% bad_ids_after_qc) 

pheno_model_ready <- pheno_model_ready %>%
  mutate(
    location = coalesce(location, Location, Location_name),
    location_id = interaction(location,lat,lon),
    location_id_original = as.character(location_id)
    ) 

pheno_model_ready <- pheno_model_ready %>%
  group_by(location, lat, lon) %>%
  mutate(Stations_id = cur_group_id()) %>%
  ungroup()


# ----------------------------
# Review tables
# ----------------------------
sparse_review <- traj_summary %>%
  filter(traj_class %in% c("too_sparse", "sparse_keep_flagged"))

sparse_review2 <- traj_summary_qc %>%
  filter(traj_class %in% c("too_sparse", "sparse_keep_flagged"))

manual_review <- bind_rows(
  manual_fix_log %>% mutate(review_source = "manual_fix"),
  manual_removed_review %>% mutate(review_source = "manual_remove")
)

# ----------------------------
# Exports
# ----------------------------
write_csv(pheno_model_ready, file.path(PEP_DIR, "pep_at_filtered_final_19_06.csv"))
write_csv(pheno_qc, file.path(PEP_DIR, "pep_at_phenology_qc_all.csv"))
write_csv(traj_summary, file.path(PEP_DIR, "pep_at_trajectory_summary_raw.csv"))
write_csv(traj_summary_qc, file.path(PEP_DIR, "pep_at_trajectory_summary_qc.csv"))
write_csv(qc_row_summary, file.path(PEP_DIR, "pep_at_qc_row_summary.csv"))
write_csv(qc_traj_summary, file.path(PEP_DIR, "pep_at_qc_traj_summary.csv"))
write_csv(manual_fix_log, file.path(PEP_DIR, "pep_at_manual_fix_log.csv"))
write_csv(candidate_decisions, file.path(PEP_DIR, "pep_at_year_offset_decisions.csv"))
write_csv(sparse_review, file.path(PEP_DIR, "pep_at_sparse_trajectories_review.csv"))
write_csv(sparse_review2, file.path(PEP_DIR, "pep_at_sparse_trajectories_after_qc_review.csv"))
write_csv(manual_review, file.path(PEP_DIR, "pep_at_manual_review_all.csv"))


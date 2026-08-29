## 01_download.R -- fetch the Movebank study, restricted to the analysis individuals.
##
## Inputs:  Movebank credentials in the OS keychain (see scripts/00_setup.R)
## Outputs: data/raw/movebank_deployments.rds  -- full deployment/reference table
##          data/raw/movebank_individuals.csv  -- selected individuals + derived cohort
##          data/raw/movebank_tracks.rds       -- GPS locations for selected individuals
##          figures/fig1_deployment_timeline.png -- the selection rule, made visible
##
## Usage:
##   Rscript scripts/01_download.R

suppressPackageStartupMessages({
  library(move2)
  library(dplyr)
  library(readr)
  library(lubridate)
})

source("config.R")
source(file.path("scripts", "utils.R"))
source(file.path("scripts", "plotting.R"))

ensure_dir(PATHS$raw_dir)

# ---------------------------------------------------------------------------------
# 1. Study metadata and license check
# ---------------------------------------------------------------------------------
# A Movebank account is required, and the license terms must be accepted once. If move2
# raises a license error, it includes the md5 hash of the terms; re-run passing that hash
# to accept them:
#
#   movebank_download_study(STUDY$movebank_id, 'license-md5' = '<hash from the error>')
#
# Accepting is a deliberate act, read the terms on movebank.org first (CC BY-NC).

message("Fetching study metadata ...")
study_info <- movebank_download_study_info(study_id = STUDY$movebank_id)
print(study_info)

if ("sensor_type_ids" %in% names(study_info)) {
  message("Sensors declared by the study: ",
          paste(unlist(study_info$sensor_type_ids), collapse = ", "))
}

# ---------------------------------------------------------------------------------
# 2. Deployment table
# ---------------------------------------------------------------------------------

message("Fetching deployment table ...")
deployments <- movebank_download_deployment(study_id = STUDY$movebank_id)
saveRDS(deployments, PATHS$deployments)
message("  ", nrow(deployments), " deployment rows.")

# Checking the presence of required columns (Movebank API field names).
required_cols <- c("individual_local_identifier", "deploy_on_timestamp",
                   "individual_comments", "tag_local_identifier")
missing_cols <- setdiff(required_cols, names(deployments))
if (length(missing_cols) > 0) {
  stop("Deployment table is missing: ", paste(missing_cols, collapse = ", "),
       "\nAvailable columns: ", paste(names(deployments), collapse = ", "), call. = FALSE)
}

# Flatten to the fields used here, dropping the sf geometry columns
# movebank_download_deployment() attaches (they break ordinary dplyr row operations).
# as.character on the identifier: move2 returns it as a factor, and the API filter below
# rejects anything that isn't character/logical/numeric. Ring is parsed from the
# identifier, not read from the inconsistent ring_id column (ring_from_animal_id()).
deploy_flat <- tibble(
  deployment_id = as.character(deployments$deployment_id),
  animal_id = as.character(deployments$individual_local_identifier),
  deploy_on = as.Date(deployments$deploy_on_timestamp),
  comment   = as.character(deployments$individual_comments),
  tag_id    = as.character(deployments$tag_local_identifier)
) |>
  filter(!is.na(animal_id), !is.na(deploy_on)) |>
  mutate(ring = ring_from_animal_id(animal_id),
         field_ring = field_ring_from_animal_id(animal_id))

# Hardcoded correction of deploy_on date for specified deployments.
for (r in names(SELECTION$deploy_on_corrections)) {
  new_date <- SELECTION$deploy_on_corrections[[r]]
  old_date <- deploy_flat$deploy_on[deploy_flat$ring == r]
  if (length(old_date) == 0) next
  message(sprintf("deploy_on correction: %s %s -> %s (config.R, SELECTION$deploy_on_corrections)",
                  r, paste(unique(old_date), collapse = "/"), new_date))
  deploy_flat$deploy_on[deploy_flat$ring == r] <- new_date
}

deploy_flat <- deploy_flat |>
  mutate(
    cohort      = cohort_from_date(deploy_on),
    is_standard = is_standard_protocol(deploy_on)
  )

# ---------------------------------------------------------------------------------
# 3. Select the analysis individuals
# ---------------------------------------------------------------------------------
# Standard protocol only: first-year crows trapped in the autumn/winter cage trap.
# Cohort is derived from the deployment date by rule. See config.R and utils.R.

individuals <- deploy_flat |>
  # A recaptured bird has several deployments; keep its first as the reference for
  # cohort and age, and carry the deployment count for the QC report.
  group_by(ring) |>
  arrange(deploy_on, .by_group = TRUE) |>
  summarize(
    animal_id     = first(animal_id),
    field_ring    = first(field_ring),
    first_deploy  = first(deploy_on),
    last_deploy   = last(deploy_on),
    n_deployments = n(),
    tag_id        = first(tag_id),
    n_tags        = n_distinct(tag_id),
    cohort        = first(cohort),
    is_standard   = first(is_standard),
    .groups = "drop"
  )

# Exclude individuals reported as rescued
rescued <- unique(deploy_flat$ring[grepl("rescued", deploy_flat$comment,
                                         ignore.case = TRUE)])
missed_by_rule <- intersect(rescued, individuals$ring[individuals$is_standard])
message(sprintf(
  "Rescued-and-released birds annotated in the metadata: %d (%d not already excluded by the date rule).",
  length(rescued), length(missed_by_rule)))
individuals <- individuals |> mutate(is_standard = is_standard & !(ring %in% rescued))

# Full cohort breakdown before filtering, so what the target-cohort restriction leaves
# behind is visible rather than implicit.
message("All cohorts present (standard protocol only):")
print(individuals |> filter(is_standard) |> count(cohort, name = "n_individuals"))

selected <- individuals |>
  filter(is_standard, cohort %in% SELECTION$cohorts_panel_a,
         !(ring %in% SELECTION$manual_exclusions))

message(sprintf("Individuals: %d total, %d standard-protocol in target cohorts.",
                nrow(individuals), nrow(selected)))
print(count(selected, cohort, name = "n_individuals"))

excluded <- individuals |> filter(!is_standard)
if (nrow(excluded) > 0) {
  message("Excluded as non-standard protocol (", nrow(excluded), "):")
  print(excluded |> select(ring, first_deploy, cohort))
}

# Manually excluded birds (config.R, SELECTION$manual_exclusions) -- either absent from
# Jiguet & Gantin (2025)'s own Table S1, or the sole 2023-cohort bird, outside this
# analysis's cohort scope. See config.R for the rationale behind each entry.
manually_excluded <- individuals |>
  filter(ring %in% SELECTION$manual_exclusions)
if (nrow(manually_excluded) > 0) {
  message("Manually excluded (", nrow(manually_excluded), "):")
  print(manually_excluded |> select(ring, first_deploy, cohort))
}

write_csv(selected, file.path(PATHS$raw_dir, "movebank_individuals.csv"))

# ---------------------------------------------------------------------------------
# 4. Figure -- Deployment timeline
# ---------------------------------------------------------------------------------
# One row per bird, one point per deployment, shaded bands for the trapping season;
# recaptures link as several points on one row. Color is per BIRD, not per deployment:
# selection is at the individual level, so a selected bird's later out-of-season
# redeployment still downloads and colors with the rest of its track.

years <- seq(year(min(deploy_flat$deploy_on)), year(max(deploy_flat$deploy_on)))

season_bands <- expand.grid(year = years, month = SELECTION$tagging_months) |>
  mutate(xmin = as.Date(sprintf("%d-%02d-01", year, month)),
         xmax = xmin + months(1))

year_starts <- as.Date(sprintf("%d-01-01", c(years, max(years) + 1)))
year_labels <- data.frame(mid = as.Date(sprintf("%d-07-01", years)), label = years)

status_levels <- c(as.character(sort(unique(selected$cohort))), "not selected")

# Any bird not in `selected` -- outside the trapping season, in a cohort this analysis
# does not cover, or individually excluded (SELECTION$manual_exclusions) -- reads the
# same way here: gray, "not selected". The reasons live in config.R, not in the figure.
timeline <- deploy_flat |>
  select(-cohort) |>
  left_join(select(individuals, ring, first_deploy, cohort), by = "ring") |>
  mutate(status = factor(
    if_else(ring %in% selected$ring, as.character(cohort), "not selected"),
    levels = status_levels
  )) |>
  arrange(first_deploy, ring) |>
  mutate(row = match(ring, unique(ring)))

p_timeline <- ggplot(timeline, aes(x = deploy_on, y = row)) +
  geom_rect(data = season_bands, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "gray93") +
  geom_vline(xintercept = year_starts, color = "gray55", linewidth = 0.35) +
  geom_line(aes(group = row), color = "gray70", linewidth = 0.3, lineend = "round") +
  geom_point(aes(color = status), size = 1.5) +
  scale_color_manual(
    name   = NULL,
    values = c(COHORT_COLORS, "not selected" = "gray60"),
    drop   = TRUE) +
  scale_x_date(breaks = year_labels$mid, labels = year_labels$label,
               minor_breaks = NULL, expand = expansion(mult = 0.01)) +
  scale_y_continuous(breaks = NULL) +
  labs(
    title    = "Figure 1: Deployment timeline and data selection",
    subtitle = sprintf("%d deployments on %d birds; %d selected",
                       nrow(deploy_flat), nrow(individuals), nrow(selected)),
    x = NULL, y = "individuals, ordered by first deployment",
    caption = paste(
      "One row per bird, one point per deployment, ordered by first deployment. Color is",
      "the cohort of the bird, taken from its first capture; gray marks birds not selected",
      "for the analysis. Vertical rules mark January 1st, shaded bands the autumn and",
      "winter trapping season.",
      sep = "\n")) +
  theme_corvus() +
  # The year boundaries are drawn as geoms.
  theme(panel.grid.major.x = element_blank()) +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 2.5)))

save_fig(p_timeline, "fig1_deployment_timeline", height_cm = 17)

# ---------------------------------------------------------------------------------
# 4. Download locations for the selected individuals
# ---------------------------------------------------------------------------------

message("Downloading GPS locations for ", nrow(selected), " individuals ",
        " ...")

tracks <- movebank_download_study(
  study_id                    = STUDY$movebank_id,
  sensor_type_id              = "gps",
  individual_local_identifier = selected$animal_id
)

if (!is.null(DATA_CUTOFF)) {
  before_n <- nrow(tracks)
  tracks <- tracks[as.Date(mt_time(tracks)) <= DATA_CUTOFF, ]
  message(sprintf("Cutoff at %s (config.R, DATA_CUTOFF): %d of %d locations kept.",
                  DATA_CUTOFF, nrow(tracks), before_n))
}

saveRDS(tracks, PATHS$raw_tracks)

message(sprintf("Saved %d locations for %d tracks to %s",
                nrow(tracks),
                length(unique(mt_track_id(tracks))),
                PATHS$raw_tracks))

# ---------------------------------------------------------------------------------
# 5. Sanity check: does deploy_on match when the tag actually started reporting?
# ---------------------------------------------------------------------------------
# Reports deployments where the tag starts reporting before deploy_on or more than 30 days after it.

first_fix_by_dep <- tibble(deployment_id = as.character(mt_track_id(tracks)),
                           timestamp = mt_time(tracks)) |>
  group_by(deployment_id) |>
  summarize(first_fix = min(timestamp), .groups = "drop")

fix_gap <- deploy_flat |>
  inner_join(first_fix_by_dep, by = "deployment_id") |>
  # deploy_on is a Date, first_fix a POSIXct; subtracting them directly silently produces
  # nonsense (a units mismatch that R does not error on). Put both on the same footing.
  mutate(gap_days = as.numeric(difftime(first_fix, as.POSIXct(deploy_on, tz = "UTC"),
                                        units = "days")))

suspect_late <- filter(fix_gap, gap_days > 30)
if (nrow(suspect_late) > 0) {
  message(sprintf(
    "\nWARNING: %d deployment(s) start transmitting far later than their recorded deploy_on date:",
    nrow(suspect_late)))
  print(as.data.frame(suspect_late |> select(ring, deployment_id, deploy_on, first_fix, gap_days)))
  message("  cohort_from_date() used deploy_on as recorded; verify manually before trusting",
          " this bird's cohort/age, and consider adding it to",
          " SELECTION$deploy_on_corrections.")
}

suspect_early <- filter(fix_gap, gap_days < 0)
if (nrow(suspect_early) > 0) {
  message(sprintf(
    "\nWARNING: %d deployment(s) start transmitting before their recorded deploy_on date:",
    nrow(suspect_early)))
  print(as.data.frame(suspect_early |> select(ring, deployment_id, deploy_on, first_fix, gap_days)))
  message("  cohort_from_date() used deploy_on as recorded; verify manually before trusting",
          " this bird's cohort/age, and consider adding it to",
          " SELECTION$deploy_on_corrections.")
}

if (nrow(suspect_late) == 0 && nrow(suspect_early) == 0) {
  message("\ndeploy_on vs first-fix check: no deployment outside the expected range.")
}

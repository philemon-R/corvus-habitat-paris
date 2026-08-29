## 02_clean_gps.R -- data quality diagnostic and track cleaning.
##
## Inputs:  data/raw/movebank_tracks.rds       (script 01)
##          data/raw/movebank_individuals.csv  (script 01)
## Outputs: data/interim/tracks_clean.rds      -- filtered locations, ring-keyed
##          data/processed/qc_summary.csv      -- per-individual diagnostic table
##          data/processed/individuals.csv     -- analysis cohort, post-filtering
##          figures/fig3_sampling_regime.png             -- realized fix interval, per bird per month
##          figures/fig4_tracking_effort_per_bird.png    -- fixes vs tracking span, and the threshold
##          figures/fig2_gps_cleaning.png                -- what each plausibility filter removed
##
## Usage:
##   Rscript scripts/02_clean_gps.R

suppressPackageStartupMessages({
  library(move2); library(dplyr); library(readr); library(lubridate)
  library(sf); library(units)
})

source("config.R")
source(file.path("scripts", "utils.R"))
source(file.path("scripts", "plotting.R"))

require_input(PATHS$raw_tracks, "scripts/01_download.R")
ensure_dir(PATHS$interim_dir)
ensure_dir(PATHS$processed_dir)

tracks <- readRDS(PATHS$raw_tracks)
n_start <- nrow(tracks)
message(sprintf("Loaded %d locations across %d tracks.",
                n_start, length(unique(mt_track_id(tracks)))))

# ---------------------------------------------------------------------------------
# 1. Re-key tracks by ring, not by deployment
# ---------------------------------------------------------------------------------
# move2 keys tracks on `deployment_id`. Birds recaptured for a tag replacement carry
# several deployments, so one individual arrives as several tracks. Collapse them onto
# the ring code, which is the stable individual identity.

track_data <- mt_track_data(tracks) |>
  mutate(ring = ring_from_animal_id(as.character(individual_local_identifier)))

n_before <- length(unique(mt_track_id(tracks)))
tracks <- tracks |>
  mt_set_track_data(track_data) |>
  mt_set_track_id("ring")

message(sprintf("Re-keyed by ring: %d deployment-tracks -> %d individuals.",
                n_before, length(unique(mt_track_id(tracks)))))

# Merging deployments can interleave records out of order; step building later assumes
# time ordering within a track.
if (!mt_is_time_ordered(tracks)) {
  tracks <- dplyr::arrange(tracks, mt_track_id(tracks), mt_time(tracks))
  message("  re-sorted records into time order")
}

# ---------------------------------------------------------------------------------
# 2. Duplicate timestamps
# ---------------------------------------------------------------------------------
# Two records at the same instant for the same bird make step length undefined. Keep the
# record with the better GPS quality where that is recorded, otherwise the first.

dup_key <- paste(as.character(mt_track_id(tracks)), as.character(mt_time(tracks)))
n_dup <- sum(duplicated(dup_key))
if (n_dup > 0) {
  hdop <- drop_units_if_any(st_drop_geometry(tracks)$gps_hdop)
  # Lower HDOP wins; NA sorts last, so a group with no quality data falls back to "first"
  # via the row-index tiebreak.
  priority <- ifelse(is.na(hdop), Inf, hdop)
  ord <- order(dup_key, priority, seq_along(dup_key))
  keep_rows <- sort(ord[!duplicated(dup_key[ord])])
  tracks <- tracks[keep_rows, ]
  message(sprintf("Removed %d duplicate (individual, timestamp) records, keeping the",
                  n_dup), " better-HDOP fix where quality is known.")
}

# ---------------------------------------------------------------------------------
# 3. Plausibility filters
# ---------------------------------------------------------------------------------
# Apply thresholds on speed and displacement as in Jiguet & Gantin (2025)

# Speed threshold.
bad_speed <- flag_both_sides(set_units(mt_speed(tracks), "km/h"),
                             QC$max_speed_kmh, mt_track_id(tracks))
message(sprintf("Speed filter (>%g km/h both sides): %d locations removed.",
                QC$max_speed_kmh, sum(bad_speed)))
tracks <- tracks[!bad_speed, ]

# Distance threshold.
bad_dist <- flag_both_sides(set_units(mt_distance(tracks), "m"),
                            QC$max_displacement_m, mt_track_id(tracks))
message(sprintf("Displacement filter (>%g km both sides): %d locations removed.",
                QC$max_displacement_m / 1000, sum(bad_dist)))
tracks <- tracks[!bad_dist, ]

# Filter GPS quality when available.
quality <- st_drop_geometry(tracks)
hdop <- drop_units_if_any(quality$gps_hdop)
sats <- drop_units_if_any(quality$gps_satellite_count)
has_quality <- !is.na(hdop)
poor_fix <- has_quality & (hdop > QC$max_hdop | sats < QC$min_satellites)
poor_fix[is.na(poor_fix)] <- FALSE

message(sprintf("GPS quality reported for %.1f%% of locations; %d flagged poor (HDOP>%g or sats<%d).",
                100 * mean(has_quality), sum(poor_fix), QC$max_hdop, QC$min_satellites))
if (QC$drop_poor_fixes) {
  tracks <- tracks[!poor_fix, ]
  message("  dropped (QC$drop_poor_fixes = TRUE)")
}

# Stationary tags: remove the tail of a track after the bird stopped moving.
#
# Coordinates are projected because the rule is in meters and the data are still in degrees at this point.
xy_m <- st_coordinates(st_transform(tracks, CRS$projected))
dead_tag <- flag_stationary_tail(xy_m[, 1], xy_m[, 2], mt_time(tracks), mt_track_id(tracks))

stationary_trimmed_rings <- unique(as.character(mt_track_id(tracks))[dead_tag])
n_dead_birds <- length(stationary_trimmed_rings)
message(sprintf("Stationary-tag tails: %d locations across %d birds.",
                sum(dead_tag), n_dead_birds))
if (sum(dead_tag) > 0) {
  tail_days <- tapply(as.Date(mt_time(tracks))[dead_tag],
                      as.character(mt_track_id(tracks))[dead_tag],
                      function(d) as.integer(max(d) - min(d)) + 1L)
  message("  longest tails (days): ",
          paste(sprintf("%s %d", names(sort(tail_days, decreasing = TRUE))[1:min(3, length(tail_days))],
                        sort(tail_days, decreasing = TRUE)[1:min(3, length(tail_days))]),
                collapse = ", "))
}
if (QC$drop_stationary_tail) {
  tracks <- tracks[!dead_tag, ]
  message("  dropped (QC$drop_stationary_tail = TRUE)")
}

# Store all filters effects
attrition <- tibble(
  step = factor(c("duplicate timestamps", "implausible speed",
                  "implausible displacement", "poor GPS quality",
                  "stationary tag"),
                levels = c("stationary tag", "poor GPS quality",
                           "implausible displacement", "implausible speed",
                           "duplicate timestamps")),
  n_removed = c(n_dup, sum(bad_speed), sum(bad_dist),
                if (QC$drop_poor_fixes) sum(poor_fix) else 0L,
                if (QC$drop_stationary_tail) sum(dead_tag) else 0L)
)

# ---------------------------------------------------------------------------------
# 4. Per-individual diagnostic
# ---------------------------------------------------------------------------------
# The realized fix interval varies with battery, season and logger model

individuals_in <- read_csv(file.path(PATHS$raw_dir, "movebank_individuals.csv"),
                           show_col_types = FALSE)

ev <- st_drop_geometry(tracks) |>
  mutate(ring = as.character(mt_track_id(tracks)),
         timestamp = mt_time(tracks))

qc <- ev |>
  arrange(ring, timestamp) |>
  group_by(ring) |>
  summarize(
    n_fixes        = n(),
    first_fix      = min(timestamp),
    last_fix       = max(timestamp),
    span_days      = as.numeric(difftime(max(timestamp), min(timestamp), units = "days")),
    days_with_data = n_distinct(as.Date(timestamp)),
    median_gap_min = median(as.numeric(diff(timestamp), units = "mins"), na.rm = TRUE),
    iqr_gap_min    = IQR(as.numeric(diff(timestamp), units = "mins"), na.rm = TRUE),
    max_gap_days   = max(as.numeric(diff(timestamp), units = "days"), na.rm = TRUE),
    hdop_ratio     = mean(!is.na(gps_hdop)),
    .groups = "drop"
  ) |>
  mutate(span_months = span_days / 30.44,
         # Share of calendar days in the tracking span that carry at least one fix.
         calendar_days = as.numeric(as.Date(last_fix) - as.Date(first_fix)) + 1,
         coverage = days_with_data / calendar_days) |>
  left_join(individuals_in |> select(ring, animal_id, field_ring, cohort, first_deploy,
                                     n_deployments, tag_id),
            by = "ring")

# GPS quality (HDOP/satellite count) is reported for all of a bird's fixes or none
# of them -- except for the two birds recaptured onto a different tag  mid-track
qc <- qc |>
  mutate(
    gps_quality = hdop_ratio > 0,
    stationary_trimmed = ring %in% stationary_trimmed_rings
  )

# ---------------------------------------------------------------------------------
# 5. Selection and derived analysis variables
# ---------------------------------------------------------------------------------

qc <- qc |>
  mutate(
    keep = n_fixes >= SELECTION$min_positions,
    independence_date = as.Date(paste0(cohort + 1, "-", AGE$independence_md)),
    # Panel B: needs enough tracking after independence to identify an age trend
    post_indep_months = as.numeric(difftime(last_fix, independence_date, units = "days")) / 30.44,
    panel_b = keep & post_indep_months >= SELECTION$min_post_indep_months_b
  )

message(sprintf("\nIndividuals: %d tracked, %d with >=%d fixes, %d in the longitudinal panel.",
                nrow(qc), sum(qc$keep), SELECTION$min_positions, sum(qc$panel_b)))
print(qc |> count(cohort, keep, name = "n"))

kept_rings <- qc$ring[qc$keep]
tracks <- tracks[as.character(mt_track_id(tracks)) %in% kept_rings, ]

# indep_years: a continuous variable relative to each bird's own independence date.
# Negative before independence, 0 exactly at independence, positive after.
track_data <- mt_track_data(tracks) |>
  left_join(qc |> select(ring, independence_date), by = "ring")

tracks <- mt_set_track_data(tracks, track_data)

ev_time <- mt_time(tracks)
ring_of <- as.character(mt_track_id(tracks))
lookup  <- track_data |> select(ring, independence_date)
idx     <- match(ring_of, lookup$ring)

tracks$indep_years <- as.numeric(difftime(
  ev_time, as.POSIXct(lookup$independence_date[idx], tz = "UTC"), units = "days")) / 365.25

message(sprintf("Years since independence: median %.2f, range %.2f to %.2f",
                median(tracks$indep_years), min(tracks$indep_years), max(tracks$indep_years)))
message(sprintf("Fixes before independence: %s (%.1f%%)",
                format(sum(tracks$indep_years < 0), big.mark = " "),
                100 * mean(tracks$indep_years < 0)))

# ---------------------------------------------------------------------------------
# 6. Figures
# ---------------------------------------------------------------------------------
# Three diagnostics on data quality, selection and cleaning.

## 6a. Sampling regime -----------------------------------------------------------------
# The realized fix interval, per bird per month. Gaps in the tiles are months with no data at all.

regime <- ev |>
  filter(ring %in% kept_rings) |>
  arrange(ring, timestamp) |>
  group_by(ring) |>
  mutate(gap_min = as.numeric(difftime(timestamp, lag(timestamp), units = "mins"))) |>
  ungroup() |>
  mutate(month = as.Date(format(timestamp, "%Y-%m-01"))) |>
  group_by(ring, month) |>
  summarize(median_gap = median(gap_min, na.rm = TRUE), .groups = "drop") |>
  left_join(qc |> select(ring, cohort, first_fix, gps_quality, stationary_trimmed),
            by = "ring") |>
  mutate(gap_class = cut(median_gap,
                         breaks = c(0, 15, 45, 90, 360, Inf),
                         labels = c("under 15 min", "15-45 min", "45-90 min",
                                    "1.5-6 h", "over 6 h")))

# Order birds within their cohort by when tracking started.
regime <- regime |>
  group_by(cohort) |>
  mutate(row = dense_rank(first_fix)) |>
  ungroup()

# Per-bird properties (GPS quality reported, stationary tag trimmed) are shown as
# markers in the margin.
margin_x <- min(regime$month) - c("stationary tag trimmed" = 145,
                                  "GPS quality reported"   = 100)

margin_marks <- bind_rows(
  regime |> filter(gps_quality) |> distinct(cohort, row) |>
    mutate(feature = "GPS quality reported"),
  regime |> filter(stationary_trimmed) |> distinct(cohort, row) |>
    mutate(feature = "stationary tag trimmed")
) |>
  mutate(month = margin_x[feature],
         feature = factor(feature, levels = names(margin_x)))

# Deployment frontiers computed for birds tagged several times.
deployment_frontiers <- ev |>
  filter(ring %in% kept_rings) |>
  group_by(ring, deployment_id) |>
  summarize(dep_start = min(timestamp), .groups = "drop") |>
  arrange(ring, dep_start) |>
  group_by(ring) |>
  filter(n() > 1) |>
  slice(-1) |>
  ungroup() |>
  # Snapped to the start of its month, matching the tiles' boundaries.
  mutate(dep_month = as.Date(format(dep_start, "%Y-%m-01"))) |>
  left_join(distinct(regime, ring, cohort, row), by = "ring")

message(sprintf("Deployment frontiers marked on the regime figure: %d bird(s).",
                n_distinct(deployment_frontiers$ring)))

# geom_rect with true month bounds, not geom_tile: a fixed-width tile centered on the 1st
# puts the 1 January gridline through the January tile instead of on its boundary.
regime <- regime |> mutate(month_end = month %m+% months(1))

p_regime <- ggplot(regime, aes(y = row)) +
  geom_rect(aes(xmin = month, xmax = month_end, ymin = row - 0.45, ymax = row + 0.45,
               fill = gap_class)) +
  geom_segment(data = deployment_frontiers,
               aes(x = dep_month, xend = dep_month, y = row - 0.45, yend = row + 0.45),
               inherit.aes = FALSE, color = "red", linewidth = 0.6) +
  geom_point(data = margin_marks, aes(x = month, y = row, shape = feature),
             color = "gray25", size = 1.3) +
  facet_grid(cohort ~ ., scales = "free_y", space = "free_y") +
  scale_fill_viridis_d(name = "median interval between GPS fixes", direction = -1,
                       na.value = "gray85") +

  scale_shape_manual(name = NULL, values = c("GPS quality reported"   = 18,
                                             "stationary tag trimmed" = 4)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_reverse(breaks = NULL) +
  labs(
    title    = "Figure 3: Sampling regime",
    subtitle = sprintf("%d individuals, one row each, one tile per month tracked",
                       length(kept_rings)),
    x = NULL, y = "individuals, by cohort, earliest first",
    caption = paste(
      "Median interval between consecutive locations, within each bird and month. A blank",
      "month is one with no location at all. A red bar marks a mid-track tag change;",
      "for the two birds this applies to, GPS quality is reported only from the second",
      "deployment onward.",
      sep = "\n")) +
  theme_corvus() +
  theme(legend.box = "vertical") +
  guides(fill = guide_legend(nrow = 1, title.position = "top", order = 1),
         shape = guide_legend(order = 2))

save_fig(p_regime, "fig3_sampling_regime", height_cm = 19)

## 6b. Tracking effort -----------------------------------------------------------------
# Plots fixes against tracking span.

p_effort <- ggplot(qc, aes(x = span_months, y = n_fixes)) +
  geom_hline(yintercept = SELECTION$min_positions, linetype = "dashed",
             color = "gray40") +
  annotate("text", x = Inf, y = SELECTION$min_positions, hjust = 1.05, vjust = -0.6,
           label = sprintf("inclusion threshold: %d locations", SELECTION$min_positions),
           color = "gray40", size = 3) +
  geom_point(aes(color = cohort_factor(cohort), shape = panel_b),
             size = 2.2, stroke = 0.8) +
  scale_color_manual(values = COHORT_COLORS, name = "Cohort", drop = TRUE) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19), name = NULL,
                     labels = c(`FALSE` = "not in the longitudinal panel",
                                `TRUE` = "longitudinal panel")) +
  scale_y_log10(labels = scales::label_comma()) +
  labs(
    title    = "Figure 4: Tracking effort per bird",
    subtitle = sprintf("%d of %d selected birds returned data; %d above the threshold, %d in the longitudinal panel",
                       nrow(qc), nrow(individuals_in), sum(qc$keep), sum(qc$panel_b)),
    x = "tracking span (months)", y = "locations retained after cleaning (log scale)",
    caption = sprintf(
      "Filled points: birds tracked at least %d months past their own independence date.",
      SELECTION$min_post_indep_months_b)) +
  theme_corvus()

save_fig(p_effort, "fig4_tracking_effort_per_bird", height_cm = 15)

## 6c. Cleaning attrition --------------------------------------------------------------
# Effect of each filter shown as horizontal bar graph.

p_attrition <- ggplot(attrition, aes(x = n_removed, y = step)) +
  geom_col(fill = "#0072B2", width = 0.6) +
  geom_text(aes(label = format(n_removed, big.mark = " ")),
            hjust = -0.15, size = 3.2, color = "gray25") +
  scale_x_continuous(labels = scales::label_comma(),
                     expand = expansion(mult = c(0, 0.18))) +
  labs(
    title    = "Figure 2: Cleaning GPS tracks",
    subtitle = sprintf("%s of %s locations removed in total (%.2f%%)",
                       format(sum(attrition$n_removed), big.mark = " "),
                       format(n_start, big.mark = " "),
                       100 * sum(attrition$n_removed) / n_start),
    x = "locations removed", y = NULL,
    caption = "Filters run top to bottom, each on the output of the last.") +
  theme_corvus()

save_fig(p_attrition, "fig2_gps_cleaning", height_cm = 11)

# ---------------------------------------------------------------------------------
# 7. Project and write
# ---------------------------------------------------------------------------------
# Everything downstream (step lengths, buffers, OSM overlays) is metric, so project once
# here rather than repeatedly later.

tracks <- st_transform(tracks, CRS$projected)

saveRDS(tracks, PATHS$clean_tracks)
write_csv(qc, PATHS$qc_summary)
write_csv(qc |> filter(keep), PATHS$individuals)

message(sprintf("\nWrote %d locations for %d individuals to %s",
                nrow(tracks), length(unique(mt_track_id(tracks))), PATHS$clean_tracks))
message("QC summary: ", PATHS$qc_summary)

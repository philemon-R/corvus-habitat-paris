## 05_settled_nights.R -- compute one settled position per bird per night.
##
## Inputs:  data/interim/tracks_clean.rds     (script 02)
## Outputs: data/interim/settled_nights.rds   -- one row per bird-night: median position,
##                                                fix count, settled-share, settled flag
##
## Usage:
##   Rscript scripts/05_settled_nights.R

suppressPackageStartupMessages({
  library(move2); library(dplyr); library(readr); library(sf); library(suncalc)
})

source("config.R")
source(file.path("scripts", "utils.R"))

require_input(PATHS$clean_tracks, "scripts/02_clean_gps.R")
ensure_dir(PATHS$interim_dir)

tracks <- readRDS(PATHS$clean_tracks)
xy <- st_coordinates(tracks)
pts <- tibble(
  ring      = as.character(mt_track_id(tracks)),
  timestamp = mt_time(tracks),
  x         = xy[, "X"],
  y         = xy[, "Y"],
  date      = as.Date(mt_time(tracks))
)
message(sprintf("Loaded %s locations for %d individuals.",
                format(nrow(pts), big.mark = " "), n_distinct(pts$ring)))

# ---------------------------------------------------------------------------------
# 1. Which fixes are nocturnal, and which night they belong to
# ---------------------------------------------------------------------------------
# night_bounds() (utils.R), shared with scripts 06/07: drops the hours containing sunrise/sunset.
# A night is labeled by the date it started on, so 23:00 and the following 01:00 are the same night.

sun <- night_bounds(pts$date)

pts <- pts |>
  left_join(sun, by = "date") |>
  mutate(is_night = timestamp >= night_starts | timestamp < night_ends,
         night    = if_else(timestamp < night_ends, date - 1, date))

message(sprintf("Nocturnal fixes: %s (%.1f%% of locations).",
                format(sum(pts$is_night), big.mark = " "), 100 * mean(pts$is_night)))

# ---------------------------------------------------------------------------------
# 2. One settled position per individual-night
# ---------------------------------------------------------------------------------
# The median of a night's fixes.
#
# A night is kept only if the bird was demonstrably settled: more than half its fixes
# within one site radius of that median.

per_night <- pts |>
  filter(is_night) |>
  group_by(ring, night) |>
  mutate(mx = median(x), my = median(y)) |>
  summarize(n_fixes = n(),
            settled_share = mean(sqrt((x - mx)^2 + (y - my)^2) <= ROOST$site_radius_m),
            x = first(mx), y = first(my),
            .groups = "drop") |>
  mutate(settled = settled_share > ROOST$settled_min_share)

message(sprintf("Bird-nights: %s, of which %.1f%% settled within %g m of their own median.",
                format(nrow(per_night), big.mark = " "),
                100 * mean(per_night$settled), ROOST$site_radius_m))
message(sprintf("  fixes per night: median %d, a quarter of nights have %d or fewer, %.1f%% have one",
                median(per_night$n_fixes), quantile(per_night$n_fixes, 0.25),
                100 * mean(per_night$n_fixes == 1)))

saveRDS(per_night, PATHS$settled_nights)
message(sprintf("\nWrote %s bird-night positions to %s",
                format(nrow(per_night), big.mark = " "), PATHS$settled_nights))

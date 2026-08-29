## 06_maps.R -- Visualizations of GPS-tracks
##
## Every map is drawn with fused_base_map() (plotting.R): OCS GE's three covariate classes
## for land cover, OSM buildings/water on top.
## Buildings are individual footprints on trajectory panels up to
## MAP$building_context_max_halfwidth_m, a density overlay beyond it..
#
## Inputs:  data/interim/tracks_clean.rds    (script 02)
##          data/interim/settled_nights.rds  (script 05)
##          data/osm/*.gpkg                  (script 03)
## Outputs: figures/fig5_day_night_positions.png       -- distinct birds per cell, day against night
##          figures/fig6_night_positions_over_year.png -- individual tracks over one calendar year
##          figures/fig7_tracks_over_month.png         -- individual tracks over one dense month
##          figures/fig8_long_dist_examples.png        -- two long-distance movement profiles
##          data/processed/spatial_coverage.csv -- per-bird share inside the study area
##          data/osm/building_density_25m.tif    -- built-up density (built on first run)
##
## Usage:
##   Rscript scripts/06_maps.R

suppressPackageStartupMessages({
  library(move2); library(dplyr); library(readr); library(sf)
  library(suncalc); library(ggplot2); library(patchwork); library(ggspatial)
})

source("config.R")
source(file.path("scripts", "utils.R"))
source(file.path("scripts", "osm.R"))
source(file.path("scripts", "ocsge.R"))
source(file.path("scripts", "plotting.R"))

require_input(PATHS$clean_tracks, "scripts/02_clean_gps.R")
require_input(PATHS$settled_nights, "scripts/05_settled_nights.R")
use_osm_cache()
ensure_dir(PATHS$figures_dir)

tracks <- readRDS(PATHS$clean_tracks)
message(sprintf("Loaded %d locations for %d individuals.",
                nrow(tracks), length(unique(mt_track_id(tracks)))))

# ---------------------------------------------------------------------------------
# 1. Flatten to a plain coordinate table
# ---------------------------------------------------------------------------------

xy <- st_coordinates(tracks)
pts <- tibble(
  ring      = as.character(mt_track_id(tracks)),
  timestamp = mt_time(tracks),
  x         = xy[, "X"],
  y         = xy[, "Y"]
)

# Extract time components.
lt <- as.POSIXlt(pts$timestamp, tz = "UTC")
pts$date  <- as.Date(pts$timestamp)
pts$year  <- lt$year + 1900L
pts$yday  <- lt$yday + 1L
pts$clock <- lt$hour + lt$min / 60    # hour of day as a fraction, for the daily figure

# First of each location's month, built by matching on the ~80 distinct year/month pairs
# rather than per row.
ym <- lt$year * 12L + lt$mon
months <- sort(unique(ym))
pts$month <- as.Date(sprintf("%d-%02d-01",
                             months %/% 12L + 1900L, months %% 12L + 1L))[match(ym, months)]
rm(lt, ym, months)

trap <- trap_point()

# ---------------------------------------------------------------------------------
# 2. Day and night
# ---------------------------------------------------------------------------------

sun <- night_bounds(pts$date)

pts <- pts |>
  left_join(sun, by = "date") |>
  mutate(is_night = timestamp >= night_starts | timestamp < night_ends,
         period   = factor(ifelse(is_night, "night", "day"), levels = c("day", "night")))

message(sprintf("Day/night split: %.1f%% of locations are nocturnal.",
                100 * mean(pts$is_night)))

# ---------------------------------------------------------------------------------
# 3. Spatial coverage against the habitat study area
# ---------------------------------------------------------------------------------

bbox_proj <- study_boundary() |> st_transform(CRS$projected) |> st_bbox()
pts <- pts |>
  mutate(inside = x >= bbox_proj[["xmin"]] & x <= bbox_proj[["xmax"]] &
                  y >= bbox_proj[["ymin"]] & y <= bbox_proj[["ymax"]],
         dist_trap_km = sqrt((x - trap$x)^2 + (y - trap$y)^2) / 1000)

coverage <- pts |>
  group_by(ring) |>
  summarize(n_fixes         = n(),
            share_inside    = mean(inside),
            median_dist_km  = median(dist_trap_km),
            max_dist_km     = max(dist_trap_km),
            .groups = "drop") |>
  arrange(share_inside)

write_csv(coverage, file.path(PATHS$processed_dir, "spatial_coverage.csv"))

message(sprintf("Inside the habitat study area: %.2f%% of locations, %d of %d birds entirely.",
                100 * mean(pts$inside),
                sum(coverage$share_inside == 1), nrow(coverage)))
leavers <- coverage |> filter(share_inside < 0.5)
if (nrow(leavers) > 0) {
  message(sprintf("  %d birds spend most of their tracked life outside it:", nrow(leavers)))
  print(as.data.frame(leavers))
}

# ---------------------------------------------------------------------------------
# 3.1. Prepare second calendar year
# ---------------------------------------------------------------------------------
# Part fixes according to calendar years

cohorts <- read_csv(PATHS$qc_summary, show_col_types = FALSE) |>
  filter(keep) |>
  transmute(ring, field_ring, cohort,
            win_start = pmax(as.Date(first_fix), as.Date(sprintf("%d-01-01", cohort + 1))),
            win_end   = as.Date(sprintf("%d-12-31", cohort + 1)),
            truncated = win_start > as.Date(sprintf("%d-01-01", cohort + 1)))

# Field-ring codes for labeling panels.
field_ring_by_ring <- setNames(cohorts$field_ring, cohorts$ring)

# The calendar year a location belongs to, counted from the hatching year: a crow hatched
# in year Y is in its first calendar year during Y, its second during Y+1, and so on.
pts <- pts |>
  left_join(select(cohorts, ring, cohort, truncated), by = "ring") |>
  mutate(calendar_year = year - cohort + 1L)

scy   <- filter(pts, calendar_year == 2)   # the second calendar year
later <- filter(pts, calendar_year >= 3)   # the third and everything after it

for (part in list(list("first calendar year", filter(pts, calendar_year <= 1)),
                  list("second calendar year", scy),
                  list("third and later", later))) {
  message(sprintf("  %-22s %10s locations, %2d birds",
                  part[[1]], format(nrow(part[[2]]), big.mark = " "),
                  n_distinct(part[[2]]$ring)))
}
message(sprintf("  %d birds were trapped inside their second calendar year, so it opens late.",
                sum(cohorts$truncated)))

# ---------------------------------------------------------------------------------
# 3.2. Hourly thinning
# ---------------------------------------------------------------------------------
# Thinning to one location per bird per clock hour allows to average tracks with different
# sampling regimes.

hourly <- thin_hourly(pts, ring, timestamp)
message(sprintf("Hourly thinning: %s of %s locations kept (%.0f%%).",
                format(nrow(hourly), big.mark = " "), format(nrow(pts), big.mark = " "),
                100 * nrow(hourly) / nrow(pts)))

# ---------------------------------------------------------------------------------
# 4. Figure: day against night
# ---------------------------------------------------------------------------------
# Distinct birds per cell, from the hourly-thinned set.

# Dense core only: half the locations sit within 2.4 km of the trap. The dispersal figure
# further down picks up the tail.
core <- square_box(trap$x, trap$y, MAP$core_halfwidth_m)

in_core <- hourly |>
  filter(x >= core[["xmin"]], x <= core[["xmax"]],
         y >= core[["ymin"]], y <= core[["ymax"]])

day_night <- in_core |>
  group_by(period, gx = round(x / MAP$cell_m), gy = round(y / MAP$cell_m)) |>
  summarize(value = n_distinct(ring), .groups = "drop") |>
  mutate(x = gx * MAP$cell_m, y = gy * MAP$cell_m)

p_daynight <- ggplot() +
  fused_base_map(core) +
  geom_tile(data = day_night, aes(x = x, y = y, fill = value),
            width = MAP$cell_m, height = MAP$cell_m, alpha = 0.6) +
  facet_wrap(~ period, nrow = 1) +
  scale_fill_viridis_c(trans = "log10", option = "magma", direction = -1,
                       name = "distinct birds",
                       guide = guide_colorbar(draw.ulim = FALSE, draw.llim = FALSE)) +
  labs(x = NULL, y = NULL)
p_daynight <- frame_map(p_daynight, core)

p_daynight <- (p_daynight / fused_legend()) +
  plot_layout(heights = c(1, 0.09)) +
  plot_annotation(
    title    = "Figure 5: Day and night positions",
    subtitle = sprintf("%s locations from %d birds, %g m cells, shared log color scale",
                       format(nrow(in_core), big.mark = " "),
                       n_distinct(in_core$ring), MAP$cell_m),
    caption  = paste(
      "At most one location per bird per hour, so a fast-sampling tag does not reach more cells",
      "than a slow one. Night runs from the first whole hour after sunset to the last whole hour",
      sprintf("before sunrise, so twilight counts as day; %.0f%% of locations are nocturnal.",
              100 * mean(hourly$is_night)),
      "The cross marks the trapping site.",
      sep = "\n"),
    theme = theme_corvus())
save_fig(p_daynight, "fig5_day_night_positions", height_cm = 18)

# ---------------------------------------------------------------------------------
# 5. Figure: individual trajectories
# ---------------------------------------------------------------------------------

pick_spread <- function(candidates, n) {
  candidates <- arrange(candidates, max_km)
  candidates$ring[round(seq(1, nrow(candidates), length.out = n))]
}

#' Half-width of a bird's panel frame. Used both to pick eligible birds and to draw them,
#' so the two never disagree.
panel_half <- function(x, y) {
  max(MAP$track_min_halfwidth_m, 0.55 * max(diff(range(x)), diff(range(y))))
}

#' Whether a panel frame lies entirely inside the map context's area. Centered on the
#' bird.
study_box <- st_bbox(st_transform(study_boundary(), CRS$projected))

frame_inside_study <- function(cx, cy, half) {
  cx - half >= study_box[["xmin"]] & cx + half <= study_box[["xmax"]] &
  cy - half >= study_box[["ymin"]] & cy + half <= study_box[["ymax"]]
}

#' One trajectory panel, framed on its own data.
track_panel <- function(d, title, color_scale, subtitle = NULL, as_path = TRUE) {
  half <- panel_half(d$x, d$y)
  cx <- mean(range(d$x)); cy <- mean(range(d$y))
  box <- square_box(cx, cy, half)

  ggplot() +
    fused_base_map(box) +
    # Low alpha, not a solid line; lineend = "round" to handle segment joints.
    (if (as_path) {
       geom_path(data = d, aes(x = x, y = y, color = t), linewidth = 0.5, alpha = 0.3,
                 lineend = "round", linejoin = "round")
     } else {
       # Points at the same alpha, shape = 16.
       geom_point(data = d, aes(x = x, y = y, color = t), size = 1.0, alpha = 0.5, shape = 16)
     }) +
    geom_point(data = trap, aes(x = x, y = y), shape = 4, size = 2, stroke = 0.9,
               color = "gray15") +
    color_scale +
    ggspatial::annotation_scale(location = "br", line_width = 0.5, height = unit(0.1, "cm"),
                                text_cex = 0.55, text_col = "gray30",
                                bar_cols = c("gray30", "white"), pad_x = unit(0.15, "cm"),
                                pad_y = unit(0.15, "cm")) +
    coord_sf(xlim = c(box[["xmin"]], box[["xmax"]]),
             ylim = c(box[["ymin"]], box[["ymax"]]),
             expand = FALSE, crs = st_crs(CRS$projected), datum = NA) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_corvus() +
    theme(axis.text = element_blank(), panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(size = rel(0.9)),
          plot.subtitle = element_text(size = rel(0.9)))
}

# ---------------------------------------------------------------------------------
# 5.1. Figure: nightly positions over one calendar year
# ---------------------------------------------------------------------------------
# One position per bird per night: the median of that night's fixes, only on nights the
# bird was demonstrably settled (script 05).

nightly <- readRDS(PATHS$settled_nights) |>
  filter(settled) |>
  left_join(select(cohorts, ring, cohort), by = "ring") |>
  mutate(calendar_year = as.integer(format(night, "%Y")) - cohort + 1L) |>
  filter(calendar_year == 2) |>
  mutate(t = as.integer(format(night, "%j")))

seasonal_candidates <- nightly |>
  group_by(ring) |>
  summarize(nights = n(),
            max_km = max(sqrt((x - trap$x)^2 + (y - trap$y)^2)) / 1000,
            # For the subtitle, not selection (still max_km above): panels are framed on
            # the bird's own range, and the trap is off-frame in several of them.
            spread_km = if (n() > 1) max(dist(cbind(x, y))) / 1000 else 0,
            half   = panel_half(x, y),
            cx     = mean(range(x)), cy = mean(range(y)),
            .groups = "drop") |>
  # Eligible if the year is well covered and the frame fits inside the study area.
  filter(nights >= MAP$track_min_days * 0.6,
         frame_inside_study(cx, cy, half))

message(sprintf("Seasonal panels: %d birds eligible (frame inside the study area), showing %d.",
                nrow(seasonal_candidates), MAP$track_n_seasonal))

seasonal_rings <- pick_spread(seasonal_candidates, MAP$track_n_seasonal)

# Day 365 is 31 December, not 1 January: labeling both ends "Jan" implied the scale
# wrapped, which a calendar year does not.
seasonal_scale <- scale_season(
  name = NULL, limits = c(1, 365),
  breaks = c(1, 91, 182, 274, 365), labels = c("Jan", "Apr", "Jul", "Oct", "Dec"))

# One shared legend, reused here and by figure 8, extracted from a throwaway reference
# plot rather than collected across panels with guides = "collect": collecting across
# independently-scaled panels forced patchwork into a much wider left margin, traced to
# that reconciliation specifically.
extract_legend <- function(p) {
  g <- ggplotGrob(p)
  g$grobs[[which(sapply(g$grobs, `[[`, "name") == "guide-box")]]
}
season_legend <- extract_legend(
  ggplot(data.frame(x = 1, y = 1, t = c(1, 365)), aes(x, y, color = t)) +
    geom_point() + seasonal_scale + theme_corvus() + theme(legend.position = "bottom"))

# Boxed in spacers so the bar keeps its native width rather than stretching full-row (same
# trick clock_strip uses below); the 1:2:1 ratio is a first guess, not computed.
season_strip <- (plot_spacer() | wrap_elements(full = season_legend) | plot_spacer()) +
  plot_layout(widths = c(1, 2, 1))

# Points, not a joined path: two consecutive nights aren't a journey (the bird moved all
# day in between), and a straight line between roosts would draw a route never flown.
p_seasonal <- wrap_plots(lapply(seasonal_rings, function(r) {
  d <- filter(nightly, ring == r) |> arrange(night)
  info <- filter(seasonal_candidates, ring == r)
  track_panel(d, field_ring_by_ring[[r]], list(seasonal_scale, guides(color = "none")),
              as_path = FALSE,
              subtitle = sprintf("%d nights, up to %.1f km apart", info$nights, info$spread_km))
}), nrow = 2)

p_seasonal <- (p_seasonal / season_strip) +
  plot_layout(heights = c(1, 0.11)) +
  plot_annotation(
    title    = "Figure 6: Nightly positions over one calendar year",
    subtitle = "One dot per bird per night, across each bird's second calendar year",
    caption  = paste(
      "Median of a bird's locations between the first whole hour after sunset and the last",
      "whole hour before sunrise, on nights it was demonstrably settled there.",
      sprintf("%d birds ranging up to %.0f km from the trap, sampled at even quantiles of that range.",
              MAP$track_n_seasonal,
              max(filter(seasonal_candidates, ring %in% seasonal_rings)$max_km)),
      "The cross marks the trapping site.",
      sep = "\n"),
    theme = theme_corvus())

save_fig(p_seasonal, "fig6_night_positions_over_year", height_cm = 21)

# ---------------------------------------------------------------------------------
# 5.2. Figure 7: individual tracks over one month
# ---------------------------------------------------------------------------------
# Hour of day needs a cyclic color scale (midnight next to both 23:00 and 01:00);
# CLOCK_COLORS (plotting.R) wraps for this.
#
# Every fix, not nightly means, over one month: a commute is a daily event, visible only
# where sampling is dense. The month shown per bird is its best-sampled one in its second
# calendar year.
best_month <- scy |>
  count(ring, month, name = "n") |>
  group_by(ring) |>
  slice_max(n, n = 1, with_ties = FALSE) |>
  ungroup()

daily <- scy |>
  inner_join(best_month, by = c("ring", "month")) |>
  mutate(t = clock)

daily_candidates <- daily |>
  group_by(ring) |>
  summarize(n = n(), month = first(month),
            max_km = max(sqrt((x - trap$x)^2 + (y - trap$y)^2)) / 1000,
            half   = panel_half(x, y),
            cx     = mean(range(x)), cy = mean(range(y)),
            .groups = "drop") |>
  filter(
    n >= 1500,   # sub-hourly sampling needed to show a rhythm at all -- a logger selection
    frame_inside_study(cx, cy, half),
    # Tighter than the seasonal figure: a commute needs the same zoom that keeps building
    # footprints legible.
    half <= MAP$building_context_max_halfwidth_m)

message(sprintf("Daily-rhythm panels: %d birds eligible, showing %d.",
                nrow(daily_candidates), MAP$track_n_daily))

daily_rings <- pick_spread(daily_candidates, MAP$track_n_daily)

daily_scale <- scale_color_gradientn(colors = CLOCK_COLORS, limits = c(0, 24),
                                      guide = "none")

# Drawn as a ring, not a bar: a bar has to cut the cyclic scale somewhere, and the cut
# would read as a boundary that isn't there. Built as its own plot, since ggplot has no
# cyclic guide.
clock_legend <- function() {
  band  <- data.frame(hour = seq(0, 24, length.out = 145)) |>
    mutate(xmin = hour, xmax = hour + 24 / 144)
  # Anchored, not centered: a polar plot draws labels horizontally, so a centered label on
  # the ring's left would grow back across the ticks.
  major <- data.frame(hour  = c(0, 6, 12, 18),
                      label = c("00:00", "06:00", "12:00", "18:00"),
                      hjust = c(0.5, 0,   0.5, 1),
                      vjust = c(0,   0.5, 1,   0.5))
  minor <- data.frame(hour = setdiff(seq(0, 23), major$hour))

  ggplot() +
    geom_rect(data = band,
              aes(xmin = xmin, xmax = xmax, ymin = 1.30, ymax = 2.00, fill = hour)) +
    # Ticks outside the band, longer where labeled, so the ring reads like an axis.
    geom_segment(data = minor, aes(x = hour, xend = hour, y = 2.06, yend = 2.20),
                 color = "gray55", linewidth = 0.25, lineend = "round") +
    geom_segment(data = major, aes(x = hour, xend = hour, y = 2.06, yend = 2.34),
                 color = "gray30", linewidth = 0.45, lineend = "round") +
    geom_text(data = major, aes(x = hour, y = 2.44, label = label,
                                hjust = hjust, vjust = vjust),
              size = 2.6, color = "gray25") +
    scale_fill_gradientn(colors = CLOCK_COLORS, limits = c(0, 24), guide = "none") +
    coord_polar(theta = "x", start = 0) +
    scale_x_continuous(limits = c(0, 24), breaks = NULL) +
    scale_y_continuous(limits = c(0, 3.1), breaks = NULL, expand = expansion(mult = 0)) +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

p_daily <- wrap_plots(lapply(daily_rings, function(r) {
  d <- filter(daily, ring == r) |> arrange(timestamp)
  info <- filter(daily_candidates, ring == r)
  track_panel(d, field_ring_by_ring[[r]], daily_scale,
              sprintf("%s %s, %s locations",
                      month.name[as.integer(format(info$month, "%m"))],
                      format(info$month, "%Y"),
                      format(info$n, big.mark = " ")))
}), nrow = 2)

# The ring sits in its own strip: spacers bound its WIDTH (a polar panel expands to fill
# whatever width it's given, then claims a matching height), and the height ratio is set
# to match, or surplus row height sits blank above/below the forced-square ring.
clock_strip <- (plot_spacer() | clock_legend() | plot_spacer()) +
  plot_layout(widths = c(2, 1.4, 2))

p_daily <- (p_daily / clock_strip) +
  plot_layout(heights = c(1, 0.30)) +
  plot_annotation(
    title    = "Figure 7: Individual tracks over one month",
    subtitle = "Every location, joined in time order, colored by time of day",
    caption  = paste(
      "Each bird's best-sampled month inside its second calendar year, among birds with at least",
      "1 500 locations in that month and a range compact enough for a commute to be visible.",
      "Every location is drawn, at the logger's own rate. Lines are translucent, so a route flown",
      "daily reads darker than one flown once. The cross marks the trapping site.",
      sep = "\n"),
    theme = theme_corvus())

save_fig(p_daily, "fig7_tracks_over_month", height_cm = 28)

# ---------------------------------------------------------------------------------
# 6. Figure: long-distance movement, two profiles
# ---------------------------------------------------------------------------------
# The two figures above stay inside the habitat study area by construction, so neither
# shows the birds that range beyond 40 km from the trap (up to 154 km). No OSM context
# exists that far out; maps::france (bundled, no download) gives a coarse alternative --
# the 96 metropolitan departments (not the old pre-2016 regions, though boundaries are
# unaffected by that merger either way), kept undissolved since a single national
# silhouette draws as a flat tint at this panel scale (40-80 km across, same order as a
# department).

departments <- function() {
  # s2 trips on a handful of self-touching polygons in the bundled dataset ("Loop 1 is not
  # valid"); switched off only for this construction.
  s2_was_on <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(s2_was_on))

  d <- maps::map("france", plot = FALSE, fill = TRUE) |> sf::st_as_sf() |> sf::st_make_valid()
  sf::st_crs(d) <- CRS$geographic
  sf::st_transform(d, CRS$projected)
}
france <- departments()

# Green outline as a size reference. The petite couronne (Paris, Hauts-de-Seine,
# Seine-Saint-Denis, Val-de-Marne) sits entirely within 27 km of the trap, so its full ring
# shows on every panel here; Ile-de-France's own boundary starts beyond 30 km and would be
# invisible on the smaller ones.
GREEN_REF_DEPARTMENTS <- c("Paris", "Hauts-de-Seine", "Seine-Saint-Denis", "Val-de-Marne")

green_ref <- local({
  s2_was_on <- sf::sf_use_s2(); sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(s2_was_on))
  present <- france$ID %in% GREEN_REF_DEPARTMENTS
  if (sum(present) != length(GREEN_REF_DEPARTMENTS)) {
    warning("Reference outline: matched ", sum(present), " of ",
            length(GREEN_REF_DEPARTMENTS), " departments by name", call. = FALSE)
  }
  sf::st_union(france[present, ])
})

# Every fix, not thinned: each panel draws one bird alone, so thin_hourly()'s
# density-weighting problem doesn't apply. Rendering up to 40 505 raw fixes looked heavy
# but the real cost was a Windows graphics-device bug, not data volume.

#' Half-width of a country-scale panel, sized on the bird's own extent AND the trap, so
#' the frame always shows the link back to Paris.
national_panel_half <- function(x, y) {
  0.58 * max(diff(range(c(x, trap$x))), diff(range(c(y, trap$y))))
}

national_panel <- function(d, title, subtitle) {
  half <- national_panel_half(d$x, d$y)
  cx <- mean(range(c(d$x, trap$x))); cy <- mean(range(c(d$y, trap$y)))
  box <- square_box(cx, cy, half)

  ggplot() +
    geom_sf(data = france, fill = "#f2f0ec", color = "gray75", linewidth = 0.2) +
    geom_sf(data = green_ref, fill = NA, color = "#4f8a53", linewidth = 0.6) +
    geom_path(data = d, aes(x = x, y = y, color = t), linewidth = 0.5, alpha = 0.3,
              lineend = "round", linejoin = "round") +
    geom_point(data = trap, aes(x = x, y = y), shape = 4, size = 2, stroke = 0.9,
               color = "gray15") +
    seasonal_scale +
    guides(color = "none") +
    ggspatial::annotation_scale(location = "br", line_width = 0.5, height = unit(0.1, "cm"),
                                text_cex = 0.55, text_col = "gray30",
                                bar_cols = c("gray30", "white"), pad_x = unit(0.15, "cm"),
                                pad_y = unit(0.15, "cm")) +
    coord_sf(xlim = c(box[["xmin"]], box[["xmax"]]), ylim = c(box[["ymin"]], box[["ymax"]]),
             expand = FALSE, crs = st_crs(CRS$projected), datum = NA) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_corvus() +
    theme(axis.text = element_blank(), panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(size = rel(0.9)),
          plot.subtitle = element_text(size = rel(0.9)))
}

# Picked by inspecting daily-max distance against day of year for a spread of shapes, not
# a mechanical quantile pick: several of the birds beyond 40 km share a similar maximum
# distance despite very different journeys.

# Two lines, ring then numbers: one line ran ~75 characters into panels that fit ~60, and
# patchwork doesn't clip a plot against its un-clipped neighbor.
dispersal_far <- tribble(
  ~ring,          ~stats,
  "FRP-EC113554", "15 days beyond 15 km, 14 of them in April-May",
  "FRP-EC112904", "17 days beyond 15 km, all in March-May"
)

# Counted, not written down, so it can't go stale.
n_far <- scy |> group_by(ring) |>
  summarize(max_km = max(dist_trap_km), .groups = "drop") |>
  summarize(n = sum(max_km > 40), total = n())
message(sprintf("Birds ranging beyond 40 km from the trap: %d of %d.", n_far$n, n_far$total))

p_dispersal_far <- wrap_plots(lapply(seq_len(nrow(dispersal_far)), function(i) {
  d <- filter(scy, ring == dispersal_far$ring[i]) |> arrange(timestamp) |>
    mutate(t = yday)
  national_panel(d, field_ring_by_ring[[dispersal_far$ring[i]]],
                 dispersal_far$stats[i])
}), nrow = 1)

p_dispersal_far <- (p_dispersal_far / plot_spacer() / season_strip) +
  plot_layout(heights = c(1, 0.0225, 0.11)) +
  plot_annotation(
    title    = "Figure 8: Examples of long-distance movement in spring",
    subtitle = "Every location in the second calendar year, joined in time order, colored by day of year",
    caption  = paste(
      "Gray lines are French department boundaries; green is the petite couronne",
      "(Paris, Hauts-de-Seine, Seine-Saint-Denis, Val-de-Marne).",
      "The cross marks the trapping site.",
      sep = "\n"),
    theme = theme_corvus())

save_fig(p_dispersal_far, "fig8_long_dist_examples", height_cm = 17)

message("\nSpatial coverage table: ", file.path(PATHS$processed_dir, "spatial_coverage.csv"))

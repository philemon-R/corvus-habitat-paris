## utils.R -- helpers shared across the pipeline. Sourced by every numbered script.

#' Assign a hatching cohort from a deployment start date: month >= cutoff_month belongs to
#' that year's cohort, earlier months to the previous year's. See config.R (SELECTION).
#'
#' @param date Date or POSIXct vector of deployment start dates.
#' @param cutoff_month Integer month at which the cohort year rolls over.
#' @return Integer vector of cohort (hatching) years.
cohort_from_date <- function(date, cutoff_month = SELECTION$cohort_cutoff_month) {
  date <- as.Date(date)
  year <- as.integer(format(date, "%Y"))
  month <- as.integer(format(date, "%m"))
  ifelse(month >= cutoff_month, year, year - 1L)
}

#' Check if a deployment is part of the standard autumn/winter trapping protocol.
#'
#' @param date Date or POSIXct vector of deployment start dates.
#' @return Logical vector.
is_standard_protocol <- function(date) {
  month <- as.integer(format(as.Date(date), "%m"))
  month %in% SELECTION$tagging_months
}

#' Extract the stable individual identifier. Movebank's `tag-id` is not stable (tags are
#' redeployed across birds); the ring bracketed in `animal-id` ("G279 [FRP-EA739150]") is
#' authoritative, since the standalone ring column has occasional typos.
#'
#' @param animal_id Character vector of Movebank animal identifiers.
#' @return Character vector of ring codes.
ring_from_animal_id <- function(animal_id) {
  bracketed <- sub("^.*\\[(.+?)\\].*$", "\\1", animal_id)
  # Normalize the hyphen/underscore inconsistency between ring code batches.
  bracketed <- gsub("_", "-", bracketed)
  trimws(bracketed)
}

#' Extract the field-readable color-ring code ("G279" in "G279 [FRP-EA739150]").
#'
#' @param animal_id Character vector of Movebank animal identifiers.
#' @return Character vector of field-ring codes.
field_ring_from_animal_id <- function(animal_id) {
  trimws(sub("^(.*?)\\s*\\[.+\\]\\s*$", "\\1", animal_id))
}

#' Drop units from a vector, if it carries any. move2 attributes arrive as `units` objects
#' wherever Movebank declares one (HDOP, speeds, heights); comparing to a bare number
#' errors rather than coercing.
#'
#' @param x A vector, with or without units.
#' @return A plain numeric vector.
drop_units_if_any <- function(x) {
  if (inherits(x, "units")) as.numeric(units::drop_units(x)) else as.numeric(x)
}

#' Flag fixes where the incoming and outgoing value both exceed a threshold.
#'
#' @param x Numeric vector, one value per fix, in track order.
#' @param threshold Numeric scalar.
#' @param track_id Vector identifying the track each fix belongs to, in the same order.
#' @return Logical vector, TRUE where the fix should be dropped.
flag_both_sides <- function(x, threshold, track_id) {
  outgoing <- as.numeric(x)
  incoming <- c(NA, head(outgoing, -1))
  incoming[!duplicated(as.character(track_id))] <- NA
  !is.na(outgoing) & !is.na(incoming) & outgoing > threshold & incoming > threshold
}

#' The capture site, projected, as a one-row table of x and y.
trap_point <- function() {
  xy <- sf::st_coordinates(sf::st_transform(
    sf::st_sfc(sf::st_point(c(STUDY$trap_lon, STUDY$trap_lat)), crs = CRS$geographic),
    CRS$projected))
  data.frame(x = xy[1], y = xy[2])
}

#' Per-calendar-date night bounds at the trap (scripts 05/06/07): first whole hour after
#' sunset to last whole hour before sunrise, twilight counts as day. Computed once per date
#' at the trap coordinates, not per location.
#'
#' @param dates Date vector; duplicates and order do not matter.
#' @return A data frame with one row per distinct date: date, night_starts, night_ends
#'   (POSIXct, UTC).
night_bounds <- function(dates) {
  suncalc::getSunlightTimes(date = sort(unique(dates)),
                            lat = STUDY$trap_lat, lon = STUDY$trap_lon,
                            keep = c("sunrise", "sunset"), tz = "UTC") |>
    dplyr::transmute(date,
              night_starts = as.POSIXct(ceiling(as.numeric(sunset) / 3600) * 3600,
                                        origin = "1970-01-01", tz = "UTC"),
              night_ends   = as.POSIXct(floor(as.numeric(sunrise) / 3600) * 3600,
                                        origin = "1970-01-01", tz = "UTC"))
}

#' Create a directory if it does not exist, silently.
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

#' Stop with a clear message if a pipeline input is missing.
require_input <- function(path, produced_by) {
  if (!file.exists(path)) {
    stop(sprintf("Missing input '%s'. Run %s first.", path, produced_by), call. = FALSE)
  }
  invisible(TRUE)
}

#' Thin a location table to at most one fix per individual per clock hour.
#' The fix kept is the one nearest the middle of its hour.
#'
#' @param df A data frame with one row per location.
#' @param id Unquoted column identifying the individual.
#' @param time Unquoted POSIXct column.
#' @return `df` with the surplus rows removed, in the original order.
thin_hourly <- function(df, id, time) {
  who  <- as.character(dplyr::pull(df, {{ id }}))
  when <- as.numeric(dplyr::pull(df, {{ time }}))
  hour <- floor(when / 3600)

  # Order by (individual, hour, distance to the middle of that hour), then keep the first
  # row of each (individual, hour) run.
  o <- order(who, hour, abs(when - (hour * 3600 + 1800)))
  w <- who[o]; h <- hour[o]
  first <- c(TRUE, w[-1] != w[-length(w)] | h[-1] != h[-length(h)])
  df[sort(o[first]), , drop = FALSE]
}

#' Pool per-individual estimates via random-effects meta-analysis (script 07's second
#' fitting stage), not a fixed-effect inverse-variance mean: fixed-effect assumes one true
#' coefficient shared by every bird and lets birds with more steps dominate the
#' aggregation.
#'
#' DerSimonian-Laird estimates the between-bird variance (tau^2) from the spread in excess
#' of sampling error, then weights each bird by 1 / (se^2 + tau^2): where birds agree,
#' tau^2 -> 0 and this reduces to fixed-effect; where they disagree, weights flatten toward
#' one-bird-one-vote and the interval widens. `metafor::rma(method = "DL")` implements this.
#'
#' @param estimate,se Per-individual coefficient and its standard error, same length.
#' @return A one-row data frame: pooled estimate and se, tau2, Cochran's Q, I^2 (the share
#'   of total variation that is between-individual rather than sampling), and n.
meta_pool <- function(estimate, se) {
  ok <- is.finite(estimate) & is.finite(se) & se > 0
  estimate <- estimate[ok]; se <- se[ok]
  n <- length(estimate)
  if (n == 0) {
    return(data.frame(estimate = NA_real_, se = NA_real_, tau2 = NA_real_,
                      q = NA_real_, i2 = NA_real_, n = 0L))
  }
  if (n == 1) {
    return(data.frame(estimate = estimate, se = se, tau2 = 0, q = 0, i2 = 0, n = 1L))
  }

  fit <- metafor::rma(yi = estimate, sei = se, method = "DL")
  data.frame(estimate = as.numeric(fit$b),
             se   = fit$se,
             tau2 = fit$tau2,
             q    = fit$QE,
             i2   = fit$I2 / 100,
             n    = n)
}

#' Flag the tail of a track after the animal last demonstrably moved.
#'
#' Needs positive evidence of stopping, not merely an absence of evidence of moving: a
#' fading battery reporting two fixes a day should not by itself end a track. So a day
#' counts only if it carries at least `min_fixes`; the tail is cut after the last such day
#' showing real movement, only once a later day positively shows the tag sitting still.
#'
#' @param x,y Projected coordinates, meters.
#' @param timestamp POSIXct.
#' @param track_id Track identifier.
#' @param still_m A day spanning less than this is "still".
#' @param min_fixes Fixes needed before a day can be judged either way.
#' @return Logical vector, TRUE for fixes to drop.
flag_stationary_tail <- function(x, y, timestamp, track_id,
                                 still_m = QC$stationary_tag_m,
                                 min_fixes = QC$stationary_tag_min_fixes) {
  d <- data.frame(id = as.character(track_id), date = as.Date(timestamp),
                  x = as.numeric(x), y = as.numeric(y))
  key <- paste(d$id, d$date)

  spread <- stats::ave(seq_len(nrow(d)), key, FUN = function(i) {
    if (length(i) < min_fixes) return(NA_real_)
    sqrt((max(d$x[i]) - min(d$x[i]))^2 + (max(d$y[i]) - min(d$y[i]))^2)
  })

  moved <- !is.na(spread) & spread >= still_m
  still <- !is.na(spread) & spread <  still_m

  # Per track: the last date with demonstrated movement, and whether anything after it is
  # positively still. Dates rather than timestamps, so a whole day goes or stays together.
  last_moved <- stats::ave(seq_len(nrow(d)), d$id, FUN = function(i) {
    m <- moved[i]
    if (!any(m)) return(-Inf)      # never demonstrated movement, nothing to cut after
    as.numeric(max(d$date[i][m]))
  })
  still_after <- stats::ave(seq_len(nrow(d)), d$id, FUN = function(i) {
    as.numeric(any(still[i] & as.numeric(d$date[i]) > last_moved[i][1]))
  })

  as.numeric(d$date) > last_moved & still_after > 0
}

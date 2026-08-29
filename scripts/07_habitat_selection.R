## 07_habitat_selection.R -- step-selection analysis
##
## Inputs:  data/interim/tracks_clean.rds       (script 02)
##          data/osm/distance_stack.tif         (script 04)
## Outputs: data/processed/ssf_coefficients.csv -- pooled coefficients, both panels
##          data/processed/ssf_individual.csv   -- the per-bird estimates behind them
##          data/interim/ssf_models.rds         -- per-bird fits, stripped of row-level parts
##          data/interim/ssf_random_steps_cache/<hash>/<ring>.rds -- cached random-step
##                                                 geometry, reused by a re-run that only
##                                                 changes SSF$covariates or a covariate's
##                                                 radius (see section 6)
##          figures/fig9_habitat_selection.png            -- cross-sectional selection, first year
##          figures/fig10_habitat_selection_with_age.png  -- how that selection changes with age
##
## The movement terms (sl_, log(sl_), cos(ta_)) make this an integrated SSF: random steps
## are drawn from a kernel fitted to the observed steps.
##
## Usage:
##   Rscript scripts/07_habitat_selection.R

suppressPackageStartupMessages({
  library(move2); library(sf); library(dplyr); library(tidyr); library(purrr)
  library(readr); library(amt); library(lubridate); library(survival); library(terra)
  library(ggplot2); library(patchwork); library(suncalc)
})

source("config.R")
source(file.path("scripts", "utils.R"))
source(file.path("scripts", "plotting.R"))

require_input(PATHS$clean_tracks, "scripts/02_clean_gps.R")
require_input(PATHS$covariate_stack, "scripts/04_covariates.R")
ensure_dir(PATHS$processed_dir)
ensure_dir(PATHS$interim_dir)
ensure_dir(PATHS$figures_dir)

set.seed(RANDOM_SEED)

covariates <- terra::rast(PATHS$covariate_stack)
missing_layers <- setdiff(SSF$covariates, names(covariates))
if (length(missing_layers)) {
  stop("Covariate stack lacks: ", paste(missing_layers, collapse = ", "),
       ". Delete ", PATHS$covariate_stack, " and re-run scripts/04_covariates.R.",
       call. = FALSE)
}

# ---------------------------------------------------------------------------------
# 1. Prepare age variable and panel info
# ---------------------------------------------------------------------------------

tracks <- readRDS(PATHS$clean_tracks)
xy <- st_coordinates(tracks)
dat <- tibble(
  ring        = as.character(mt_track_id(tracks)),
  t           = mt_time(tracks),
  x           = xy[, "X"],
  y           = xy[, "Y"],
  indep_years = tracks$indep_years
)
individuals <- read_csv(PATHS$individuals, show_col_types = FALSE) |>
  select(ring, cohort, panel_b, span_months)

message(sprintf("Loaded %s locations for %d individuals.",
                format(nrow(dat), big.mark = " "), n_distinct(dat$ring)))

# ---------------------------------------------------------------------------------
# 2. Daytime restriction
# ---------------------------------------------------------------------------------
# night_bounds() (utils.R), shared with scripts 05/06. Applied to raw fixes before
# track_resample()/steps_for().

if (SSF$restrict_daytime) {
  n_before <- nrow(dat)
  sun <- night_bounds(as.Date(dat$t))
  dat <- dat |>
    mutate(date = as.Date(t)) |>
    left_join(sun, by = "date") |>
    filter(t >= night_ends, t < night_starts) |>
    select(-date, -night_starts, -night_ends)

  message(sprintf("Daytime restriction: %s of %s locations kept (%.1f%%).",
                  format(nrow(dat), big.mark = " "), format(n_before, big.mark = " "),
                  100 * nrow(dat) / n_before))
}

# ---------------------------------------------------------------------------------
# 3. Steps
# ---------------------------------------------------------------------------------
# Resampling is stricter than the maps' hourly thinning (thin_hourly(), utils.R): a step
# needs a known duration and an unbroken predecessor, not just a comparable count. Two
# further filters forced by the model's own terms: a zero-length step has log(sl_) = -Inf,
# and a burst's first step has no turning angle (cos(ta_) is a parameter).

steps_for <- function(bird) {
  d <- dat[dat$ring == bird, ]
  trk <- make_track(d, x, y, t, indep_years = indep_years,
                    crs = CRS$projected)
  s <- trk |>
    track_resample(rate = minutes(SSF$resample_rate_min),
                   tolerance = minutes(SSF$resample_tolerance_min)) |>
    filter_min_n_burst(min_n = SSF$min_burst_length)
  if (nrow(s) == 0) return(NULL)
  s <- steps_by_burst(s, keep_cols = "start")
  s <- s[s$sl_ > 0 & !is.na(s$ta_), ]
  if (nrow(s) == 0) NULL else s
}

birds <- sort(unique(dat$ring))
observed <- map(set_names(birds), steps_for)
observed <- compact(observed)
message(sprintf("Steps: %s from %d birds.",
                format(sum(map_int(observed, nrow)), big.mark = " "), length(observed)))

# ---------------------------------------------------------------------------------
# 4. Scaling constants, computed once over every bird
# ---------------------------------------------------------------------------------
# Center/scale come from the pooled observed endpoints, applied unchanged to every bird
# including its random steps -- per-bird scaling would give each bird's coefficients their
# own units, and section 7's pooling would then add up different quantities.

#' Distance covariates are negated after logging so every term reads the same way:
#' positive means more selected. A "dist_" layer's raw value stays available (e.g. for a
#' map), but the model sees -log1p(distance) -- proximity, agreeing with a density
#' covariate's sign.
prepare_covariates <- function(d) {
  for (v in SSF$covariates) {
    x <- d[[v]]
    if (SSF$log_distances && startsWith(v, "dist_")) x <- -log1p(x)
    d[[paste0("z_", v)]] <- x
  }
  d
}

endpoints <- bind_rows(map(observed, ~ tibble(x = .x$x2_, y = .x$y2_)))
ref <- terra::extract(covariates[[SSF$covariates]],
                      cbind(endpoints$x, endpoints$y)) |>
  as_tibble() |>
  prepare_covariates()

scaling <- tibble(
  covariate = SSF$covariates,
  center = map_dbl(SSF$covariates, ~ mean(ref[[paste0("z_", .x)]], na.rm = TRUE)),
  scale  = map_dbl(SSF$covariates, ~ sd(ref[[paste0("z_", .x)]], na.rm = TRUE))
)
message("\nScaling constants (-log1p meters, i.e. proximity, for distances; fraction for density):")
for (i in seq_len(nrow(scaling))) {
  message(sprintf("  %-20s center %6.2f  sd %5.2f",
                  scaling$covariate[i], scaling$center[i], scaling$scale[i]))
}
rm(ref, endpoints)

# ---------------------------------------------------------------------------------
# 5. Compute available steps, covariates and fit for each bird
# ---------------------------------------------------------------------------------
#
# The movement kernel is fitted per bird on the bird's whole track.
#
# A stratum with any missing covariate is dropped whole, not row by row: dropping only the
# row for a random step outside the study area would redefine availability as
# "reachable and inside the box".

TERMS_MOVEMENT <- c("sl_", "log_sl_", "cos_ta_")

# ---------------------------------------------------------------------------------
# 6. Random-step cache
# ---------------------------------------------------------------------------------
# random_steps() depends only on steps_for()'s resampled geometry and SSF$n_random_steps.
# A covariate-only re-run doesn't have to refit and redraw the kernel.
#
# The cache key HASHES what steps_for() depends on (tracks_clean.rds's content, its own
# deparsed source, the resampling params), plus SSF$restrict_daytime since that changes `dat`
# upstream of steps_for() in a way neither hash sees.
#
# One file per bird: the per-bird loop below never holds more than one bird's steps at
# once, and a single combined file would force the whole several-GB table into memory.
ssf_cache_key <- local({
  tracks_hash <- unname(tools::md5sum(PATHS$clean_tracks))
  key_material <- c(
    tracks_hash, deparse(steps_for),
    as.character(SSF$resample_rate_min), as.character(SSF$resample_tolerance_min),
    as.character(SSF$min_burst_length), as.character(SSF$n_random_steps),
    as.character(SSF$restrict_daytime)
  )
  tmp <- tempfile()
  writeLines(key_material, tmp)
  key <- unname(tools::md5sum(tmp))
  file.remove(tmp)
  key
})
ssf_cache_dir <- file.path(PATHS$interim_dir, "ssf_random_steps_cache", ssf_cache_key)
ensure_dir(ssf_cache_dir)

#' Observed steps plus SSF$n_random_steps random ones, drawn from the kernel fitted to
#' `s`. Cached to ssf_cache_dir/<bird>.rds -- a cache hit reuses the same draws, which
#' keeps a bird's random steps stable across covariate-only re-runs.
random_steps_cached <- function(bird, s) {
  path <- file.path(ssf_cache_dir, paste0(bird, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  rs <- random_steps(s, n_control = SSF$n_random_steps)
  saveRDS(rs, path)
  rs
}

#' Read covariates at every step endpoint, scale them, and add the terms fit_one() needs.
#' Split from random_steps_cached() so a covariate- or radius-only re-run can skip straight
#' here with a cached `rs`.
attach_covariates <- function(rs) {
  rs <- extract_covariates(rs, covariates[[SSF$covariates]], where = "end")
  rs <- prepare_covariates(rs)

  z <- paste0("z_", SSF$covariates)
  for (i in seq_along(z)) {
    rs[[z[i]]] <- if (SSF$scale_covariates) {
      (rs[[z[i]]] - scaling$center[i]) / scaling$scale[i]
    } else rs[[z[i]]]
  }

  bad <- unique(rs$step_id_[rowSums(is.na(as.data.frame(rs[z]))) > 0])
  rs <- rs[!rs$step_id_ %in% bad, ]

  rs$log_sl_ <- log(rs$sl_)
  rs$cos_ta_ <- cos(rs$ta_)
  attr(rs, "n_strata_dropped") <- length(bad)
  rs
}

#' Fit one conditional logistic model and return its coefficients, or NULL if it fails.
fit_one <- function(rs, terms, bird, panel) {
  # Built by pasting, not reformulate(), which back-quotes non-syntactic terms.
  f <- stats::as.formula(
    paste("case_ ~", paste(c(terms, "strata(step_id_)"), collapse = " + ")))
  m <- tryCatch(clogit(f, data = rs), error = function(e) e)
  if (inherits(m, "error")) {
    message(sprintf("  %s: %s fit failed (%s)", bird, panel, conditionMessage(m)))
    return(NULL)
  }
  co <- summary(m)$coefficients
  if (any(!is.finite(co[, "coef"])) || any(!is.finite(co[, "se(coef)"]))) {
    message(sprintf("  %s: %s fit did not converge to finite estimates", bird, panel))
    return(NULL)
  }
  # Trim unused row-level components.
  m$residuals <- NULL; m$linear.predictors <- NULL; m$y <- NULL
  # Remove pointers to the fitting environment to save memory and space.
  if (!is.null(m$formula)) environment(m$formula) <- baseenv()
  if (!is.null(m$terms)) attr(m$terms, ".Environment") <- baseenv()
  if (!is.null(m$call$formula)) environment(m$call$formula) <- baseenv()
  list(
    coefficients = tibble(ring = bird, panel = panel, term = rownames(co),
                          estimate = co[, "coef"], se = co[, "se(coef)"],
                          n_steps = sum(rs$case_), n_rows = nrow(rs)),
    model = m
  )
}

panel_b_birds <- individuals$ring[individuals$panel_b %in% TRUE]

results <- list()
models  <- list()
skipped <- list()
dropped_strata <- 0L

for (bird in names(observed)) {
  s <- observed[[bird]]
  rs <- attach_covariates(random_steps_cached(bird, s))
  dropped_strata <- dropped_strata + attr(rs, "n_strata_dropped")

  # -- Panel A: the year of independence only --------------------------------------
  # indep_years is 0 exactly at independence (config.R, AGE), so this window is a full
  # year starting there.
  a <- rs[rs$indep_years >= 0 & rs$indep_years < 1, ]
  if (n_distinct(a$step_id_) >= SSF$min_steps_per_bird) {
    fit <- fit_one(a, c(paste0("z_", SSF$covariates), TERMS_MOVEMENT), bird, "A")
    if (!is.null(fit)) {
      results[[length(results) + 1L]] <- fit$coefficients
      models[[paste0("A_", bird)]] <- fit$model
    }
  } else {
    skipped[[length(skipped) + 1L]] <-
      tibble(ring = bird, panel = "A", n_steps = n_distinct(a$step_id_),
             reason = "too few year-of-independence steps")
  }

  # -- Panel B: from independence onward, with age interacting with habitat --------
  # Age is constant within a stratum (an observed step and its random steps share a start
  # instant), so clogit can't estimate an age main effect -- only the interaction.
  #
  # Restricted to indep_years >= 0: pre-independence is a different behavioral regime
  # (still dependent on parents). Centered at indep_years == 1 so the habitat term reads as
  # selection at the exact instant Panel A's window ends, and the interaction as change per
  # further year.
  if (bird %in% panel_b_birds) {
    b <- rs[rs$indep_years >= 0, ]
    age_span <- diff(range(b$indep_years))
    if (n_distinct(b$step_id_) >= SSF$min_steps_per_bird &&
        age_span >= SSF$min_age_span_years) {
      b$age_c <- b$indep_years - SSF$age_center_years
      z <- paste0("z_", SSF$covariates)
      fit <- fit_one(b, c(z, paste0(z, ":age_c"), TERMS_MOVEMENT), bird, "B")
      if (!is.null(fit)) {
        results[[length(results) + 1L]] <- fit$coefficients
        models[[paste0("B_", bird)]] <- fit$model
      }
    } else {
      skipped[[length(skipped) + 1L]] <-
        tibble(ring = bird, panel = "B", n_steps = n_distinct(b$step_id_),
               reason = if (age_span < SSF$min_age_span_years)
                 sprintf("age span %.2f yr", age_span) else "too few steps")
    }
  }

  rm(rs, s)
}

per_bird <- bind_rows(results)
skipped  <- bind_rows(skipped)

message(sprintf("\nStrata dropped for a covariate outside the study area: %s.",
                format(dropped_strata, big.mark = " ")))
message(sprintf("Panel A: %d birds fitted. Panel B: %d birds fitted.",
                n_distinct(per_bird$ring[per_bird$panel == "A"]),
                n_distinct(per_bird$ring[per_bird$panel == "B"])))
if (nrow(skipped)) {
  message("Not fitted:")
  for (i in seq_len(nrow(skipped))) {
    message(sprintf("  %s panel %s: %s (%d steps)", skipped$ring[i], skipped$panel[i],
                    skipped$reason[i], skipped$n_steps[i]))
  }
}

per_bird <- per_bird |> left_join(individuals |> select(ring, cohort), by = "ring")
write_csv(per_bird, PATHS$ssf_individual)
saveRDS(models, PATHS$ssf_models)

# ---------------------------------------------------------------------------------
# 7. Pooling
# ---------------------------------------------------------------------------------
#
# Pooling is random-effects (DerSimonian-Laird, meta_pool() in utils.R): a bird with many
# steps can't dominate, and real between-bird differences widen the interval instead of
# being absorbed as noise. I^2 reports how much of the spread is between birds.

pooled <- per_bird |>
  group_by(panel, term) |>
  group_modify(~ meta_pool(.x$estimate, .x$se)) |>
  ungroup() |>
  mutate(
    conf_low   = estimate - 1.96 * se,
    conf_high  = estimate + 1.96 * se,
    # Relative selection strength: multiplicative change in the odds a step ends here, per
    # SD of the covariate, relative to the stratum's other endpoints -- never an absolute
    # probability of use.
    rss        = exp(estimate),
    rss_low    = exp(conf_low),
    rss_high   = exp(conf_high),
    p_value    = 2 * pnorm(-abs(estimate / se))
  )

write_csv(pooled, PATHS$ssf_coefficients)

show <- pooled |>
  filter(!term %in% TERMS_MOVEMENT) |>
  arrange(panel, term)
message("\nPooled coefficients (log relative selection per SD):")
for (i in seq_len(nrow(show))) {
  message(sprintf("  %s  %-28s %+6.3f [%+6.3f, %+6.3f]  I2 = %2.0f%%  n = %d",
                  show$panel[i], show$term[i], show$estimate[i],
                  show$conf_low[i], show$conf_high[i], 100 * show$i2[i], show$n[i]))
}

# ---------------------------------------------------------------------------------
# 8. Figures
# ---------------------------------------------------------------------------------

#' "Proximity", not "distance", because prepare_covariates() negates dist_-prefixed
#' covariates so every term shares one sign convention.
COVARIATE_LABELS <- c(
  dist_agricultural          = "Proximity to farmland",
  dist_water                 = "Proximity to water",
  woody_density_200m         = "Woody canopy density, 200 m",
  nonagri_herb_density_200m  = "Non-agricultural herbaceous density, 200 m",
  built_density_200m         = "Built-up density, 200 m"
)

PLOT_COVARIATE_ORDER <- c("woody_density_200m", "nonagri_herb_density_200m",
                           "built_density_200m", "dist_agricultural", "dist_water")

# Term names come back from the model as "z_dist_agricultural" and
# "z_dist_agricultural:age_c"; this strips both wrappers so a label can be looked up.
covariate_of <- function(term) sub("^z_", "", sub(":age_c$", "", term))
label_factor <- function(term) {
  factor(COVARIATE_LABELS[covariate_of(term)],
         levels = rev(unname(COVARIATE_LABELS[PLOT_COVARIATE_ORDER])))
}

# Points beyond the range are dropped before plotting, not cropped, so none spill 
# into the neighboring panel; captions report how many.
FIG_A_XLIM <- c(-3, 3)
FIG_B_XLIM <- c(-2, 2)

# Space above the top row for the direction cue. Cue sits at CUE_FRACTION up that gap.
TOP_PAD <- 0.825
BOTTOM_PAD <- 0.6
CUE_FRACTION <- 0.75

n_outside <- function(bird_rows, lim) {
  sum(bird_rows$estimate < lim[1] | bird_rows$estimate > lim[2], na.rm = TRUE)
}

#' Add "avoidance <- -> preference" cue above the panel.
add_direction_cue <- function(p, lim, y0) {
  gap <- diff(lim) * 0.05
  arrow_end <- diff(lim) * 0.11
  p +
    annotate("segment", x = -gap, xend = -arrow_end, y = y0, yend = y0,
             arrow = arrow(length = unit(0.16, "cm"), type = "open"),
             linewidth = 0.4, color = "gray30") +
    annotate("segment", x = gap, xend = arrow_end, y = y0, yend = y0,
             arrow = arrow(length = unit(0.16, "cm"), type = "open"),
             linewidth = 0.4, color = "gray30") +
    annotate("text", x = -arrow_end * 1.15, y = y0, label = "avoidance", hjust = 1,
             vjust = 0.35, fontface = "bold", size = 3.1, color = "gray20") +
    annotate("text", x = arrow_end * 1.15, y = y0, label = "preference", hjust = 0,
             vjust = 0.35, fontface = "bold", size = 3.1, color = "gray20")
}

coefficient_panel <- function(pooled_rows, bird_rows, xlim, title = NULL,
                               direction_cue = FALSE) {
  lim <- xlim
  bird_rows <- bird_rows[bird_rows$estimate >= lim[1] & bird_rows$estimate <= lim[2], ]
  n_levels <- nlevels(pooled_rows$label)

  p <- ggplot(pooled_rows, aes(x = estimate, y = label)) +
    geom_vline(xintercept = 0, color = "gray45", linewidth = 0.4) +
    # The per-bird cloud with jitter.
    geom_point(data = bird_rows, aes(x = estimate, y = label), color = "gray60",
               size = 0.7, alpha = 0.45,
               position = position_jitter(height = 0.16, width = 0, seed = RANDOM_SEED)) +
    # Line + white-rimmed point for the estimate.
    geom_linerange(aes(xmin = conf_low, xmax = conf_high),
                   linewidth = 0.9, color = "#0072B2") +
    geom_point(shape = 21, fill = "#0072B2", color = "white", stroke = 0.7, size = 2.6) +
    # Breaks from the clipped range
    scale_x_continuous(breaks = pretty(lim, n = 5)) +
    # Explicit expansion
    scale_y_discrete(expand = expansion(add = c(BOTTOM_PAD, TOP_PAD))) +
    labs(title = title, x = NULL, y = NULL) +
    theme_corvus() +
    theme(panel.grid.major.y = element_blank())

  if (!is.null(title)) {
    # theme_corvus() aligns plot.title to the whole figure, right for fig 9's headline but
    # wrong for fig 10's two side-by-side panel titles, which need to center on THIS
    # panel's x = 0.
    p <- p + theme(plot.title.position = "panel", plot.title = element_text(hjust = 0.5))
  }

  # Same top pad whether or not this panel draws the cue: fig 10's two panels are
  # patchworked side by side and need equal row heights.
  p <- p + coord_cartesian(xlim = lim)
  if (direction_cue) p <- add_direction_cue(p, lim, n_levels + TOP_PAD * CUE_FRACTION)
  p
}

habitat_terms <- paste0("z_", SSF$covariates)
steps_in <- function(p) sum(per_bird$n_steps[per_bird$panel == p &
                                             per_bird$term == habitat_terms[1]])

pa <- pooled |> filter(panel == "A", term %in% habitat_terms) |>
  mutate(label = label_factor(term))
ba <- per_bird |> filter(panel == "A", term %in% habitat_terms) |>
  mutate(label = label_factor(term))

p1 <- coefficient_panel(pa, ba, FIG_A_XLIM, direction_cue = TRUE) +
  labs(
    title = "Figure 9: Habitat selection in the year of independence",
    subtitle = sprintf("%d birds, %s steps, each against %d random steps from the same point",
                       n_distinct(ba$ring), format(steps_in("A"), big.mark = " "),
                       SSF$n_random_steps),
    x = "Change in log relative selection per standard deviation",
    caption = sprintf(paste(
      "Blue: the population estimate and its 95%% interval, pooled across all birds.",
      "Gray: one point per bird, axis fixed to %g to %g and %d of %d birds outside it",
      "not shown.",
      sep = "\n"), FIG_A_XLIM[1], FIG_A_XLIM[2], n_outside(ba, FIG_A_XLIM), nrow(ba)))

save_fig(p1, "fig9_habitat_selection", height_cm = 13)

# Panel B: the same six covariates twice over. The left column is where selection sits one
# year after independence, the right is how it moves per further year.
pb_main <- pooled |> filter(panel == "B", term %in% habitat_terms) |>
  mutate(label = label_factor(term))
bb_main <- per_bird |> filter(panel == "B", term %in% habitat_terms) |>
  mutate(label = label_factor(term))
age_terms <- paste0(habitat_terms, ":age_c")
pb_age <- pooled |> filter(panel == "B", term %in% age_terms) |>
  mutate(label = label_factor(term))
bb_age <- per_bird |> filter(panel == "B", term %in% age_terms) |>
  mutate(label = label_factor(term))

left  <- coefficient_panel(pb_main, bb_main, FIG_B_XLIM, "One year after independence",
                            direction_cue = TRUE) +
  labs(x = "Log relative selection per SD")
right <- coefficient_panel(pb_age, bb_age, FIG_B_XLIM, "Change per additional year") +
  labs(x = "Change in that value per year") +
  theme(axis.text.y = element_blank())

p2 <- (left | right) +
  plot_annotation(
    title = "Figure 10: Does habitat selection change as a crow ages?",
    subtitle = sprintf("%d birds tracked at least %g months past independence, %s steps",
                       n_distinct(bb_main$ring), SELECTION$min_post_indep_months_b,
                       format(steps_in("B"), big.mark = " ")),
    caption = sprintf(paste(
      "Age is continuous and centered one year after independence, so the left panel reads",
      "as selection at that point and the right as its change per further year.",
      "Blue: the population estimate and its 95%% interval, pooled across all birds.",
      "Gray: one point per bird, both axes fixed to %g to %g, %d and %d of %d birds outside",
      "them not shown.",
      sep = "\n"), FIG_B_XLIM[1], FIG_B_XLIM[2],
      n_outside(bb_main, FIG_B_XLIM), n_outside(bb_age, FIG_B_XLIM), nrow(bb_main)),
    theme = theme_corvus())

save_fig(p2, "fig10_habitat_selection_with_age", height_cm = 13)

message("\nCoefficients:     ", PATHS$ssf_coefficients)
message("Per individual:   ", PATHS$ssf_individual)

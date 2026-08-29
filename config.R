## config.R -- central configuration for the corvus-habitat-paris pipeline. Every scope
## choice lives here so no script duplicates a constant.

# ---------------------------------------------------------------------------------
# Study
# ---------------------------------------------------------------------------------

STUDY <- list(
  # "Corvus corone [ID_PROG 883]", PI F. Jiguet
  movebank_id = 1266784970,
  species     = "Corvus corone",
  # Capture site (Jardin des Plantes).
  trap_lon    = 2.361,
  trap_lat    = 48.844
)

# Movebank serves WGS84 lon/lat (4326). Metric work (distances, step lengths, buffers, OSM
# overlays) uses Lambert-93 (2154).
CRS <- list(
  geographic = 4326,
  projected  = 2154
)

# ---------------------------------------------------------------------------------
# Individual selection
# ---------------------------------------------------------------------------------
# Restricted to first-year crows caught in the autumn/winter cage trap -- the same
# population the reference paper analyses (Jiguet & Gantin 2025). Cohort is derived from
# deployment date.

SELECTION <- list(
  cohort_cutoff_month = 7,           # July: after fledging, before autumn trapping

  # Standard tagging season; deployments starting outside it are excluded.
  tagging_months = c(9, 10, 11, 12, 1, 2, 3),

  min_positions = 300,   # minimum GPS fixes to enter the analysis

  # Panel A -- cross-sectional: all cohorts, year-of-independence window.
  cohorts_panel_a = c(2020, 2021, 2022, 2024),

  # Panel B -- longitudinal: needs enough tracking after independence to fit an age trend,
  # measured directly against each bird's own independence.
  min_post_indep_months_b = 12,

  # Manual correction to deploy_on applied in script 01 before the
  # cohort/protocol rules run.
  deploy_on_corrections = list(
    "FRP-EC111894" = as.Date("2021-10-15")
  ),

  # Birds the date/season rule would otherwise select but that are excluded for a reason
  # Jiguet & Gantin (2025).
  manual_exclusions = c(
    "FRP-EA739139",  # G272 -- not in Table S1
    "FRP-EA739150",  # G279 -- not in Table S1
    "FRP-EA739222",  # G353 -- not in Table S1
    "FRP-EC113511",  # G770 -- not in Table S1
    "FRP-EC115807"   # R160 -- sole 2023-cohort bird
  )
)

# ---------------------------------------------------------------------------------
# Quality control of GPS-tracks (script 02)
# ---------------------------------------------------------------------------------
# Thresholds mirror the reference paper.

QC <- list(
  # Thresholds mirror the reference paper.
  max_speed_kmh      = 65,     # above this, non-physiological
  max_displacement_m = 50000,  # isolated outliers

  # GPS quality
  max_hdop       = 5,
  min_satellites = 4,
  drop_poor_fixes = TRUE,

  # Stationary tag detection
  stationary_tag_m = 100,
  stationary_tag_min_fixes = 5,  # fixes a day needs before it counts as judged
  drop_stationary_tail = TRUE,

  # Target interval for hourly thinning
  nominal_fix_interval_min = 60
)

# ---------------------------------------------------------------------------------
# Age (scripts 02, 07)
# ---------------------------------------------------------------------------------
# Age is measured relative to INDEPENDENCE (31 March following hatching, when adults evict
# the brood).

AGE <- list(
  # Nominal date within the cohort year, shared by every bird in that cohort.
  independence_md = "03-31"
)

# ---------------------------------------------------------------------------------
# Settled-night detection (script 05)
# ---------------------------------------------------------------------------------

ROOST <- list(
  # Two nocturnal fixes within this distance count as the same site; used
  # only to test whether a night's fixes were tight around their own median.
  site_radius_m = 100,

  # A night counts as settled only if >50% of its fixes fall within site_radius_m of the
  # night's median.
  settled_min_share = 0.50
)

# ---------------------------------------------------------------------------------
# Habitat covariate stack (script 04): shared between the OSM and OCS GE sources below
# ---------------------------------------------------------------------------------

COVARIATES <- list(
  # Focal density radius, shared by every focal-density covariate: built-up (OSM), and
  # woody/non-agricultural-herbaceous (OCSGE$classes below).
  focal_density_radius_m = 200,

  # Cell size of the covariate stack, set by the building-density raster and shared by
  # every other layer once they're aligned to its grid.
  density_res_m = 25
)

# ---------------------------------------------------------------------------------
# OpenStreetMap habitat covariates (script 04)
# ---------------------------------------------------------------------------------
# Feature classes extracted for the study area.

OSM <- list(
  # "latest" (Geofabrik's nightly rebuild) lets a first run work with no date to look up.
  # 03_geodata_download.R records the served snapshot's URL/md5/date into
  # data/processed/osm_provenance.csv (committed) -- paste that date here to pin a re-run.
  version = "latest",

  # 62x62 km box centered on the trap.
  bbox = c(xmin = 1.94, ymin = 48.57, xmax = 2.78, ymax = 49.12),

  # Classes feeding basemap and covariates: water and building.
  features = list(
    water = list(list(key = "natural", value = "water")),
    building = list(list(key = "building", value = NULL))   # NULL = any tagged building
  )
)

# ---------------------------------------------------------------------------------
# IGN OCS GE habitat covariates (script 03)
# ---------------------------------------------------------------------------------
# Second geodata source alongside OSM: OCS GE's `couverture` dimension separates woody
# canopy (CS2.1) from herbaceous (CS2.2) and bare ground (CS1), a split OSM's tags can't
# make. Etalab-licensed, distributed per-department.

OCSGE <- list(
  # Most recent millesime available per department when this was set up; close to
  # the tracking period (2018-01-01 and 2021-01-01 also exist).
  millesime = "2024-01-01",

  # Departments the study bbox touches.
  departments = c("075", "077", "078", "091", "092", "093", "094", "095", "060"),

  # code_cs (couverture) / code_us (usage) WHERE clauses against OCCUPATION_SOL. `type`
  # picks distance-to-nearest or focal density at COVARIATES$focal_density_radius_m.
  #
  # Woody couverture codes used to define `woody` and `agricultural`.
  classes = local({
    woody_codes <- sprintf("'%s'", c("CS2.1.1.1", "CS2.1.1.2", "CS2.1.1.3", "CS2.1.2", "CS2.1.3"))
    woody_in <- paste(woody_codes, collapse = ",")
    list(
      woody = list(
        type = "density",
        where = sprintf("code_cs IN (%s)", woody_in)),
      nonagri_herb = list(
        type = "density",
        where = "code_cs = 'CS2.2.1' AND (code_us IS NULL OR code_us <> 'US1.1')"),
      agricultural = list(
        type = "distance",
        where = sprintf("code_us = 'US1.1' AND code_cs NOT IN (%s)", woody_in))
    )
  })
)

# ---------------------------------------------------------------------------------
# Step-selection analysis (script 07)
# ---------------------------------------------------------------------------------

SSF <- list(
  # Restricts to daytime fixes (same night-hour rule as scripts 05/06) before any step is
  # built.
  restrict_daytime = TRUE,

  # Tracks are regularized to the nominal fix interval before steps are built; bursts break
  # on gaps beyond tolerance.
  resample_rate_min    = 60,
  resample_tolerance_min = 15,
  min_burst_length     = 3,    # steps needed to fit turning angles

  # Random steps per observed step.
  n_random_steps = 10,

  # Movement kernel fitted from the observed steps (amt defaults).
  step_dist = "gamma",
  angle_dist = "vonmises",

  # From the script 04 stack: distance-to-nearest for sparse features (agricultural,
  # water), focal density within 200 m for denser ones (woody canopy, non-agricultural
  # herbaceous, built-up).
  covariates = c("dist_agricultural", "dist_water", "woody_density_200m",
                 "nonagri_herb_density_200m", "built_density_200m"),

  # Distances are log1p-transformed (most mass sits in the first few hundred meters).
  log_distances = TRUE,

  # Centered/scaled once over pooled step endpoints of all birds. A coefficient then
  # reads as change in log relative selection per SD of that covariate.
  scale_covariates = TRUE,

  # Minimum number of steps to fit the clogit model.
  min_steps_per_bird = 100,

  # Panel B only: minimum age spread, or the interaction slope is fit from a single point.
  min_age_span_years = 0.5,

  # Centers Panel B's `indep_years` so the habitat main effect reads as selection exactly
  # one year after independence (where Panel A's window ends), and the interaction as
  # change per further year.
  age_center_years = 1
)

# ---------------------------------------------------------------------------------
# Maps
# ---------------------------------------------------------------------------------

MAP <- list(
  # Cell size in meters for location density display.
  cell_m = 250,

  # Define core of the study area were observations are concentrated.
  core_halfwidth_m = 16000,

  # Birds for the trajectory panels are picked by rule, not eye: an even spread across
  # range distance among those with good second-calendar-year coverage.
  track_min_days   = 300,   # days of the second calendar year that must carry a fix
  track_n_seasonal = 6,     # panels in the seasonal (day-of-year) figure
  track_n_daily    = 4,     # panels in the daily-rhythm (hour-of-day) figure

  # Floor on a trajectory panel's half-width.
  track_min_halfwidth_m = 1500,

  # Buildings are drawn as individual footprints up to this half-width; wider frames use
  # the density raster instead.
  building_context_max_halfwidth_m = 4000
)

# Freezes GPS timestamps at this date (applied once, in 01_download.R) so span, panel-B
# membership and every derived count don't drift with the calendar on re-run. NULL takes
# every fix up to run time.
DATA_CUTOFF <- as.Date("2026-08-25")

RANDOM_SEED <- 42

# ---------------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------------

PATHS <- list(
  # Gitignored: large, regenerable.
  raw_dir     = file.path("data", "raw"),
  interim_dir = file.path("data", "interim"),
  osm_dir     = file.path("data", "osm"),
  ocsge_dir   = file.path("data", "ocsge"),

  # Committed: small derived artifacts and figures.
  processed_dir = file.path("data", "processed"),
  figures_dir   = "figures",

  # Named outputs.
  raw_tracks        = file.path("data", "raw", "movebank_tracks.rds"),
  deployments       = file.path("data", "raw", "movebank_deployments.rds"),
  clean_tracks      = file.path("data", "interim", "tracks_clean.rds"),
  settled_nights    = file.path("data", "interim", "settled_nights.rds"),
  qc_summary        = file.path("data", "processed", "qc_summary.csv"),
  individuals       = file.path("data", "processed", "individuals.csv"),
  covariate_stack   = file.path("data", "osm", "distance_stack.tif"),
  covariate_correlation = file.path("data", "processed", "covariate_correlation.csv"),
  ssf_coefficients  = file.path("data", "processed", "ssf_coefficients.csv"),
  ssf_individual    = file.path("data", "processed", "ssf_individual.csv"),
  ssf_models        = file.path("data", "interim", "ssf_models.rds")
)

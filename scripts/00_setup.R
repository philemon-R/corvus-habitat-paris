## 00_setup.R -- install the packages the pipeline needs.
##
## Run once, interactively, before script 01. Installs only what is missing.
##
##   Rscript scripts/00_setup.R

# When run non-interactively (Rscript), no CRAN mirror is set by default and
# install.packages() fails outright. Set one unless the user already has a preference.
repos <- getOption("repos")
if (is.null(repos[["CRAN"]]) || repos[["CRAN"]] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

PACKAGES <- c(
  # Data access
  "move2",        # Movebank download
  "keyring",      # stores Movebank credentials in the OS keychain
  # Spatial
  "sf",           # vector geometry
  "terra",        # raster / density surfaces
  "osmextract",   # OpenStreetMap via bulk regional extracts (Geofabrik)
  "archive",      # extracts OCS GE's .7z distribution
  "units",
  # Movement analysis
  "amt",          # step-selection functions
  "suncalc",      # sunrise/sunset, to isolate nocturnal fixes
  "purrr",
  # Modeling
  "survival",     # conditional logistic regression underlying the SSF
  "metafor",      # random-effects pooling of per-bird coefficients (utils.R, meta_pool())
  # Data handling and output
  "dplyr", "tidyr", "readr", "lubridate", "ggplot2",
  "patchwork",    # composes plots that cannot share a scale (see 06_maps.R)
  "ggspatial",    # scale bars on the map figures
  "maps",         # coarse national outline for the dispersal map
  "ragg"          # PNG device save_fig() uses (plotting.R)
)

missing <- PACKAGES[!vapply(PACKAGES, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) == 0) {
  message("All required packages are already installed.")
} else {
  message("Installing ", length(missing), " package(s): ", paste(missing, collapse = ", "))
  install.packages(missing)
}

## Movebank credentials -- run this ONCE, interactively, then never again:
##
##   move2::movebank_store_credentials("your_movebank_username")
##
## It prompts for the password and stores it in the OS keychain via `keyring`. The
## password is never written to this repository and never passed on a command line.
## Script 01 picks it up automatically.

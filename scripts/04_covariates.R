## 04_covariates.R -- build the habitat covariate layers, from OSM and OCS GE.
##
## Inputs:  config.R (OSM$bbox, OSM$features, OCSGE$classes), the cached extract from
##          script 03
## Outputs: data/osm/osm_layers.rds          -- cleaned sf layers, Lambert-93
##          data/osm/building_density_25m.tif -- built-up density, the first covariate
##          data/osm/distance_stack.tif      -- distance/density rasters from OSM and OCS GE
##
## Everything here is offline: it needs only config.R and the cached extracts from script
## 03. The stack takes a few minutes to build and is cached afterwards.
##
##   Rscript scripts/04_covariates.R

suppressPackageStartupMessages({
  library(sf); library(terra); library(dplyr); library(readr)
})

source("config.R")
source(file.path("scripts", "utils.R"))
source(file.path("scripts", "osm.R"))
source(file.path("scripts", "ocsge.R"))

ensure_dir(PATHS$osm_dir)
use_osm_cache()

# ---------------------------------------------------------------------------------
# 1. Extract the feature layers
# ---------------------------------------------------------------------------------

layers <- list()

# Buildings excluded: the covariate is a density, built below.
for (nm in setdiff(names(OSM$features), "building")) {
  message("Extracting '", nm, "' ...")
  layers[[nm]] <- osm_layer(nm)
  message("  ", nrow(layers[[nm]]), " features")
}

saveRDS(layers, file.path(PATHS$osm_dir, "osm_layers.rds"))
message("Saved ", length(layers), " layers to ", file.path(PATHS$osm_dir, "osm_layers.rds"))

# ---------------------------------------------------------------------------------
# 2. Built-up density
# ---------------------------------------------------------------------------------
# Fraction of each 25 m cell covered by a building footprint.

built <- building_density()
message("Built-up density: ", terra::ncell(built), " cells at ",
        terra::res(built)[1], " m, mean cover ",
        sprintf("%.3f", terra::global(built, "mean", na.rm = TRUE)[1, 1]))

# ---------------------------------------------------------------------------------
# 3. Distance to the nearest feature of each class
# ---------------------------------------------------------------------------------
# Distance computed on rasterised presence and not a unioned polygon layer.
# This quantisation loses information but saves computational cost.

STACK_PATH   <- file.path(PATHS$osm_dir, "distance_stack.tif")
DIST_CLASSES <- setdiff(names(OSM$features), "building")
FOCAL_NAME   <- sprintf("built_density_%dm", COVARIATES$focal_density_radius_m)

#' The stack layer name for one OCSGE$classes entry, matching what SSF$covariates names it.
ocsge_layer_name <- function(class_name, spec) {
  if (spec$type == "distance") paste0("dist_", class_name)
  else paste0(class_name, "_density_", COVARIATES$focal_density_radius_m, "m")
}
OCSGE_LAYER_NAMES <- mapply(ocsge_layer_name, names(OCSGE$classes), OCSGE$classes)

LAYER_NAMES  <- c(paste0("dist_", DIST_CLASSES), FOCAL_NAME, unname(OCSGE_LAYER_NAMES))

crs_wkt <- st_crs(CRS$projected)$wkt

#' Distance to the nearest cell holding any feature of one class.
distance_to_class <- function(name, grid) {
  v <- terra::project(
    terra::vect(osm_gpkg(), query = build_query(OSM$features[[name]], geometry_only = TRUE)),
    crs_wkt)
  # No `background`, so absent cells stay NA: terra::distance() measures from NA cells to
  # the nearest non-NA one, and returns 0 inside the features themselves.
  present <- terra::rasterize(v, grid, field = 1, touches = TRUE)
  list(n = nrow(v), d = terra::distance(present))
}

#' Distance-to-nearest or focal density for one OCSGE$classes entry. Unlike OSM, OCS GE
#' ships natively in Lambert-93, so no terra::project() step is needed, and a class is one
#' flat WHERE clause rather than osm_layer()'s tag-key/value rules.
ocsge_class_layer <- function(spec, grid) {
  q <- sprintf("SELECT * FROM %s WHERE %s", ocsge_schema()$layer, spec$where)
  v <- terra::vect(ocsge_gpkg(), query = q)
  if (spec$type == "distance") {
    present <- terra::rasterize(v, grid, field = 1, touches = TRUE)
    r <- terra::distance(present)
  } else {
    presence <- terra::rasterize(v, grid, cover = TRUE, background = 0)
    r <- terra::focal(presence, w = terra::focalMat(presence, COVARIATES$focal_density_radius_m, "circle"),
                      fun = "sum", na.rm = TRUE)
  }
  list(n = nrow(v), r = r)
}

# Cached like the extract itself: extent and layer names are both compared, so adding or
# removing a class rebuilds it. Doesn't catch a redefinition of an existing class: delete
# the file by hand for that.
stack_is_current <- function() {
  if (!file.exists(STACK_PATH)) return(FALSE)
  r <- try(terra::rast(STACK_PATH), silent = TRUE)
  if (inherits(r, "try-error")) return(FALSE)
  identical(names(r), LAYER_NAMES) &&
    max(abs(as.vector(terra::ext(r)) - as.vector(terra::ext(built)))) < COVARIATES$density_res_m
}

if (stack_is_current()) {
  covariates <- terra::rast(STACK_PATH)
  message("Covariate stack already matches config.R; reusing ", STACK_PATH)
} else {
  message("Building the covariate stack; cached afterwards.")

  # Reuse the grid of the building density
  grid <- terra::rast(built)

  layers <- list()
  for (nm in DIST_CLASSES) {
    t0 <- Sys.time()
    res <- distance_to_class(nm, grid)
    layers[[paste0("dist_", nm)]] <- res$d
    message(sprintf("  dist_%-9s %7s features, %5.0f s",
                    nm, format(res$n, big.mark = " "),
                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }

  # Compute built-up density at 200 m.
  t0 <- Sys.time()
  layers[[FOCAL_NAME]] <- terra::focal(
    built, w = terra::focalMat(built, COVARIATES$focal_density_radius_m, "circle"),
    fun = "sum", na.rm = TRUE)
  message(sprintf("  %-20s (focal mean, %d m), %5.0f s", FOCAL_NAME,
                  COVARIATES$focal_density_radius_m,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  # OCS GE classes, mirroring the OSM loop above.
  for (nm in names(OCSGE$classes)) {
    t0 <- Sys.time()
    res <- ocsge_class_layer(OCSGE$classes[[nm]], grid)
    layer_name <- ocsge_layer_name(nm, OCSGE$classes[[nm]])
    layers[[layer_name]] <- res$r
    message(sprintf("  %-20s %7s features, %5.0f s", layer_name,
                    format(res$n, big.mark = " "),
                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }

  covariates <- c(terra::rast(layers[paste0("dist_", DIST_CLASSES)]),
                  layers[[FOCAL_NAME]],
                  terra::rast(layers[unname(OCSGE_LAYER_NAMES)]))
  names(covariates) <- LAYER_NAMES

  terra::writeRaster(covariates, STACK_PATH, overwrite = TRUE,
                     datatype = "FLT4S",
                     gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES"))
  covariates <- terra::rast(STACK_PATH)
  message(sprintf("Wrote %s (%d layers, %.0f MB)", STACK_PATH, terra::nlyr(covariates),
                  file.size(STACK_PATH) / 1e6))
}

# ---------------------------------------------------------------------------------
# 4. Collinearity
# ---------------------------------------------------------------------------------
# Computed on each covariate's final form (-log1p(distance), i.e. proximity, matching
# prepare_covariates() in script 07).

# Cells drawn by number and read with extract(), not terra::spatSample(): same result
# but faster. Subset to SSF$covariates explicitly SSF$covariates is left out.
set.seed(RANDOM_SEED)
samp <- terra::extract(covariates, sample.int(terra::ncell(covariates), 200000))
samp <- samp[, SSF$covariates, drop = FALSE]
samp <- samp[complete.cases(samp), , drop = FALSE]
message(sprintf("\nCollinearity on %s randomly sampled cells, final (modeled) form:",
                format(nrow(samp), big.mark = " ")))

is_dist <- startsWith(names(samp), "dist_")
final   <- samp
final[is_dist] <- -log1p(samp[is_dist])
names(final)[is_dist] <- paste0("prox_", sub("^dist_", "", names(samp)[is_dist]))

cor_final <- round(stats::cor(final), 2)
print(cor_final)

# Tabular form, so the file stays readable when a class is added and the matrix changes shape.
idx <- which(upper.tri(cor_final), arr.ind = TRUE)
cor_long <- data.frame(a = rownames(cor_final)[idx[, "row"]],
                       b = colnames(cor_final)[idx[, "col"]],
                       r = cor_final[idx])
readr::write_csv(cor_long, PATHS$covariate_correlation)

worst <- slice_max(cor_long, abs(r), n = 3)
message("\nStrongest pairs:")
for (i in seq_len(nrow(worst))) {
  message(sprintf("  %-28s %-28s r = %+.2f", worst$a[i], worst$b[i], worst$r[i]))
}
message("\nCorrelation table: ", PATHS$covariate_correlation)


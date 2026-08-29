## osm.R -- shared OpenStreetMap extraction, used by the map and covariate scripts.
##
## Sourced, not run. Assumes 03_geodata_download.R has filled data/osm/ with a converted
## extract, clipped to the study area; every call here is local and offline.

suppressPackageStartupMessages({
  library(osmextract)
  library(sf)
  library(dplyr)
})

OSM_PLACE <- "Ile-de-France"

#' The study bounding box as an sf geometry, in geographic coordinates.
study_boundary <- function() {
  st_as_sfc(st_bbox(
    c(xmin = OSM$bbox[["xmin"]], ymin = OSM$bbox[["ymin"]],
      xmax = OSM$bbox[["xmax"]], ymax = OSM$bbox[["ymax"]]),
    crs = st_crs(CRS$geographic)))
}

# ---------------------------------------------------------------------------------
# The cached extract and its manifest
# ---------------------------------------------------------------------------------
#
# The GeoPackage script 03 writes is clipped to OSM$bbox, and nothing inside the file says
# which box that was. The manifest is one line written at conversion time, checked before any
# read.

MANIFEST_PATH <- function() file.path(PATHS$osm_dir, "extract_manifest.csv")

write_extract_manifest <- function(gpkg_path, pbf_md5 = NA_character_) {
  utils::write.csv(
    data.frame(gpkg = basename(gpkg_path),
               xmin = OSM$bbox[["xmin"]], ymin = OSM$bbox[["ymin"]],
               xmax = OSM$bbox[["xmax"]], ymax = OSM$bbox[["ymax"]],
               source_md5 = pbf_md5,
               converted_at = format(Sys.time(), tz = "UTC", usetz = TRUE)),
    MANIFEST_PATH(), row.names = FALSE)
  invisible(MANIFEST_PATH())
}

#' Path to the cached, study-area GeoPackage, or a clear error explaining what to run.
osm_gpkg <- function() {
  gpkg <- list.files(PATHS$osm_dir, pattern = "[.]gpkg$", full.names = TRUE)
  if (length(gpkg) == 0) {
    stop("No converted OSM extract in ", PATHS$osm_dir,
         ". Run scripts/03_geodata_download.R first (316 MB download).", call. = FALSE)
  }
  if (!file.exists(MANIFEST_PATH())) {
    stop("The OSM extract in ", PATHS$osm_dir, " has no manifest, so there is no way to ",
         "tell which study area it covers. Re-run scripts/03_geodata_download.R; the download ",
         "is cached, only the conversion repeats.", call. = FALSE)
  }
  m <- utils::read.csv(MANIFEST_PATH())
  want <- c(OSM$bbox[["xmin"]], OSM$bbox[["ymin"]], OSM$bbox[["xmax"]], OSM$bbox[["ymax"]])
  have <- c(m$xmin[1], m$ymin[1], m$xmax[1], m$ymax[1])
  if (max(abs(want - have)) > 1e-9) {
    stop("The cached OSM extract covers ", paste(sprintf("%.3f", have), collapse = ", "),
         " but OSM$bbox is now ", paste(sprintf("%.3f", want), collapse = ", "),
         ".\nRe-run scripts/03_geodata_download.R; the download is cached, only the ",
         "conversion repeats.", call. = FALSE)
  }
  gpkg[1]
}

#' Check the cache before a script starts reading from it.
use_osm_cache <- function() {
  Sys.setenv(OSMEXT_DOWNLOAD_DIRECTORY = normalizePath(PATHS$osm_dir, mustWork = TRUE))
  invisible(osm_gpkg())
}

# ---------------------------------------------------------------------------------
# Reading feature classes
# ---------------------------------------------------------------------------------

#' Turn a config.R OSM$features entry into an OGR SQL query. A feature class is a list of
#' rules, each naming one key and the values that count, ORed together (woodland has no
#' single OSM tag, so it's `natural` = wood OR `landuse` = forest). Keys are
#' double-quoted; NULL means "any non-null value of this key".
#'
#' `geometry` must appear in the SELECT list explicitly, or OGR SQL silently returns a
#' plain data frame that only fails on the first spatial operation.
build_query <- function(spec, geometry_only = FALSE) {
  clause <- function(rule) {
    key <- sprintf('"%s"', rule$key)
    if (is.null(rule$value)) {
      sprintf("%s IS NOT NULL", key)
    } else {
      sprintf("%s IN (%s)", key, paste(sprintf("'%s'", rule$value), collapse = ", "))
    }
  }
  where <- paste(vapply(spec, clause, character(1)), collapse = " OR ")

  # One informational column saying which tag matched. COALESCE over the keys of the class,
  # which can keep the wrong one for a feature carrying two of them (a park that is also
  # tagged landuse=grass), but nothing downstream keys on it.
  keys <- unique(vapply(spec, function(rule) rule$key, character(1)))
  tag  <- if (length(keys) == 1) sprintf('"%s"', keys) else
    sprintf("COALESCE(%s)", paste(sprintf('"%s"', keys), collapse = ", "))

  cols <- if (geometry_only) "geometry" else sprintf("osm_id, %s AS tag, geometry", tag)
  sprintf("SELECT %s FROM multipolygons WHERE %s", cols, where)
}

#' Every key a feature class queries, for the schema check in script 03.
feature_keys <- function(features = OSM$features) {
  unique(unlist(lapply(features, function(spec)
    vapply(spec, function(rule) rule$key, character(1)))))
}

#' Read one feature class from the cached extract, optionally cropped, and projected.
#'
#' Reads the cached GeoPackage directly rather than calling `oe_get()`: a `boundary`
#' argument to `oe_get()` re-converts the whole extract on every call. The clip happens
#' once in script 03; reads go straight at the GeoPackage with a spatial filter, answered
#' by the spatial index in seconds.
#'
#' OSM polygons overlap, self-intersect and are frequently mapped in several pieces, so
#' geometry is repaired on the way out (callers needing "distance to nearest X" should
#' union within the class too). Only broken geometry is repaired, not every geometry:
#' checking first with `st_is_valid()` and repairing only the handful that fail saves
#' computational cost.
#'
#' @param name A name in OSM$features.
#' @param box Optional `st_bbox` in the PROJECTED CRS, to read only what one panel needs.
#' @return An sf data frame in the projected CRS.
osm_layer <- function(name, box = NULL) {
  spec <- OSM$features[[name]]
  if (is.null(spec)) stop("No such OSM feature class: ", name, call. = FALSE)

  # The extract is in geographic coordinates; a projected frame has to come back to them
  # before filtering.
  # character(0), not "": sf hands the string to the WKT parser unconditionally, and an
  # empty string returns a parse error.
  filter_wkt <- character(0)
  if (!is.null(box)) {
    filter_wkt <- st_as_text(st_transform(st_as_sfc(box), CRS$geographic))
  }

  x <- st_read(osm_gpkg(), query = build_query(spec), wkt_filter = filter_wkt, quiet = TRUE)

  broken <- !st_is_valid(x)
  broken[is.na(broken)] <- TRUE
  if (any(broken)) x[broken, ] <- st_make_valid(x[broken, ])

  x |>
    filter(!st_is_empty(geometry)) |>
    st_transform(CRS$projected)
}

# ---------------------------------------------------------------------------------
# Building density
# ---------------------------------------------------------------------------------

#' Building footprint density over the study area, as a cached raster. Two uses:
#' as a basemap it's the only workable way to draw the urban fabric on a wide frame
#' (the 2.4M-polygon footprint layer flattens into a gray wash past a few km, a raster
#' draws in constant time); as a covariate it's the built-up density script 07 reads at
#' every step endpoint, which is why it lives here and why its resolution follows the
#' modeling need (COVARIATES$density_res_m).
#'
#' `cover = TRUE` returns the fraction of each cell covered, not a presence flag.
#'
#' @param res_m Cell size in meters.
#' @param rebuild Rebuild even if the cached file exists.
#' @return A terra SpatRaster of footprint area fraction, 0 to 1, in the projected CRS.
building_density <- function(res_m = COVARIATES$density_res_m, rebuild = FALSE) {
  path <- file.path(PATHS$osm_dir, sprintf("building_density_%dm.tif", res_m))
  bb   <- st_bbox(st_transform(study_boundary(), CRS$projected))

  # Staleness trap: compare the extents of the raster and extract.
  if (!rebuild && file.exists(path)) {
    cached <- terra::rast(path)
    want <- c(bb[["xmin"]], bb[["xmax"]], bb[["ymin"]], bb[["ymax"]])
    if (max(abs(as.vector(terra::ext(cached)) - want)) < res_m) return(cached)
    message("Cached building density does not match OSM$bbox; rebuilding.")
  }

  message("Building the ", res_m, " m building-density raster. First run only.")

  # Read with terra, not osm_layer(): terra::vect() to avoid sf polygon rebuild.
  #
  # CRS as WKT, not "EPSG:2154": on some Windows/PostGIS setups, a stale PROJ database
  # shadows the real one and the EPSG lookup fails with "empty srs". WKT needs no lookup.
  crs_wkt <- st_crs(CRS$projected)$wkt
  footprints <- terra::project(
    terra::vect(osm_gpkg(), query = build_query(OSM$features$building, geometry_only = TRUE)),
    crs_wkt)
  message("  ", format(nrow(footprints), big.mark = " "), " footprints")

  grid <- terra::rast(terra::ext(bb[["xmin"]], bb[["xmax"]], bb[["ymin"]], bb[["ymax"]]),
                      resolution = res_m, crs = crs_wkt)

  # cover = TRUE, not rasterizeGeom(fun = "area"): the two aren't the same quantity.
  # rasterizeGeom sums each footprint's intersection area separately; cover
  # returns the fraction covered by the UNION, bounded at 1, faster, and tolerant of
  # OSM's self-intersecting geometry.
  r <- terra::rasterize(footprints, grid, cover = TRUE, background = 0,
                        filename = path, overwrite = TRUE,
                        wopt = list(names = "built_fraction",
                                    gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3")))
  message("  wrote ", path, " (", round(file.size(path) / 1e6), " MB)")
  r
}

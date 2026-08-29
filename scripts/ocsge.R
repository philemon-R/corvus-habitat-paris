## ocsge.R -- shared IGN OCS GE reading, used alongside osm.R by the covariate scripts.
##
## Sourced, not run. Assumes 03_geodata_download.R has filled data/ocsge/ with a converted,
## clipped GeoPackage; every call here is local and offline. Mirrors osm.R's shape: same
## manifest-staleness gate, same wkt_filter convention, same validity repair on the way out.

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
})

# ---------------------------------------------------------------------------------
# The cached extract and its manifest
# ---------------------------------------------------------------------------------
#
# As in osm.R: the GeoPackage is clipped to OSM$bbox with no record of that inside the
# file. Also records which layer/field script 03 found the couverture code in, since that
# schema wasn't knowable until a real file could be extracted.

OCSGE_MANIFEST_PATH <- function() file.path(PATHS$ocsge_dir, "ocsge_manifest.csv")

write_ocsge_manifest <- function(gpkg_path, layer, code_field, departments) {
  utils::write.csv(
    data.frame(gpkg = basename(gpkg_path),
               xmin = OSM$bbox[["xmin"]], ymin = OSM$bbox[["ymin"]],
               xmax = OSM$bbox[["xmax"]], ymax = OSM$bbox[["ymax"]],
               millesime = OCSGE$millesime,
               departments = paste(departments, collapse = ";"),
               layer = layer, code_field = code_field,
               converted_at = format(Sys.time(), tz = "UTC", usetz = TRUE)),
    OCSGE_MANIFEST_PATH(), row.names = FALSE)
  invisible(OCSGE_MANIFEST_PATH())
}

#' Path to the cached, study-area OCS GE GeoPackage, or a clear error explaining what to run.
ocsge_gpkg <- function() {
  gpkg <- file.path(PATHS$ocsge_dir, "ocsge_couverture.gpkg")
  if (!file.exists(gpkg)) {
    stop("No converted OCS GE extract in ", PATHS$ocsge_dir,
         ". Run scripts/03_geodata_download.R first.", call. = FALSE)
  }
  if (!file.exists(OCSGE_MANIFEST_PATH())) {
    stop("The OCS GE extract in ", PATHS$ocsge_dir, " has no manifest, so there is no way ",
         "to tell which study area it covers. Re-run scripts/03_geodata_download.R.",
         call. = FALSE)
  }
  m <- utils::read.csv(OCSGE_MANIFEST_PATH())
  want <- c(OSM$bbox[["xmin"]], OSM$bbox[["ymin"]], OSM$bbox[["xmax"]], OSM$bbox[["ymax"]])
  have <- c(m$xmin[1], m$ymin[1], m$xmax[1], m$ymax[1])
  if (max(abs(want - have)) > 1e-9) {
    stop("The cached OCS GE extract covers ", paste(sprintf("%.3f", have), collapse = ", "),
         " but OSM$bbox is now ", paste(sprintf("%.3f", want), collapse = ", "),
         ".\nRe-run scripts/03_geodata_download.R.", call. = FALSE)
  }
  gpkg
}

#' The layer and code-field name script 03 discovered in the OCS GE schema, from the manifest.
ocsge_schema <- function() {
  ocsge_gpkg()   # trigger the staleness check first
  m <- utils::read.csv(OCSGE_MANIFEST_PATH())
  list(layer = m$layer[1], code_field = m$code_field[1])
}

# ---------------------------------------------------------------------------------
# Reading couverture classes
# ---------------------------------------------------------------------------------

#' Read one couverture/usage class from the cached, clipped OCS GE extract. Takes a raw SQL
#' WHERE clause, not a codes-in-list, since a flat equality/IN filter can't express
#' `nonagri_herb` (which also excludes an agricultural usage code).
#'
#' @param where A SQL WHERE clause against the couverture layer, e.g.
#'   `"code_cs IN ('CS2.1.1.1','CS2.1.1.2')"` or one of `OCSGE$classes[[x]]$where`.
#' @param box Optional `st_bbox` in the PROJECTED CRS, to read only what one panel needs.
#' @return An sf data frame in the projected CRS.
ocsge_layer <- function(where, box = NULL) {
  schema <- ocsge_schema()

  # OCS GE ships natively in Lambert-93, unlike the OSM extract, so `box` needs no
  # transform here (unlike osm_layer()'s), since it's already in this layer's own CRS.
  filter_wkt <- character(0)
  if (!is.null(box)) {
    filter_wkt <- st_as_text(st_as_sfc(box))
  }

  query <- sprintf("SELECT * FROM %s WHERE %s", schema$layer, where)
  x <- st_read(ocsge_gpkg(), query = query, wkt_filter = filter_wkt, quiet = TRUE)

  broken <- !st_is_valid(x)
  broken[is.na(broken)] <- TRUE
  if (any(broken)) x[broken, ] <- st_make_valid(x[broken, ])

  # Not filter(!st_is_empty(geometry)): OCS GE's geometry column isn't named "geometry".
  # st_geometry() finds it by its sf attribute instead.
  x[!st_is_empty(st_geometry(x)), ] |>
    st_transform(CRS$projected)
}

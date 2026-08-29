## 03_geodata_download.R -- fetch and prepare both land-cover sources: the Geofabrik
## OpenStreetMap extract and IGN's OCS GE.
##
## Inputs:  none beyond config.R (OSM$version, OCSGE$millesime, OCSGE$departments)
## Outputs: data/osm/geofabrik_ile-de-france-*.osm.pbf   -- the raw OSM extract
##          data/osm/geofabrik_ile-de-france-*.gpkg      -- converted, clipped to OSM$bbox
##          data/osm/extract_manifest.csv               -- which box that GeoPackage covers
##          data/processed/osm_provenance.csv           -- which OSM snapshot, committed
##          data/ocsge/*.7z                    -- one raw OCS GE archive per department
##          data/ocsge/*.gpkg                  -- each archive's extracted GeoPackage
##          data/ocsge/ocsge_couverture.gpkg   -- all departments, clipped and merged
##          data/ocsge/ocsge_manifest.csv      -- which box/schema that GeoPackage covers
##
## Nothing here depends on the tracking data, so it can run at any point before scripts 04
## or 05.
##
##   Rscript scripts/03_geodata_download.R
##
## OpenStreetMap data is ODbL; IGN's OCS GE is Etalab-licensed (open). Both are regenerated
## by this script rather than redistributed: data/osm/ and data/ocsge/ are gitignored.

suppressPackageStartupMessages({
  library(osmextract)
  library(sf)
  library(archive)
})

source("config.R")
source(file.path("scripts", "utils.R"))
source(file.path("scripts", "osm.R"))     # OSM_PLACE, study_boundary(), the OSM manifest
source(file.path("scripts", "ocsge.R"))   # write_ocsge_manifest()

ensure_dir(PATHS$osm_dir)
ensure_dir(PATHS$ocsge_dir)

# Cache inside the repo's gitignored data/osm/, not osmextract's default per-user app-data
# directory, so a reader can see the extract, check its date, and delete it. Must be set
# before any oe_* call; script 04 sets the same variable to find what this script cached.
Sys.setenv(OSMEXT_DOWNLOAD_DIRECTORY = normalizePath(PATHS$osm_dir, mustWork = TRUE))

# ===================================================================================
# 1. OpenStreetMap (Geofabrik)
# ===================================================================================

# ---------------------------------------------------------------------------------
# 1.1. Resolve the extract
# ---------------------------------------------------------------------------------

match <- oe_match(OSM_PLACE, version = OSM$version, quiet = TRUE)
message("Geofabrik extract for '", OSM_PLACE, "', version '", OSM$version, "'")
message("  url: ", match$url)

# Resolved before the download so it's on record even if the download fails. "latest" is
# rebuilt nightly, so without the server's Last-Modified date there would be no way to say
# afterward which state of the map was used.
published_md5 <- tryCatch(
  sub("\\s.*$", "", readLines(paste0(match$url, ".md5"), warn = FALSE)[1]),
  error = function(e) NA_character_)
last_modified <- tryCatch({
  h <- curlGetHeaders(match$url)
  trimws(sub("^[^:]+:", "", grep("^Last-Modified:", h, value = TRUE, ignore.case = TRUE)[1]))
}, error = function(e) NA_character_)

message("  snapshot (Last-Modified): ", last_modified)
message("  published md5: ", published_md5)
if (!is.na(match$file_size)) {
  message("  size: ", round(match$file_size / 1024^2), " MB")
}

# Recorded before oe_download() to distinguish whether a checksum mismatch
# means a corrupt transfer or just a new "latest" version.
was_cached <- length(list.files(PATHS$osm_dir, pattern = "[.]osm[.]pbf$")) > 0

pbf_path <- oe_download(
  file_url         = match$url,
  file_size        = match$file_size,
  download_directory = PATHS$osm_dir,
  max_file_size    = 6e8,
  quiet            = FALSE
)
message("Extract at: ", pbf_path)

# Verified against Geofabrik's published md5: a truncated download is otherwise silent
# until GDAL fails mid-conversion, or converts a partial file without complaint..
local_md5 <- unname(tools::md5sum(pbf_path))
md5_matches <- !is.na(published_md5) && identical(local_md5, published_md5)

if (!md5_matches && !was_cached) {
  stop("Checksum mismatch on a file just downloaded: ", basename(pbf_path),
       "\n  published: ", published_md5,
       "\n  local:     ", local_md5,
       "\nThe transfer was corrupt. Delete the file and re-run.", call. = FALSE)
}
if (!md5_matches) {
  message("  NOTE: the cached extract is not the latest Geofabrik snapshot.")
  message("        published ", published_md5, ", local ", local_md5)
  message("        'latest' is rebuilt nightly. Using the cached file, which is the one")
  message("        recorded in osm_provenance.csv. Delete it to fetch the current build,")
  message("        or set OSM$version to a dated snapshot to stop the question arising.")
} else {
  message("  checksum verified")
}

# Written to data/processed/ (committed, unlike data/osm/): the record of which OSM
# snapshot produced a given result. To reproduce a past run, set OSM$version to this date.
ensure_dir(PATHS$processed_dir)
write.csv(
  data.frame(
    place         = OSM_PLACE,
    version       = OSM$version,
    url           = match$url,
    md5           = local_md5,
    published_md5 = published_md5,
    is_current    = md5_matches,
    snapshot      = if (md5_matches) last_modified else
                      format(file.mtime(pbf_path), tz = "UTC", usetz = TRUE),
    size_bytes    = file.size(pbf_path),
    recorded_at   = format(Sys.time(), tz = "UTC", usetz = TRUE)
  ),
  file.path(PATHS$processed_dir, "osm_provenance.csv"), row.names = FALSE)
message("  provenance: ", file.path(PATHS$processed_dir, "osm_provenance.csv"))

# ---------------------------------------------------------------------------------
# 1.2. Convert to GeoPackage
# ---------------------------------------------------------------------------------
# GDAL cannot query a .pbf efficiently, so osmextract converts it once to a GeoPackage.
#
# The conversion is clipped to OSM$bbox. Passing `boundary` to oe_get() also works, but
# forces a full re-conversion on every call. The box used is recorded in extract_manifest.csv.
# scripts/osm.R refuses to read the extract once OSM$bbox has changed.

# Verify whether conversion is required
gpkg_existing <- list.files(PATHS$osm_dir, pattern = "[.]gpkg$", full.names = TRUE)
gpkg_current  <- tryCatch(!is.null(osm_gpkg()), error = function(e) FALSE)

if (gpkg_current) {
  gpkg_path <- gpkg_existing[1]
  message("GeoPackage already built for this study area: ", gpkg_path)
} else {
  message("Converting to GeoPackage, clipped to the study area. Fifteen minutes or so.")
  gpkg_path <- oe_vectortranslate(
    file_path     = pbf_path,
    layer         = "multipolygons",
    boundary      = study_boundary(),
    boundary_type = "spat",     # bbox intersection; "clipsrc" would cut geometries
    quiet         = FALSE
  )
  write_extract_manifest(gpkg_path, pbf_md5 = local_md5)
  message("GeoPackage at: ", gpkg_path)
  message("  manifest: ", MANIFEST_PATH())
}

# ---------------------------------------------------------------------------------
# 1.3. Verify the keys this pipeline queries are actually present
# ---------------------------------------------------------------------------------
# An OSM query against a missing column returns an empty layer, not an error.
# Check the columns exist here, where the failure is still legible.

needed <- feature_keys()
# No `layer` argument: GDAL ignores it when a query is given, and warns about it.
present <- names(sf::st_read(gpkg_path, quiet = TRUE,
                             query = "SELECT * FROM multipolygons LIMIT 0"))
missing <- setdiff(needed, present)

if (length(missing) > 0) {
  stop("Keys missing from the converted layer: ", paste(missing, collapse = ", "),
       "\nRe-run oe_vectortranslate() with extra_tags = c(",
       paste(sprintf('"%s"', missing), collapse = ", "), ")", call. = FALSE)
}
message("All queried keys present: ", paste(needed, collapse = ", "))

for (f in c(pbf_path, gpkg_path)) {
  message(sprintf("  %-60s %6.0f MB", basename(f), file.size(f) / 1024^2))
}
message("OpenStreetMap part done. Those layers can now be read offline.\n")

# ===================================================================================
# 2. IGN OCS GE
# ===================================================================================
# OCS GE ships per department, not as one regional file like Geofabrik's OSM extract, so
# this downloads several small archives rather than one large one (~440 MB total).

# ---------------------------------------------------------------------------------
# 2.1. Download each department's archive
# ---------------------------------------------------------------------------------
# Departments come from OCSGE$departments in config.R -- hardcoded, not re-derived from
# OSM$bbox here. If OSM$bbox changes, the list might need to be updated.
#
# Idempotent: skipped if already on disk. No API key -- the Géoplateforme download API is
# open.

message("Departments to fetch from OCSGE$departments in config.R: ",
       paste(OCSGE$departments, collapse = ", "))
message("Make sure this list still covers OSM$bbox before relying on the result.")

title_for <- function(dept) sprintf("OCS-GE_2-0__GPKG_LAMB93_D%s_%s", dept, OCSGE$millesime)
url_for   <- function(dept) {
  t <- title_for(dept)
  sprintf("https://data.geopf.fr/telechargement/download/OCSGE/%s/%s.7z", t, t)
}

archives <- character(0)
for (dept in OCSGE$departments) {
  archive_path <- file.path(PATHS$ocsge_dir, paste0(title_for(dept), ".7z"))
  if (file.exists(archive_path)) {
    message("D", dept, ": already downloaded (", round(file.size(archive_path) / 1e6), " MB)")
  } else {
    message("D", dept, ": downloading ", url_for(dept))
    utils::download.file(url_for(dept), archive_path, mode = "wb", quiet = FALSE)
  }
  archives <- c(archives, archive_path)
}

# ---------------------------------------------------------------------------------
# 2.2. Extract
# ---------------------------------------------------------------------------------
# The .7z bundles two GeoPackages (OCCUPATION_SOL.gpkg and ZONE_CONSTRUITE.gpkg), not one
# -- confirmed by extracting Paris's archive, so this doesn't assume a single output file.
# `archive::archive_extract()` handles the extraction directly, without depending on a
# system `tar`/`7z` binary supporting `.7z`'s LZMA codec, which isn't guaranteed to be
# present or capable on every platform.
#
# The couverture code lives in OCCUPATION_SOL.gpkg, layer OCCUPATION_SOL, field code_cs.

schema <- list(layer = "OCCUPATION_SOL", code_field = "code_cs")
gpkg_paths <- character(0)

for (archive_path in archives) {
  dept_dir <- file.path(PATHS$ocsge_dir, tools::file_path_sans_ext(basename(archive_path)))
  gpkg_path <- list.files(dept_dir, pattern = "^OCCUPATION_SOL[.]gpkg$",
                          full.names = TRUE, recursive = TRUE)

  if (length(gpkg_path) == 0) {
    ensure_dir(dept_dir)
    message("Extracting ", basename(archive_path))
    archive::archive_extract(archive_path, dir = dept_dir)
    gpkg_path <- list.files(dept_dir, pattern = "^OCCUPATION_SOL[.]gpkg$",
                            full.names = TRUE, recursive = TRUE)
  }
  if (length(gpkg_path) == 0) {
    stop("No OCCUPATION_SOL.gpkg found after extracting ", basename(archive_path), call. = FALSE)
  }
  gpkg_paths <- c(gpkg_paths, gpkg_path[1])
}
message("Couverture layer: ", schema$layer, ", code field: ", schema$code_field)

# ---------------------------------------------------------------------------------
# 2.3. Clip each department once, merge, cache
# ---------------------------------------------------------------------------------

message("Clipping ", length(gpkg_paths), " department(s) to the study area...")
# OCS GE ships natively in Lambert-93, unlike the OSM extract; study_boundary() returns
# geographic coordinates, so it must be transformed before use as a wkt_filter.
filter_wkt <- st_as_text(st_transform(study_boundary(), CRS$projected))

clips <- lapply(gpkg_paths, function(p) {
  x <- st_read(p, layer = schema$layer, wkt_filter = filter_wkt, quiet = TRUE)
  broken <- !st_is_valid(x)
  broken[is.na(broken)] <- TRUE
  if (any(broken)) x[broken, ] <- st_make_valid(x[broken, ])
  x
})
combined <- do.call(rbind, clips)
# Not dplyr::filter(!st_is_empty(geometry)): OCS GE's active geometry column isn't named
# "geometry". st_geometry() finds it by its sf attribute instead.
combined <- combined[!st_is_empty(st_geometry(combined)), ]
message("  ", format(nrow(combined), big.mark = " "), " couverture polygons after clipping")

out_gpkg <- file.path(PATHS$ocsge_dir, "ocsge_couverture.gpkg")
st_write(combined, out_gpkg, layer = schema$layer, delete_dsn = TRUE, quiet = TRUE)

# ---------------------------------------------------------------------------------
# 2.4. Manifest
# ---------------------------------------------------------------------------------

write_ocsge_manifest(out_gpkg, layer = schema$layer, code_field = schema$code_field,
                     departments = OCSGE$departments)
message("  manifest: ", OCSGE_MANIFEST_PATH())

for (f in c(archives, out_gpkg)) {
  message(sprintf("  %-60s %6.0f MB", basename(f), file.size(f) / 1024^2))
}
message("OCS GE part done. Those layers can now be read offline.")
message("\nBoth parts done. scripts/04_covariates.R and 06_maps.R can now run offline.")

## plotting.R -- shared figure style and output helper.
##
## Sourced by any script that writes a figure.

suppressPackageStartupMessages(library(ggplot2))

#' Project theme: light, readable at README width.
theme_corvus <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title      = element_text(face = "bold", size = rel(1.1)),
      plot.subtitle   = element_text(color = "gray30", margin = margin(b = 10)),
      plot.caption    = element_text(color = "gray45", hjust = 0, size = rel(0.8),
                                     margin = margin(t = 10)),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      # Extra right margin: the outermost x-axis tick label overhangs the panel and gets
      # clipped otherwise.
      plot.margin = margin(t = 5.5, r = 18, b = 5.5, l = 5.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray92"),
      strip.text      = element_text(face = "bold", size = rel(0.9)),
      legend.position = "bottom",
      legend.title    = element_text(size = rel(0.9))
    )
}

# One stable color per cohort across every figure. Okabe-Ito, distinguishable under common
# forms of color blindness.
COHORT_LEVELS  <- 2019:2025
COHORT_COLORS <- setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442"),
  COHORT_LEVELS
)
GRAY_EXCLUDED <- "gray65"

# Palette for the maps
scale_season <- function(...) {
  ggplot2::scale_color_viridis_c(
    option = "plasma", end = 0.92,
    guide = ggplot2::guide_colorbar(draw.ulim = FALSE, draw.llim = FALSE), ...)
}

# CLOCK is cyclic: returns to its starting color. 
# Navy at midnight, violet at dawn, amber at midday, brick at dusk.
CLOCK_COLORS <- c("#2b2d5c", "#8c4bb0", "#f0a202", "#c9403a", "#2b2d5c")

#' Cohort as a factor with fixed levels, so colors do not shift between figures that
#' contain different cohorts.
cohort_factor <- function(x) factor(as.integer(x), levels = COHORT_LEVELS)

scale_cohort <- function(aesthetics = "color", name = "Cohort", ...) {
  ggplot2::scale_color_manual(aesthetics = aesthetics, name = name,
                               values = COHORT_COLORS, drop = TRUE, ...)
}

# Every figure is written at the same width, only height varies: Markdown scales an image
# to the column width regardless of its pixel size, so a fixed width keeps axis text a
# consistent size across the README.
FIG_WIDTH_CM <- 24

# 500 DPI for legibility of small map elements.
FIG_DPI      <- 500

cm_to_in <- function(x) x / 2.54

#' Write a figure to figures/ as PNG, and report where it went.
#'
#' @param plot A ggplot object.
#' @param name File name without extension, prefixed by the script number.
#' @param height_cm Figure height in centimeters. The width is fixed, see above.
#'
#' `ragg::agg_png()`, not the base PNG device: the base Windows device is
#' slow on many overlapping translucent segments and cairo produces artifacts.
save_fig <- function(plot, name, height_cm = 13) {
  ensure_dir(PATHS$figures_dir)
  path <- file.path(PATHS$figures_dir, paste0(name, ".png"))
  ggsave(path, plot, width = cm_to_in(FIG_WIDTH_CM), height = cm_to_in(height_cm),
         dpi = FIG_DPI, bg = "white", device = ragg::agg_png)
  message(sprintf("  figure: %-34s %2.0f x %2.0f cm", path, FIG_WIDTH_CM, height_cm))
  invisible(path)
}

# Land cover colors.
LAND_COLORS <- c(
  farmland = "#e4d6ab",
  green    = "#cfdfc7",
  wood     = "#8fb37e",
  water    = "#bcd6e8"
)

# Which OCS GE covariate class (config.R, OCSGE$classes) draws in which LAND_COLORS tone.
OCSGE_COVARIATE_COLORS <- c(
  woody        = unname(LAND_COLORS[["wood"]]),
  nonagri_herb = unname(LAND_COLORS[["green"]]),
  agricultural = unname(LAND_COLORS[["farmland"]])
)

# Building color for the fused basemap's overlay (coverage encoded as alpha).
BUILT_COLOR <- "#c6c3bd"

# Footprint fraction at which built_overlay() saturates to fully opaque..
BUILT_SATURATION <- 0.6

# ---------------------------------------------------------------------------------
# Basemap
# ---------------------------------------------------------------------------------
# The layers a map is drawn on, and the frame it is drawn in.

#' The building-density raster, read once per session.
built_raster <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- building_density()
    cache
  }
})

#' The three OCS GE covariate classes (config.R, OCSGE$classes), read once per session for
#' the whole study area as a named list of sf frames: several figures draw the same class
#' at more than one frame.
ocsge_context_layers <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      classes <- OCSGE$classes[names(OCSGE_COVARIATE_COLORS)]
      cache <<- lapply(classes, function(spec) ocsge_layer(spec$where))
    }
    cache
  }
})

#' Select the features of a layer that touch one frame.
#'
#' Faster to select features and clip them at drawing time than to cut them to the frame.
#'
#' @param layer An `sf` data frame or an `sfc`, in the projected CRS.
#' @param box An `st_bbox` in the projected CRS.
#' @return `layer`, keeping only the features whose bounding box meets `box`.
clip_to_box <- function(layer, box) {
  hit <- sf::st_intersects(sf::st_as_sfc(box), layer)[[1]]
  if (inherits(layer, "sfc")) layer[hit] else layer[hit, ]
}

#' Every basemap layer for one frame, in drawing order: OCS GE's three covariate classes
#' first, OSM buildings and water on top.
#'
#' Buildings switch between individual footprints and the density overlay at
#' MAP$building_context_max_halfwidth_m..
#'
#' @param box An `st_bbox` in the projected CRS.
#' @return A list of ggplot layers.
fused_base_map <- function(box) {
  half   <- unname((box[["xmax"]] - box[["xmin"]]) / 2)
  layers <- ocsge_context_layers()

  building_layer <- if (half <= MAP$building_context_max_halfwidth_m) {
    ggplot2::geom_sf(data = osm_layer("building", box), fill = "gray90", color = NA)
  } else {
    built_overlay(box)
  }

  c(
    lapply(names(layers), function(nm)
      ggplot2::geom_sf(data = clip_to_box(layers[[nm]], box),
                       fill = OCSGE_COVARIATE_COLORS[[nm]], color = NA)),
    list(building_layer),
    list(ggplot2::geom_sf(data = osm_layer("water", box),
                          fill = LAND_COLORS[["water"]], color = NA))
  )
}

#' Finalize a map: mark the capture site, add a scale bar, fix the frame, drop the graticule.
#'
#' `expand = FALSE` to enforce axis limits.
#'
#' @param p A ggplot object with its data layers already added.
#' @param box An `st_bbox` in the projected CRS.
#' @param trap One-row data frame of x and y; the capture site by default.
#' @return The ggplot, framed.
frame_map <- function(p, box, trap = trap_point()) {
  p +
    ggplot2::geom_point(data = trap, ggplot2::aes(x = x, y = y), shape = 4, size = 2.4,
                        stroke = 1, color = "gray15") +
    ggspatial::annotation_scale(location = "br", line_width = 0.5,
                                height = ggplot2::unit(0.12, "cm"),
                                text_cex = 0.7, text_col = "gray30",
                                bar_cols = c("gray30", "white")) +
    ggplot2::coord_sf(xlim = c(box[["xmin"]], box[["xmax"]]),
                      ylim = c(box[["ymin"]], box[["ymax"]]),
                      expand = FALSE, crs = sf::st_crs(CRS$projected), datum = NA) +
    theme_corvus() +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(fill = "white", color = NA))
}

#' A square frame of a given half-width around a point, in the projected CRS.
square_box <- function(x, y, half_m) {
  sf::st_bbox(c(xmin = x - half_m, xmax = x + half_m,
                ymin = y - half_m, ymax = y + half_m),
              crs = sf::st_crs(CRS$projected))
}

#' The building-density raster as a ggplot annotation, cropped to one frame -- an ALPHA
#' overlay, not a color ramp: transparent at zero density, opaque at BUILT_SATURATION.
#'
#' `annotation_raster`, not `geom_raster`: draws an image without claiming the fill scale
#' needed for the density maps.
#'
#' @param box An `st_bbox` in the projected CRS.
#' @param r A SpatRaster of footprint fraction; defaults to the cached one.
#' @param color The building color at full saturation; an RGB hex, no alpha of its own.
#' @return A ggplot layer, or NULL where the frame and the raster do not overlap.
built_overlay <- function(box, r = built_raster(), color = BUILT_COLOR) {
  e <- as.vector(terra::ext(r))
  xmin <- max(e[1], box[["xmin"]]); xmax <- min(e[2], box[["xmax"]])
  ymin <- max(e[3], box[["ymin"]]); ymax <- min(e[4], box[["ymax"]])
  if (xmax <= xmin || ymax <= ymin) return(NULL)

  cropped <- terra::crop(r, terra::ext(xmin, xmax, ymin, ymax))
  m <- terra::as.matrix(cropped, wide = TRUE)   # first row is the northern edge, as needed

  # Square root, not linear. Footprint fraction is heavily skewed: outside the dense core
  # most cells sit at a few percent, and on a linear ramp the whole suburban ring would be
  # indistinguishable from open ground.
  shade <- sqrt(pmin(pmax(m, 0), BUILT_SATURATION) / BUILT_SATURATION)
  alpha <- as.integer(round(pmin(pmax(shade, 0), 1) * 255))
  alpha[is.na(alpha)] <- 0L

  rgb <- grDevices::col2rgb(color)
  hex <- sprintf("#%02X%02X%02X", rgb[1], rgb[2], rgb[3])
  cols <- matrix(sprintf("%s%02X", hex, alpha), nrow = nrow(m), ncol = ncol(m))

  e <- terra::ext(cropped)
  ggplot2::annotation_raster(cols, xmin = e[1], xmax = e[2], ymin = e[3], ymax = e[4],
                             interpolate = TRUE)
}

#' A compact legend for what the fused basemap means: buildings, woody canopy,
#' herbaceous, agricultural, water. Drawn once and reused rather than
#' spelled out in every caption. Not a ggplot scale, since the basemap is drawn with
#' `annotation_raster()`/`geom_sf()` fills rather than a mapped aesthetic -- built as its
#' own tiny plot instead, same approach as `clock_legend()`.
#'
#' @return A ggplot object, meant to be stacked under a figure with patchwork's `/`.
fused_legend <- function() {
  discrete <- data.frame(
    x     = 2:5,
    label = c("Woody canopy", "Herbaceous", "Agricultural", "Water"),
    fill  = c(unname(LAND_COLORS[["wood"]]), unname(LAND_COLORS[["green"]]),
             unname(LAND_COLORS[["farmland"]]), unname(LAND_COLORS[["water"]]))
  )
  # Gradient for built-up density.
  grad <- data.frame(x = seq(1 - 0.35, 1 + 0.35, length.out = 24)) |>
    dplyr::mutate(fill = scales::gradient_n_pal(c("white", BUILT_COLOR))(
      seq(0, 1, length.out = 24)))
  labels <- data.frame(x = 1:5, label = c("Buildings", discrete$label))

  ggplot() +
    geom_tile(data = grad, aes(x = x, y = 1, fill = fill), width = 0.033, height = 0.6) +
    geom_tile(data = discrete, aes(x = x, y = 1, fill = fill), width = 0.6, height = 0.6,
              color = "gray65", linewidth = 0.2) +
    geom_text(data = labels, aes(x = x, y = 0.55, label = label),
              size = 2.5, color = "gray30", vjust = 1) +
    scale_fill_identity() +
    scale_x_continuous(limits = c(0.4, 5.6)) +
    scale_y_continuous(limits = c(0, 1.35)) +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(plot.margin = margin(2, 4, 0, 4))
}

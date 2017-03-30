# Figure 2: Rank priority maps for different methods (RWR, ZON and ILP) in
# columns, and for different datasets in rows (ALL, ES, BD).

library(dplyr)
library(grid)
library(lazyeval)
library(magick)
library(maptools)
library(raster)
library(RColorBrewer)
library(rgdal)
library(sp)
library(tmap)

data(Europe, land)

# Helper functions --------------------------------------------------------


get_map_params <- function(legend.reversed=FALSE) {
  # Define a suitable bounding box
  bbox <- matrix(c(2635899, 1386018, 6084606, 5307234),
                         nrow = 2, ncol = 2, dimnames = list(c("x", "y"),
                                                             c("min", "max")))

  breaks <- c(0, 0.2, 0.5, 0.75, 0.9, 0.95, 0.98, 1)
  colors <- rev(RColorBrewer::brewer.pal(length(breaks) - 1, "RdYlBu"))
  labels <- (100 - breaks * 100)
  labels <- cbind(labels[1:(length(labels) - 1)], labels[2:length(labels)])
  labels[,2] <- paste(labels[,2], "%")
  labels[7,2] <- ""
  labels <- apply(labels, 1, paste, collapse = " - ")
  labels[7] <- gsub(" - ", " %", labels[7])

  params <- list()

  params$bbox <- bbox

  if (legend.reversed) {
    params$breaks <- rev(breaks)
    params$colors <- rev(colors)
    params$labels <- rev(labels)
  } else {
    params$breaks <- breaks
    params$colors <- colors
    params$labels <- labels
  }
  return(params)
}

create_raster_levels <- function(raster) {

  params <- get_map_params()

  # Create a RasterLayer with a RAT
  rat_raster <- raster::ratify(raster)
  rat <- levels(rat_raster)[[1]]
  rat$priorities_cat <- cut(rat$ID, breaks = params$breaks)
  rat$priorities_cat <- factor(rat$priorities_cat,
                               levels = rev(levels(rat$priorities_cat)))
  levels(rat_raster) <- rat
  return(rat_raster)
}

create_legend <- function(raster, title, legend.reversed = FALSE) {

  params <- get_map_params(legend.reversed)

  if (legend.reversed) {
    raster <- create_raster_levels(raster)
  }

  raster_legend <- tm_shape(raster, is.master = TRUE) +
    tm_raster(title = title, palette = params$colors,
              labels = params$labels, breaks = params$breaks,
              auto.palette.mapping = FALSE,
              legend.show = TRUE) +
    tm_format_Europe(legend.only = TRUE, legend.position = c("left", "center"))
  return(raster_legend)
}

create_map <- function(raster, title) {

  params <- get_map_params()

  raster_map <- tm_shape(Europe, bbox = params$bbox, is.master = TRUE) +
    tm_fill("lightgrey") +
    tm_shape(raster) +
    tm_raster(palette = params$colors, labels = params$labels,
              breaks = params$breaks, auto.palette.mapping = FALSE,
              legend.show = FALSE) +
    tm_shape(Europe, bbox = params$bbox) +
    tm_borders(col = "black", lwd = 0.3) +
    tm_format_Europe(title = title, title.size = 4.0)
  return(raster_map)
}

create_map_list <- function(input_rasters) {
  map_list <- list()

  for (i in 1:length(input_rasters)) {
    raster_name <- names(input_rasters[i])[1]
    raster_path <- input_rasters[[raster_name]]
    map_list[[raster_name]] <- create_map(raster::raster(raster_path),
                                          title = LETTERS[i])
  }
  return(map_list)
}

print_map_list <- function(x, filepath, width = 1800, height = 1800,
                           nrow = 3, ncol = 3) {

  if ((nrow * ncol) != length(x)) {
    stop("The multiple of nrow and ncol must be the same as ",
         "length of x")
  }

  # Generate row/col indexes for the grid layout
  indexes <- mapply(c, rep(1:nrow, each = 3, times = 1),
                    rep(1:ncol, times = 3), SIMPLIFY = FALSE)

  message("Creating composite image ", filepath)
  # Create a new image
  png(filepath, width = width, height = height)
  # Establish a grid layout to which the image panels are printed
  # to
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow, ncol)))

  for (i in 1:length(x)) {
    message("Printing ", names(x[i])[1], " in grid position (",
            indexes[[i]][1], ", ", indexes[[i]][2], ")...")
    print(x[[i]], vp = viewport(layout.pos.row = indexes[[i]][1],
                                layout.pos.col = indexes[[i]][2]))
  }
  dev.off()
  message("All done")
  return(invisible(TRUE))
}

# Read in pixel-based rank data -------------------------------------------

# NOTE: the order matters -> the order is retained in the composite image
# layout.

input_rasters <- list(
  "rwr_raster_all" = "analyses/RWR/rwr_all_weights.tif",
  "zon_raster_all" = "analyses/zonation/priocomp/04_abf_all_wgt/04_abf_all_wgt_out/04_abf_all_wgt.rank.compressed.tif",
  "ilp_raster_all" = "analyses/ILP/ilp_all_weights.tif",
  "rwr_raster_es" = "analyses/RWR/rwr_es.tif",
  "zon_raster_es" = "analyses/zonation/priocomp/08_abf_es/08_abf_es_out/08_abf_es.rank.compressed.tif",
  "ilp_raster_es" = "analyses/ILP/ilp_es.tif",
  "rwr_raster_bd" = "analyses/RWR/rwr_bd.tif",
  "zon_raster_bd" = "analyses/zonation/priocomp/12_abf_bd/12_abf_bd_out/12_abf_bd.rank.compressed.tif",
  "ilp_raster_bd" = "analyses/ILP/ilp_bd.tif")

input_rasters_costs <- list(
  "rwr_raster_all_cst" = "analyses/RWR/rwr_all_weights_costs.tif",
  "zon_raster_all_cst" = "analyses/zonation/priocomp/06_abf_all_wgt_cst/06_abf_all_wgt_cst_out/06_abf_all_wgt_cst.rank.compressed.tif",
  "ilp_raster_all_cst" = "analyses/ILP/ilp_all_weights_costs.tif",
  "rwr_raster_es_cst" = "analyses/RWR/rwr_es_costs.tif",
  "zon_raster_es_cst" = "analyses/zonation/priocomp/10_abf_es_cst/10_abf_es_cst_out/10_abf_es_cst.rank.compressed.tif",
  "ilp_raster_es_cst" = "analyses/ILP/ilp_es_costs.tif",
  "rwr_raster_bd_cst" = "analyses/RWR/rwr_bd_costs.tif",
  "zon_raster_bd_cst" = "analyses/zonation/priocomp/14_abf_bd_cst/14_abf_bd_cst_out/14_abf_bd_cst.rank.compressed.tif",
  "ilp_raster_bd_cst" = "analyses/ILP/ilp_bd_costs.tif")

# Create maps -------------------------------------------------------------

maps <- create_map_list(input_rasters)
maps_costs <- create_map_list(input_rasters_costs)

# Save plots --------------------------------------------------------------

file_legend <- "reports/figures/figure03/01_figure_03_legend.png"
file_main <- "reports/figures/figure03/02_figure_03_main_nocosts.png"
file_main_costs <- "reports/figures/figure03/03_figure_03_main_costs.png"

print_map_list(maps, file_main, width = 1800, height = 1800,
               nrow = 3, ncol = 3)
print_map_list(maps_costs, file_main_costs, width = 1800, height = 1800,
               nrow = 3, ncol = 3)

# Create legend separately and only once based on a single raster map
rastermap_legend <- create_legend(raster::raster(input_rasters[[1]]),
                                  legend.reversed = TRUE,
                                  title = "Best X% of \nthe solution")

save_tmap(rastermap_legend, file_legend, width = 400, height = 600)

# Crop legend
img_legend <- magick::image_read(file_legend)
img_legend <- magick::image_crop(img_legend, geometry = "236x354+30+123")
magick::image_write(img_legend, path = file_legend)


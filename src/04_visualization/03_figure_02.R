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

# Read in pixel-based rank data -------------------------------------------

## RWR rank rasters
rwr_raster_all <- raster::raster("analyses/RWR/rwr_all_weights.tif")
rwr_raster_es <- raster::raster("analyses/RWR/rwr_es.tif")
rwr_raster_bd <- raster::raster("analyses/RWR/rwr_bd.tif")

## ZON rank rasters
zon_raster_all <- raster::raster("analyses/zonation/priocomp/04_abf_wgt/04_abf_wgt_out/04_abf_wgt.rank.compressed.tif")
zon_raster_es <- raster::raster("analyses/zonation/priocomp/06_abf_es/06_abf_es_out/06_abf_es.rank.compressed.tif")
zon_raster_bd <- raster::raster("analyses/zonation/priocomp/08_abf_bd/08_abf_bd_out/08_abf_bd.rank.compressed.tif")

## ILP rank rasters
ilp_raster_all <- raster::raster(x = "analyses/ILP/ilp_all_weights.tif")
ilp_raster_es <- raster::raster(x = "analyses/ILP/ilp_es.tif")
ilp_raster_bd <- raster::raster(x = "analyses/ILP/ilp_bd.tif")


# Set parameters ----------------------------------------------------------

inner_margins <- c(0.02, 0.02, 0.02, 0)
title_size <- 4.0

# Make a background map for all panels
tm_eur <- tm_shape(Europe) +
  tm_fill("lightgrey") +
  tm_format_Europe(inner.margins = inner_margins)

# Define a suitable bounding box
project_bbox <- matrix(c(2635899, 1386018, 6084606, 5307234),
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

# Make plots --------------------------------------------------------------

# ES

rwr_rastermap_es <- tm_eur +
  tm_shape(rwr_raster_es, bbox = project_bbox, is.master = TRUE) +
  tm_raster(palette = colors, labels = labels,
            breaks = breaks, auto.palette.mapping = FALSE,
            legend.show = FALSE) +
  tm_shape(Europe) +
  tm_borders(col = "black", lwd = 0.3) +
  tm_format_Europe(title = "A", title.size = title_size)

zon_rastermap_es <- tm_eur +
  tm_shape(zon_raster_es, bbox = project_bbox, is.master = TRUE) +
  tm_raster(palette = colors, labels = labels,
            breaks = breaks, auto.palette.mapping = FALSE,
            legend.show = FALSE) +
  tm_shape(Europe) +
  tm_borders(col = "black", lwd = 0.3) +
  tm_format_Europe(title = "B", title.size = title_size)

ilp_rastermap_es <- tm_eur +
  tm_shape(ilp_raster_es, bbox = project_bbox, is.master = TRUE) +
  tm_raster(palette = colors, labels = labels,
            breaks = breaks, auto.palette.mapping = FALSE,
            legend.show = FALSE) +
  tm_shape(Europe) +
  tm_borders(col = "black", lwd = 0.3) +
  tm_format_Europe(title = "C", title.size = title_size)

# ALL

rwr_rastermap_all <- tm_eur +
  tm_shape(rwr_raster_all, bbox = project_bbox, is.master = TRUE) +
    tm_raster(palette = colors, labels = labels,
              breaks = breaks, auto.palette.mapping = FALSE,
              legend.show = FALSE) +
  tm_shape(Europe) +
    tm_borders(col = "black", lwd = 0.3) +
  tm_format_Europe(title = "D", title.size = title_size)

rwr_rastermap_legend <- tm_eur +
  tm_shape(rwr_raster_all, bbox = project_bbox, is.master = TRUE) +
  tm_raster(title = "Top fraction", palette = colors,
            labels = labels, breaks = breaks,
            auto.palette.mapping = FALSE,
            legend.show = TRUE) +
  tm_format_Europe(legend.only = TRUE, legend.position = c("left", "center"))

zon_rastermap_all <- tm_eur +
  tm_shape(zon_raster_all, bbox = project_bbox, is.master = TRUE) +
  tm_raster(palette = colors, labels = labels,
            breaks = breaks, auto.palette.mapping = FALSE,
            legend.show = FALSE) +
  tm_shape(Europe) +
  tm_borders(col = "black", lwd = 0.3) +
  tm_format_Europe(title = "E", title.size = title_size)

ilp_rastermap_all <- tm_eur +
  tm_shape(ilp_raster_all, bbox = project_bbox, is.master = TRUE) +
  tm_raster(palette = colors, labels = labels,
            breaks = breaks, auto.palette.mapping = FALSE,
            legend.show = FALSE) +
  tm_shape(Europe) +
  tm_borders(col = "black", lwd = 0.3) +
  tm_format_Europe(title = "F", title.size = title_size)

# BD

rwr_rastermap_bd <- tm_eur +
  tm_shape(rwr_raster_bd, bbox = project_bbox, is.master = TRUE) +
  tm_raster(palette = colors, labels = labels,
            breaks = breaks, auto.palette.mapping = FALSE,
            legend.show = FALSE) +
  tm_shape(Europe) +
  tm_borders(col = "black", lwd = 0.3) +
  tm_format_Europe(title = "G", title.size = title_size)

zon_rastermap_bd <- tm_eur +
  tm_shape(zon_raster_bd, bbox = project_bbox, is.master = TRUE) +
  tm_raster(palette = colors, labels = labels,
            breaks = breaks, auto.palette.mapping = FALSE,
            legend.show = FALSE) +
  tm_shape(Europe) +
  tm_borders(col = "black", lwd = 0.3) +
  tm_format_Europe(title = "H", title.size = title_size)

ilp_rastermap_bd <- tm_eur +
  tm_shape(ilp_raster_bd, bbox = project_bbox, is.master = TRUE) +
  tm_raster(palette = colors, labels = labels,
            breaks = breaks, auto.palette.mapping = FALSE,
            legend.show = FALSE) +
  tm_shape(Europe) +
  tm_borders(col = "black", lwd = 0.3) +
  tm_format_Europe(title = "I", title.size = title_size)

# Save plots --------------------------------------------------------------

file_legend <- "reports/figures/02_figure_02_legend.png"
file_main <- "reports/figures/04_figure_02_main.png"

png(file_main, width = 1800, height = 1800)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3,3)))
print(rwr_rastermap_es, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(zon_rastermap_es, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(ilp_rastermap_es, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(rwr_rastermap_all, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(zon_rastermap_all, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(ilp_rastermap_all, vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
print(zon_rastermap_bd, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(rwr_rastermap_bd, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
print(ilp_rastermap_bd, vp = viewport(layout.pos.row = 3, layout.pos.col = 3))
dev.off()

save_tmap(rwr_rastermap_legend, file_legend, width = 400, height = 600)

# Combine images

# Combine images using magick (couldn't figure a better way...)
#img_main <- magick::image_read(file_main)
#img_legend <- magick::image_read(file_legend)
#img_legend <- magick::image_crop(img_legend, geometry = "400x600+0+100")
#img_composite <- magick::image_append(c(img_main, img_legend), stack = TRUE)
#magick::image_write(img_composite, path = file_composite)


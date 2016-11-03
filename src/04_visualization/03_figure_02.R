# Figure 2: Rank priority maps for different methods (RWR, ZON and ILP) in
# columns, and for different datasets in rows (ALL, ES, BD).

library(dplyr)
library(grid)
library(lazyeval)
library(maptools)
library(raster)
library(RColorBrewer)
library(rgdal)
library(sp)
library(tmap)
library(viridis)

data(Europe, land)

read_data <- function(x, value_field, group) {
  # Define CRS (ETRS89 / ETRS-LAEA), http://spatialreference.org/ref/epsg/3035/
  p4s <- "+init=epsg:3035"
  # Read in the data
  sp_obj <- rgdal::readOGR(dsn = x, layer = rgdal::ogrListLayers(x), p4s = p4s)
  if (names(sp_obj)[1] == "id") {
    names(sp_obj)[1] <- "ID"
  }
  # Filter out polygons that have NA for value_field
  sp_obj <- subset(sp_obj, !is.na(eval(parse(text = value_field))))
  # Rename value field
  new_value_field <- paste0(group, "_mean_rank")
  names(sp_obj)[names(sp_obj) == value_field] <- new_value_field
  # Select only certain columns
  sp_obj@data <- sp_obj@data[,c("ID", "NUTS_ID", "country", new_value_field)]
  return(sp_obj)
}

# Go through all the analyses results as aggregated to NUTS2 units

## RWR
rwr_all <- read_data(x = "analyses/RWR/rwr_eu26_all_weights_stats.geojson",
                     value_field = "X_mean", group = "rwr_all")
rwr_es <- read_data(x = "analyses/RWR/rwr_eu26_es_stats.geojson",
                    value_field = "X_mean", group = "rwr_es")
rwr_bd <- read_data(x = "analyses/RWR/rwr_eu26_bd_stats.geojson",
                    value_field = "X_mean", group = "rwr_bd")
rwr_combo <- sp::merge(rwr_all, rwr_es, by.x = "id", by.y = "id")
rwr_combo <- sp::merge(rwr_combo, rwr_bd, by.x = "id", by.y = "id")

## ZON
zon_all <- read_data(x = "analyses/zonation/priocomp/04_abf_wgt/04_abf_wgt_out/04_abf_wgt_nwout1.shp",
                     value_field = "Men_rnk", group = "zon_all")
zon_es <- read_data(x = "analyses/zonation/priocomp/06_abf_es/06_abf_es_out/06_abf_es_nwout1.shp",
                     value_field = "Men_rnk", group = "zon_es")
zon_bd <- read_data(x = "analyses/zonation/priocomp/08_abf_bd/08_abf_bd_out/08_abf_bd_nwout1.shp",
                     value_field = "Men_rnk", group = "zon_bd")
zon_combo <- sp::merge(zon_all, zon_es, by.x = "id", by.y = "id")
zon_combo <- sp::merge(zon_combo, zon_bd, by.x = "id", by.y = "id")

## ILP
ilp_all <- read_data(x = "analyses/ILP/ilp_eu26_all_weights_stats.geojson",
                     value_field = "X_mean", group = "ilp_all")
ilp_es <- read_data(x = "analyses/ILP/ilp_eu26_es_stats.geojson",
                    value_field = "X_mean", group = "ilp_es")
ilp_bd <- read_data(x = "analyses/ILP/ilp_eu26_bd_stats.geojson",
                    value_field = "X_mean", group = "ilp_bd")
ilp_combo <- sp::merge(ilp_all, ilp_es, by.x = "id", by.y = "id")
ilp_combo <- sp::merge(ilp_combo, ilp_bd, by.x = "id", by.y = "id")

inner_margins <- c(0.02, 0.02, 0.02, 0)
title_size <- 4.0
colors <- rev(RColorBrewer::brewer.pal(10, "RdYlBu"))
breaks <- seq(0, 1, 1 / length(colors))
labels <- format(round((100 - breaks * 100), 1), nsmall = 1)
labels <- cbind(labels[1:(length(labels) - 1)], labels[2:length(labels)])
labels <- apply(labels, 1, paste, collapse = " -")

# Make a background map for all panels
tm_eur <- tm_shape(Europe) +
  tm_fill("lightgrey") +
  tm_format_Europe(inner.margins = inner_margins)

rwr_maps <- tm_eur + tm_shape(rwr_combo, is.master = TRUE) +
  tm_polygons(c("rwr_es_mean_rank", "rwr_all_mean_rank", "rwr_bd_mean_rank"),
              style = "fixed", palette = colors, labels = labels,
              breaks = breaks, border.col = "lightgrey", lwd = 0.3,
              auto.palette.mapping = FALSE) +
  tm_layout(title.size = title_size) +
  tm_format_Europe(title = c("A", "B", "C"),
                   inner.margins = inner_margins,
                   legend.show = FALSE)

rwr_legend <- tm_shape(rwr_combo) +
  tm_polygons("rwr_all_mean_rank",
              style = "fixed", palette = colors, labels = labels,
              breaks = breaks, title = "Top fraction (%)",
              auto.palette.mapping = FALSE) +
  tm_format_Europe(legend.only = TRUE, legend.position = c("left", "center"))

zon_maps <- tm_eur + tm_shape(zon_combo, is.master = TRUE) +
  tm_polygons(c("zon_es_mean_rank", "zon_all_mean_rank", "zon_bd_mean_rank"),
              style = "fixed", palette = colors, labels = labels,
              breaks = breaks, border.col = "lightgrey", lwd = 0.3,
              auto.palette.mapping = FALSE) +
  tm_layout(title.size = title_size) +
  tm_format_Europe(title = c("D", "E", "F"),
                   inner.margins = inner_margins,
                   legend.show = FALSE)

ilp_maps <- tm_eur + tm_shape(ilp_combo, is.master = TRUE) +
  tm_polygons(c("ilp_es_mean_rank", "ilp_all_mean_rank", "ilp_bd_mean_rank"),
              style = "fixed", palette = colors, labels = labels,
              breaks = breaks, border.col = "lightgrey", lwd = 0.3,
              auto.palette.mapping = FALSE) +
  tm_layout(title.size = title_size) +
  tm_format_Europe(title = c("G", "H", "I"),
                   inner.margins = inner_margins,
                   legend.show = FALSE)

# FIXME: get the legend into the same image

png("reports/figures/02_figure_02.png", width = 1800, height = 1800)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3,1)))
print(rwr_maps, vp = viewport(layout.pos.row = 1))
print(zon_maps, vp = viewport(layout.pos.row = 2))
print(ilp_maps, vp = viewport(layout.pos.row = 3))
dev.off()

save_tmap(rwr_legend,"reports/figures/02_figure_02_legend.png", width = 400,
          height = 600)


# Variation ---------------------------------------------------------------

all_combo <- sp::merge(rwr_all, zon_all, by.x = "id", by.y = "id")
all_combo <- sp::merge(all_combo, ilp_all, by.x = "id", by.y = "id")

es_combo <- sp::merge(rwr_es, zon_es, by.x = "id", by.y = "id")
es_combo <- sp::merge(es_combo, ilp_es, by.x = "id", by.y = "id")

bd_combo <- sp::merge(rwr_bd, zon_bd, by.x = "id", by.y = "id")
bd_combo <- sp::merge(bd_combo, ilp_bd, by.x = "id", by.y = "id")

combo <- sp::merge(all_combo, es_combo, by.x = "id", by.y = "id")
combo <- sp::merge(combo, bd_combo, by.x = "id", by.y = "id")

combo$agg_mean <- apply(combo@data[, 4:12], 1, mean)
combo$agg_median <- apply(combo@data[, 4:12], 1, median)
combo$agg_sd <- apply(combo@data[, 4:12], 1, sd)
maptools::writePolyShape(combo, "analyses/comparison/variation.shp")


# Pixel-based rank maps ---------------------------------------------------

## RWR rank rasters
rwr_raster_all <- raster::raster("analyses/RWR/rwr_eu26_all_weights.tif")
rwr_raster_es <- raster::raster("analyses/RWR/rwr_eu26_es.tif")
rwr_raster_bd <- raster::raster("analyses/RWR/rwr_eu26_bd.tif")

## ZON rank rasters
zon_raster_all <- raster::raster("analyses/zonation/priocomp/04_abf_wgt/04_abf_wgt_out/04_abf_wgt.rank.compressed.tif")
zon_raster_es <- raster::raster("analyses/zonation/priocomp/06_abf_es/06_abf_es_out/06_abf_es.rank.compressed.tif")
zon_raster_bd <- raster::raster("analyses/zonation/priocomp/08_abf_bd/08_abf_bd_out/08_abf_bd.rank.compressed.tif")

## ILP rank rasters
ilp_raster_all <- raster::raster(x = "analyses/ILP/ilp_eu26_all_weights.tif")
ilp_raster_es <- raster::raster(x = "analyses/ILP/ilp_eu26_es.tif")
ilp_raster_bd <- raster::raster(x = "analyses/ILP/ilp_eu26_bd.tif")

breaks <- c(0, 0.2, 0.5, 0.75, 0.9, 0.95, 0.98, 1)
colors <- rev(RColorBrewer::brewer.pal(length(breaks) - 1, "RdYlBu"))
labels <- (100 - breaks * 100)
labels <- cbind(labels[1:(length(labels) - 1)], labels[2:length(labels)])
labels[,2] <- paste(labels[,2], "%")
labels[7,2] <- ""
labels <- apply(labels, 1, paste, collapse = " - ")
labels[7] <- gsub(" - ", " %", labels[7])

# Get proper bbox from rwr_all, otherwise too wide
project_bbox <- bbox(rwr_all)

rwr_rastermap_all <- tm_eur +
  tm_shape(rwr_raster_all, bbox = project_bbox, is.master = TRUE) +
    tm_raster(palette = colors, labels = labels,
              breaks = breaks, auto.palette.mapping = FALSE,
              legend.show = FALSE) +
  tm_shape(Europe) +
    tm_borders(col = "black", lwd = 0.3) +
  tm_format_Europe(title = "A", title.size = title_size)

rwr_rastermap_legend <- tm_eur +
  tm_shape(rwr_raster_all, bbox = project_bbox, is.master = TRUE) +
  tm_raster(title = "Top fraction", palette = colors,
            labels = labels, breaks = breaks,
            auto.palette.mapping = FALSE,
            legend.show = TRUE) +
  tm_format_Europe(legend.only = TRUE, legend.position = c("left", "center"))

rwr_rastermap_es <- tm_eur +
  tm_shape(rwr_raster_es, bbox = project_bbox, is.master = TRUE) +
    tm_raster(palette = colors, labels = labels,
              breaks = breaks, auto.palette.mapping = FALSE,
              legend.show = FALSE) +
  tm_shape(Europe) +
    tm_borders(col = "black", lwd = 0.3) +
  tm_format_Europe(title = "B", title.size = title_size)

rwr_rastermap_bd <- tm_eur +
  tm_shape(rwr_raster_bd, bbox = project_bbox, is.master = TRUE) +
  tm_raster(palette = z_legend$colors, labels = z_legend$labels,
            breaks = z_legend$values, auto.palette.mapping = FALSE,
            legend.show = FALSE) +
  tm_shape(Europe) +
  tm_borders(col = "black", lwd = 0.3) +
  tm_format_Europe(title = "C", title.size = title_size)

zon_rastermap_all <- tm_eur +
  tm_shape(zon_raster_all, bbox = project_bbox, is.master = TRUE) +
  tm_raster(palette = z_legend$colors, labels = z_legend$labels,
            breaks = z_legend$values, auto.palette.mapping = FALSE,
            legend.show = FALSE) +
  tm_shape(Europe) +
  tm_borders(col = "black", lwd = 0.3) +
  tm_format_Europe(title = "D", title.size = title_size)

zon_rastermap_es <- tm_eur +
  tm_shape(zon_raster_es, bbox = project_bbox, is.master = TRUE) +
  tm_raster(palette = z_legend$colors, labels = z_legend$labels,
            breaks = z_legend$values, auto.palette.mapping = FALSE,
            legend.show = FALSE) +
  tm_shape(Europe) +
  tm_borders(col = "black", lwd = 0.3) +
  tm_format_Europe(title = "E", title.size = title_size)

zon_rastermap_bd <- tm_eur +
  tm_shape(zon_raster_bd, bbox = project_bbox, is.master = TRUE) +
  tm_raster(palette = z_legend$colors, labels = z_legend$labels,
            breaks = z_legend$values, auto.palette.mapping = FALSE,
            legend.show = FALSE) +
  tm_shape(Europe) +
  tm_borders(col = "black", lwd = 0.3) +
  tm_format_Europe(title = "F", title.size = title_size)

ilp_rastermap_all <- tm_eur +
  tm_shape(ilp_raster_all, bbox = project_bbox, is.master = TRUE) +
  tm_raster(palette = z_legend$colors, labels = z_legend$labels,
            breaks = z_legend$values, auto.palette.mapping = FALSE,
            legend.show = FALSE) +
  tm_shape(Europe) +
  tm_borders(col = "black", lwd = 0.3) +
  tm_format_Europe(title = "G", title.size = title_size)

ilp_rastermap_es <- tm_eur +
  tm_shape(ilp_raster_es, bbox = project_bbox, is.master = TRUE) +
  tm_raster(palette = z_legend$colors, labels = z_legend$labels,
            breaks = z_legend$values, auto.palette.mapping = FALSE,
            legend.show = FALSE) +
  tm_shape(Europe) +
  tm_borders(col = "black", lwd = 0.3) +
  tm_format_Europe(title = "H", title.size = title_size)

ilp_rastermap_bd <- tm_eur +
  tm_shape(ilp_raster_bd, bbox = project_bbox, is.master = TRUE) +
  tm_raster(palette = z_legend$colors, labels = z_legend$labels,
            breaks = z_legend$values, auto.palette.mapping = FALSE,
            legend.show = FALSE) +
  tm_shape(Europe) +
  tm_borders(col = "black", lwd = 0.3) +
  tm_format_Europe(title = "I", title.size = title_size)

png("reports/figures/09_figure_02_rank_raster.png", width = 1800, height = 1800)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3,3)))
print(rwr_rastermap_all, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(rwr_rastermap_es, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(rwr_rastermap_bd, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(zon_rastermap_all, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(zon_rastermap_es, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(zon_rastermap_bd, vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
print(ilp_rastermap_all, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(ilp_rastermap_es, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
print(ilp_rastermap_bd, vp = viewport(layout.pos.row = 3, layout.pos.col = 3))
dev.off()

save_tmap(rwr_rastermap_legend,"reports/figures/09_figure_02_rank_raster_legend.png", width = 400,
          height = 600)

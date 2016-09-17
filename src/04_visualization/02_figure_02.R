# Figure 2: Rank priority maps for different methods (RWR, ZON and ILP) in
# columns, and for different datasets in rows (ALL, ES, BD).

library(dplyr)
library(grid)
library(lazyeval)
library(maptools)
library(sp)
library(tmap)
library(rgdal)
library(viridis)

read_data <- function(x, value_field, group) {
  # Define CRS (ETRS89 / ETRS-LAEA), http://spatialreference.org/ref/epsg/3035/
  p4s <- "+init=epsg:3035"
  # Read in the data
  sp_obj <- rgdal::readOGR(dsn = x, layer = rgdal::ogrListLayers(x), p4s = p4s)
  if (names(sp_obj)[1] == "id") {
    names(sp_obj)[1] <- "ID"
  }
  # Filter out polygons that have NA for value_field
  sp_obj <- subset(sp_obj, !is.na(value_field))
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

breaks <- 7

rwr_maps <- tm_shape(rwr_combo) +
  tm_polygons(c("rwr_es_mean_rank", "rwr_all_mean_rank", "rwr_bd_mean_rank"),
              palette = viridis::viridis(breaks),
              style = "fixed", breaks = seq(0, 1, 0.1),
              legend.is.portrait = TRUE) +
  tm_format_Europe(title = c("RWR ES", "RWR ALL", "RWR BD"),
                   inner.margins = c(0.02, 0.02, 0.02, 0),
                   legend.show = FALSE) +
  tm_style_grey()

zon_maps <- tm_shape(zon_combo) +
  tm_polygons(c("zon_es_mean_rank", "zon_all_mean_rank", "zon_bd_mean_rank"),
              palette = viridis::viridis(breaks),
              style = "fixed", breaks = seq(0, 1, 0.1),
              legend.is.portrait = TRUE) +
  tm_format_Europe(title = c("ZON ES", "ZON ALL", "ZON BD"),
                   inner.margins = c(0.02, 0.02, 0.02, 0),
                   legend.show = FALSE) +
  tm_style_grey()

ilp_maps <- tm_shape(ilp_combo) +
  tm_polygons(c("ilp_es_mean_rank", "ilp_all_mean_rank", "ilp_bd_mean_rank"),
              palette = viridis::viridis(breaks),
              style = "fixed", breaks = seq(0, 1, 0.1),
              legend.is.portrait = TRUE) +
  tm_format_Europe(title = c("ILP ES", "ILP ALL", "ILP BD"),
                   inner.margins = c(0.02, 0.02, 0.02, 0),
                   legend.show = FALSE) +
  tm_style_grey()

png("reports/figures/02_figure_02.png", width = 1800, height = 1800)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3,1)))
print(rwr_maps, vp = viewport(layout.pos.row = 1))
print(zon_maps, vp = viewport(layout.pos.row = 2))
print(ilp_maps, vp = viewport(layout.pos.row = 3))
dev.off()

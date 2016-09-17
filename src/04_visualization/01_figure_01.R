# Figure 1: Area of interest spanning 26 EU Member States.

library(maptools)
library(sp)
library(tmap)
library(viridis)

# Define CRS (ETRS89 / ETRS-LAEA), http://spatialreference.org/ref/epsg/3035/
crs_etrslaea <- CRS("+init=epsg:3035")
# Read in area-of-interest, a selceted subset of NUTS level 0 regions
shp_ds <- "data/processed/eurostat/nuts_level0/NUTS_RG_01M_2013_level0.shp"
aoi <- maptools::readShapePoly(shp_ds, proj4string = crs_etrslaea)

# Define color palette using viridis
col_pal <- viridis::viridis(3)

tm <- tm_shape(aoi) +
        tm_fill(col_pal[2]) +
        tm_borders("#4D4D4D", lwd = 0.6) +
        tm_format_Europe(title = "EU-26",
                         inner.margins = c(0, 0.05, 0.02, 0)) +
        tm_style_grey()

save_tmap(tm, "reports/figures/01_figure_01.png", width = 1200, height = 1400)

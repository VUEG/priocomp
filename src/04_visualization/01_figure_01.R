# Figure 1: Area of interest spanning 26 EU Member States.

library(maptools)
library(sp)
library(tmap)
library(viridis)

data(Europe)

# Define CRS (ETRS89 / ETRS-LAEA), http://spatialreference.org/ref/epsg/3035/
crs_etrslaea <- CRS("+init=epsg:3035")
# Read in area-of-interest, a selceted subset of NUTS level 0 regions
europe_ds <- "data/processed/eurostat/nuts_level0/NUTS_RG_01M_2013_level0.shp"
aoi_ds <- "data/processed/eurostat/nuts_level0/NUTS_RG_01M_2013_level0.shp"
aoi <- maptools::readShapePoly(shp_ds, proj4string = crs_etrslaea)

# Remove Liechtenstein
aoi <- subset(aoi, aoi$NUTS_ID != "LI")

# Define color palette using viridis
col_pal <- viridis::viridis(3)

tm_eur <- tm_shape(Europe) +
  tm_fill("lightgrey") +
  tm_format_Europe(inner.margins = c(0, 0.05, 0.02, 0))

tm_aoi <- tm_shape(aoi, is.master = TRUE) +
  tm_fill(col_pal[2]) +
  tm_text(text = "NUTS_ID", size = 0.6, col = "black", shadow = FALSE) +
  tm_borders("lightgrey", lwd = 0.3) +
  tm_format_Europe(inner.margins = c(0, 0.05, 0.02, 0))

tm <- tm_eur + tm_aoi + tm_scale_bar(position = c("left", "top"))

save_tmap(tm, "reports/figures/01_figure_01.png", width = 1200, height = 1400)

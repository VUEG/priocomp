# Figure S2: Cost raters based on European population density.
library(maptools)
library(raster)
library(sp)
library(tmap)
library(viridis)

data(Europe)

# Define CRS (ETRS89 / ETRS-LAEA), http://spatialreference.org/ref/epsg/3035/
crs_etrslaea <- CRS("+init=epsg:3035")
# Read in area-of-interest, a selceted subset of NUTS level 0 regions
europe_ds <- "data/processed/eurostat/nuts_level0/NUTS_RG_01M_2013_level0.shp"
aoi_ds <- "data/processed/eurostat/nuts_level0/NUTS_RG_01M_2013_level0.shp"
aoi <- maptools::readShapePoly(aoi_ds, proj4string = crs_etrslaea)

# Read in the cost data rasters
cost_raster <- raster::raster("data/processed/features/eea/pop_density/pop_density_v5.tif")

# Remove Liechtenstein
aoi <- subset(aoi, aoi$NUTS_ID != "LI")

tm_aoi <- tm_shape(Europe) +
  tm_fill("lightgrey") +
  tm_shape(cost_raster) +
  tm_raster(title = "Cost proxy", palette = viridis(8)) +
  tm_shape(aoi, is.master = TRUE) +
  tm_borders("lightgrey", lwd = 0.3) +
  tm_format_Europe(inner.margins = c(0, 0.05, 0.02, 0)) +
  tm_legend(position = c("left", "top"))

save_tmap(tm_aoi, "reports/figures/figureS02/01_figure_S02.png", width = 1200, height = 1200)

library(RColorBrewer)
library(htmlwidgets)
library(leaflet)
library(raster)
library(rgdal)
library(zonator)

source("src/utils.R")

priocomp_project <- load_zproject(root = "analyses/zonation/priocomp/",
                                  debug = TRUE)

output_dir <- "/home/jlehtoma/dev/git-projects/priocomp/analyses/zonation/priocomp/01_caz/01_caz_out"
dsn_geojson <- file.path(output_dir, "01_caz_nwout1.geojson")
dsn_rank_tif <- file.path(output_dir, "01_caz.rank.compressed.tif")
variant1_LSM <- readOGR(dsn_geojson, layer = ogrListLayers(dsn_geojson),
                        p4s = "+init=epsg:3035")
# Re-prject to WGS84
variant1_LSM_wgs84 <- spTransform(variant1_LSM, CRS("+init=epsg:4326"))

variant1_rank <- raster(dsn_rank_tif)

zpal <- colorBin(
  palette = zlegend("spectral")$colors,
  domain = c(0, 1)
)

prefecture_popup <- paste0("<strong>NUTS ID: </strong>",
                           variant1_LSM_wgs84$NUTS_ID,
                           "<br><strong>Mean rank: </strong>",
                           variant1_LSM_wgs84$Men_rnk,
                           "<br><strong>Spp dist. sum: </strong>",
                           variant1_LSM_wgs84$Spp_ds_,
                           "<br><strong>Spp dist. > 10%:    </strong>",
                           variant1_LSM_wgs84$Plus_10,
                           "<br><strong>Spp dist. > 1%: </strong>",
                           variant1_LSM_wgs84$Plus_1,
                           "<br><strong>Spp dist. > 0.1%: </strong>",
                           variant1_LSM_wgs84$Plus_01)

base_map <- leaflet(variant1_LSM_wgs84) %>% addTiles()

poly_map <- base_map %>%
  addPolygons(
    stroke = FALSE, fillOpacity = 0.8, smoothFactor = 0.5,
    color = ~zpal(Men_rnk), popup = prefecture_popup
  ) %>%
  addLegend("bottomright", pal = zpal, values = c(0, 1),
            title = "Mean rank",
            #labFormat = labelFormat(prefix = "$"),
            opacity = 1
  )

html_output <- gsub(".geojson", ".html", normalizePath(dsn_geojson))
saveWidget(poly_map, file = html_output)

raster_map <- base_map %>%
  addRasterImage(variant1_rank, colors = zpal, opacity = 0.8) %>%
  addLegend(pal = zpal, values = c(0, 1),
            title = "Rank priority")

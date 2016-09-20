library(maptools)
library(sp)
library(tmap)
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
  sp_obj <- subset(sp_obj, !is.na(eval(parse(text = value_field))))
  # Rename value field
  new_value_field <- paste0(group, "_mean_rank")
  names(sp_obj)[names(sp_obj) == value_field] <- new_value_field
  # Select only certain columns
  sp_obj@data <- sp_obj@data[,c("ID", "NUTS_ID", "country", new_value_field)]
  return(sp_obj)
}

# RWR
sp_rwr_all <- read_data("analyses/RWR/rwr_eu26_all_stats.geojson",
                        "X_mean", group = "rwr_all")
sp_rwr_all_wgt <- read_data("analyses/RWR/rwr_eu26_all_weights_stats.geojson",
                            "X_mean", group = "rwr_all_wgt")
sp_rwr_es <- read_data("analyses/RWR/rwr_eu26_es_stats.geojson",
                       "X_mean", group = "rwr_es")
sp_rwr_es <- read_data("analyses/RWR/rwr_eu26_bd_stats.geojson",
                       "X_mean", group = "rwr_bd")

# ZON

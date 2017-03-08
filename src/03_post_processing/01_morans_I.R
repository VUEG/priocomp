library(rgdal)
library(raster)
library(SpatialPack) # compute crh stats

source("../00_lib/utils.R")

# Setup -------------------------------------------------------------------

# Provide data
provide_file_paths <- list.files("../../data/processed/features/provide", pattern = "\\.tif$",
                                 full.names = TRUE, recursive = TRUE)
# Datadryad data
datadryad_file_paths <- list.files("../../data/processed/features/datadryad",
                                   pattern = "\\.tif$", full.names = TRUE,
                                   recursive = TRUE)

# UDR
udr_file_paths <- list.files("../../data/processed/features/udr", pattern = "\\.tif$",
                             full.names = TRUE, recursive = TRUE)

# ES and BD files
es_file_paths <- c(provide_file_paths, datadryad_file_paths)
bd_file_paths <- udr_file_paths
all_file_paths <- c(es_file_paths, bd_file_paths)

# Read into raster --------------------------------------------------------

morans <- {}
n_files <- length(all_file_paths)
for (i in 1:n_files) {
  infile <- all_file_paths[i]
  message(" [", i, "/", n_files, "] processing ", infile)
  morans[[i]] <- data.frame(feature = c(basename(infile)),
                            morans_i = c(raster::Moran(raster::raster(infile))))
}

morans <- do.call("rbind", morans)

filename <- file.path("..", "..", "data",
                      paste0("morans_I_values_",
                             nrow(morans), "_features_",
                             format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
                             ".csv"))

write.csv(morans, file = filename, row.names = FALSE)

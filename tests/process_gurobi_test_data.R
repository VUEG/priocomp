#!/usr/bin/env r

library(raster)

normalize <- function(x) {
  min <- raster::minValue(x)
  max <- raster::maxValue(x)
  return((x - min) / (max - min))
}

process_gurobi_results <- function(results, template_raster) {
  result_raster <- template_raster
  result_raster[] <- results$x
  # Fill in the nodata
  result_raster[is.na(template_raster)] <- NA
  # Make this a categorical raster
  result_raster <- raster::ratify(result_raster)
  rat <- levels(result_raster)[[1]]
  # In some cases (most notably when budget is 100%) all cells might be
  # selected. Check for this.
  if (nrow(rat) == 2) {
    rat$status <- c("Not Selected", "Selected")
  } else {
    rat$status <- c("Selected")
  }
  levels(result_raster) <- rat
  return(result_raster)
}

load("results_mc.RData")

data_dir <- "data/processed/features_ol_normalized"
if (!file.exists(data_dir)) {
  stop("Dir ", data_dir, " does not exist")
}

message("Reading in data...")
raster_files <- list.files(path = data_dir, pattern = ".+\\.tif$",
                           full.names = TRUE, recursive = TRUE)
message("Found ", length(raster_files), " rasters")

if (length(raster_files) == 0) {
  stop("No input files found in ", data_dir)
}

rasters <- raster::stack(raster_files)

template <- extent(rasters) %>%
  raster(nrows = nrow(rasters), ncols = ncol(rasters), vals = 1)
template[rasters[[1]] == 0] <- NA

mcp_ilp <- raster::stack(lapply(results_mc, process_gurobi_results, template))

# Sum up all the layers in the stack -> result is a selection frequency
mcp_ilp_hier <- sum(mcp_ilp, na.rm = TRUE)
# Replace 0s with NAs
mcp_ilp_hier[mcp_ilp_hier == 0] <- NA
# Normalize value into scale [0, 1]
mcp_ilp_hier <- normalize(mcp_ilp_hier)

raster::writeRaster(mcp_ilp_hier, "mcp_ilp_hier.tif")

#!/usr/bin/env r

library(magrittr)
library(protectr)
library(raster)
#library(rasterVis)
#library(viridis)

# Set raster options
rasterOptions(tmpdir = "/data/tmp/raster",
              progress = "text",
              chunksize = 10000000,
              maxmemory = 200000000)

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

message("Creating cost layer...")
cost <- extent(rasters) %>%
  raster(nrows = nrow(rasters), ncols = ncol(rasters), vals = 1)

# Fill areas normally NA (NoData) with cost of 0
cost[is.na(rasters[[1]])] <- 0

message("Converting NoData to zeros...")
# NAs (NoData) must be raplaced wiht 0s for GUROBI
#browser()
rasters_filled <- rasters
rasters_filled[is.na(rasters_filled)] <- 0
rasters_filled <- raster::stack(rasters_filled)

# Solve the maximum coverage problem for a range of target budgets
budgets <- c(0.05, 0.10)
results_mc <- list()
for (b in budgets) {
  message("Optimizing with target budget ", b)
  b_cells <- b * raster::cellStats(cost, "sum")
  results <- protectr::gurobi_maxcoverage(cost, rasters_filled, budget = b_cells)
  results_mc[[as.character(b)]] <- results
  message("Saving results...")
  save(results_mc, file = "results_mc.RData")
}

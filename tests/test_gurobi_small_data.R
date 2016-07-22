#!/usr/bin/env r

library(magrittr)
library(prioritizr)
library(raster)
#library(rasterVis)
#library(viridis)

# Set raster options
#rasterOptions(tmpdir = "/data/tmp/raster",
#              progress = "text",
#              chunksize = 100000000,
#              maxmemory = 200000000)

# Start load timer
load_start <- Sys.time()

data_dir <- "/home/jlehtoma/dev/git-projects/priocomp/data/processed/features/"
if (!file.exists(data_dir)) {
  stop("Dir ", data_dir, " does not exist")
}

message("Reading in data...")
raster_files <- list.files(path = data_dir, pattern = ".+\\.tif$",
                           full.names = TRUE, recursive = TRUE)[1:11]
message("Found ", length(raster_files), " rasters")
for (raster_file in raster_files) {
  message(raster_file)
}
if (length(raster_files) == 0) {
  stop("No input files found in ", data_dir)
}

rasters <- raster::stack(raster_files)

message("Creating cost layer...")
cost <- extent(rasters) %>%
  raster(nrows = nrow(rasters), ncols = ncol(rasters), vals = 1)

# Fill areas marked 0 with 0 in the cost data as well
cost[rasters[[1]] == 0] <- 0

#message("Converting NoData to zeros...")
# NAs (NoData) must be raplaced wiht 0s for GUROBI
#browser()
#rasters_filled <- rasters
#rasters_filled[is.na(rasters_filled)] <- 0
#rasters_filled <- raster::stack(rasters_filled)

# Stop load timer
load_time <- Sys.time() - load_start
# Start solver timer
solver_start <- Sys.time()

# Solve the maximum coverage problem for a range of target budgets
budgets <- seq(0.05, 1, 0.05)
results_mc <- list()
for (b in budgets) {
  message("Optimizing with target budget ", b)
  b_cells <- b * raster::cellStats(cost, "sum")
  mc_model <- prioritizr::maxcover_model(x = cost, features = rasters,
                                         budget = b_cells)
  mc_results <- prioritizr::prioritize(mc_model, gap = 0.001)
  results_mc[[as.character(b)]] <- mc_results
  message("Saving results...")
  save(results_mc, file = "results_mc.RData")
}
solver_time <- Sys.time() - solver_start

message("Results for ", length(raster_files), " rasters")
message("Load time: ")
print(load_time)
message("Solver time: ")
print(solver_time)


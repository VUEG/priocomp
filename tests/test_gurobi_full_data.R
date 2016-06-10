library(protectr)
library(raster)
library(rasterVis)
library(viridis)

data_dir <- "data/processed/features/ol_normalized/provide/"
raster_files <- list.files(path = data_dir, pattern = ".+\\.tif$",
                           full.names = TRUE, recursive = TRUE)

es_rasters <- raster::stack(raster_files)

cost <- extent(es_rasters) %>%
  raster(nrows = nrow(es_rasters), ncols = ncol(es_rasters), vals = 1)

# Alternatively, cost data can be simulated
# cost <- gaussian_field(r, 20, mean = 1000, variance = 500) #%>%
#  setNames("cost")

# Fill areas normalle NA (NoData) with cost of 0
cost[is.na(es_rasters[[1]])] <- 0

levelplot(cost, main = "Cost", margin = FALSE, scales = list(draw = FALSE),
          col.regions = viridis)

# NAs (NoData) must be raplaced wiht 0s for GUROBI
es_rasters_filled <- es_rasters
es_rasters_filled[is.na(es_rasters_filled)] <- 0
es_rasters_filled <- raster::stack(es_rasters_filled)

# Solve the maximum coverage problem for a range of target budgets
budgets <- c(0.05)
results_mc <- list()
for (b in budgets) {
  b_cells <- b * raster::cellStats(cost, "sum")
  results <- protectr::gurobi_maxcoverage(cost, es_rasters_filled, budget = b_cells)
  results_mc[[as.character(b)]] <- results
}

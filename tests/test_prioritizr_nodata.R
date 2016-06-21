library(prioritizr)
library(raster)

z_tutorial_dir <- "~/dev/git-data/zonation-tutorial"
z_tutorial_data_dir <- file.path(z_tutorial_dir, "data")
tutorial_files <- list.files(path = z_tutorial_data_dir,
                             pattern = "species.+\\.tif",
                             full.names = TRUE)
sp_rasters <- raster::stack(tutorial_files)
# For some mysterious reason, minmax isn't set automatically for sp_rasters[[1]]
# (species1)
sp_rasters[[1]] <- setMinMax(sp_rasters[[1]])

cost <- extent(sp_rasters) %>%
  raster(nrows = nrow(sp_rasters), ncols = ncol(sp_rasters), vals = 1)

cost[is.na(sp_rasters[[1]])] <- NA

b_cells <- 0.10 * raster::cellStats(cost, "sum")
mc_model <- prioritizr::maxcover_model(x = cost, features = sp_rasters,
                                       budget = b_cells)
mc_results <- prioritizr::prioritize(mc_model, gap = 0.001)
prioritizr::plot_selection(cost, mc_results$x, title = "Maximum Cover Solution")

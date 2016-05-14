library(ggplot2)
library(raster)
library(rasterVis)
options(scipen = 9999)

# Search for datasets. HACK: dealing with path depedning on whether we're running
# in the RStudio project root "priocomp" or Knitr folder "reports".

if (basename(getwd()) == "priocomp") {
  data_path <- "data/interim/rescaled"
} else {
  data_path <- "../../data/interim"
}

rescaled_files <- list.files(data_path, pattern = ".\\.tif$", recursive = TRUE,
                             full.names = TRUE)

# Create rasterstacks
rescaled_rasters <- raster::stack(rescaled_files)

notna_rescaled <- !is.na(rescaled_rasters)
aggregate_na <- raster::calc(notna_rescaled, sum)
raster::writeRaster(aggregate_na, filename = "data/interim/common/contextual_data_extent.tif")

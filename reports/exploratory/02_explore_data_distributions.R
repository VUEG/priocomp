library(ggplot2)
library(grid)
library(gridExtra)
library(raster)
library(rasterVis)
options(scipen = 9999)

# Search for datasets. HACK: dealing with path depedning on whether we're running
# in the RStudio project root "priocomp" or Knitr folder "reports".

if (basename(getwd()) == "priocomp") {
  data_path <- "data/interim/"
} else {
  data_path <- "../../data/interim/"
}

harmonized_files <- list.files(file.path(data_path, "harmonized"),
                               pattern = ".\\.tif$", recursive = TRUE,
                               full.names = TRUE)
rescaled_files <- list.files(file.path(data_path, "rescaled"),
                             pattern = ".\\.tif$", recursive = TRUE,
                             full.names = TRUE)

# Create rasterstacks
harmonized_rasters <- raster::stack(harmonized_files)
rescaled_rasters <- raster::stack(rescaled_files)

# Histograms

plots1 <- lapply(harmonized_files, function(x) {
  x <- raster(x)
  ggplot2::qplot(getValues(x), geom = "histogram", bins = 50, xlab = "Value",
                 main = names(x))
})

arr_plots1 <- do.call(grid.arrange, plots)
ggsave("reports/figures/harmonized_rasters.png", arr_plots1, width = 16, height = 16)

plots2 <- lapply(rescaled_files, function(x) {
  x <- raster(x)
  ggplot2::qplot(getValues(x), geom = "histogram", bins = 50, xlab = "Value",
                 main = names(x))
})

arr_plots2 <- do.call(grid.arrange, plots2)
ggsave("reports/figures/rescaled_rasters.png", arr_plots2, width = 16, height = 16)

# Boxplots

png(filename = "reports/figures/rescaled_rasters_boxplot.png", width = 1200,
    height = 1000)
par(mar = c(3, 18, 3, 3))
boxplot(rescaled_rasters, maxpixels = 1000000, horizontal = TRUE, las = 2,
        main = "Rescaled (normalized) values", ylab = "Value")
dev.off()

library(ggplot2)
library(grid)
library(gridExtra)
library(raster)
library(rasterVis)
options(scipen = 9999)

# Search for datasets. HACK: dealing with path depedning on whether we're running
# in the RStudio project root "priocomp" or Knitr folder "reports".

if (basename(getwd()) == "priocomp") {
  data_path <- "data/"
} else {
  data_path <- "../../data/"
}

harmonized_files <- list.files(file.path(data_path, "interim/harmonized"),
                               pattern = ".\\.tif$", recursive = TRUE,
                               full.names = TRUE)
processed_files <- list.files(file.path(data_path, "processed"),
                             pattern = ".\\.tif$", recursive = TRUE,
                             full.names = TRUE)

# Remove data_coverage.tif
processed_files <- processed_files[2:length(processed_files)]

# Histograms

plots1 <- lapply(harmonized_files, function(x) {
  x <- raster(x)
  ggplot2::qplot(getValues(x), geom = "histogram", bins = 50, xlab = "Value",
                 main = names(x))
})

arr_plots1 <- do.call(grid.arrange, plots1)
ggsave("reports/figures/harmonized_rasters_histogram.png", arr_plots1, width = 16, height = 16)

plots2 <- lapply(processed_files, function(x) {
  x <- raster(x)
  ggplot2::qplot(getValues(x), geom = "histogram", bins = 50, xlab = "Value",
                 main = names(x))
})

arr_plots2 <- do.call(grid.arrange, plots2)
ggsave("reports/figures/processed_rasters_histogram.png", arr_plots2, width = 16, height = 16)

# Boxplots

png(filename = "reports/figures/processed_rasters_boxplot.png", width = 1200,
    height = 1000)
par(mar = c(3, 18, 3, 3))
boxplot(raster::stack(processed_files), maxpixels = 1000000, horizontal = TRUE, las = 2,
        main = "Processed (harmonized + normalized) values", ylab = "Value")
dev.off()

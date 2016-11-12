library(gridExtra)
library(magick)
library(maptools)
library(RColorBrewer)
library(sp)
library(tmap)

data(Europe)

# Load data ---------------------------------------------------------------

nuts2_var_ds <- "analyses/comparison/nuts2_rank_variation.shp"
nuts2_var <- maptools::readShapePoly(nuts2_var_ds,
                                   proj4string = CRS("+init=epsg:3035"))

# Plot data ---------------------------------------------------------------

img_width <- 3000
img_height <- 1800
inner_margins <- c(0.02, 0.02, 0.02, -0.05)

# Make a background map for all panels
tm_eur <- tm_shape(Europe) +
  tm_fill("lightgrey") +
  tm_format_Europe(inner.margins = inner_margins)

title_size <- 1.5

mean_colors <- rev(RColorBrewer::brewer.pal(10, "RdYlBu"))
mean_breaks <- seq(0, 1, 1 / length(mean_colors))
mean_labels <- format((100 - mean_breaks * 100), nsmall = 0)
mean_labels <- cbind(mean_labels[1:(length(mean_labels) - 1)],
                     gsub(" ", "", mean_labels[2:length(mean_labels)]))
mean_labels[,2] <- paste(mean_labels[,2], "%")
mean_labels <- apply(mean_labels, 1, paste, collapse = " - ")

tm_mean_top <- tm_eur + tm_shape(nuts2_var, is.master = TRUE) +
  tm_polygons("agg_mean", title = "Mean best X% of \nthe solutions", style = "fixed",
              palette = mean_colors, labels = mean_labels, breaks = mean_breaks,
              border.col = "lightgrey", lwd = 0.3,
              auto.palette.mapping = FALSE) +
  tm_layout(title.size = title_size) +
  tm_format_Europe(title = "A", inner.margins = inner_margins,
                   legend.position = c("left", "top"),
                   legend.bg.color = "white",
                   title.position = c("right", "top"))

sd_colors <- rev(RColorBrewer::brewer.pal(7, "RdYlBu"))
sd_breaks <- seq(0, 0.35, 0.05)
sd_labels <- format(sd_breaks, nsmall = 0)
sd_labels <- cbind(sd_labels[1:(length(sd_labels) - 1)],
                     gsub(" ", "", sd_labels[2:length(sd_labels)]))
sd_labels[,2] <- sd_labels[,2]
sd_labels <- apply(sd_labels, 1, paste, collapse = " - ")

tm_sd_top <- tm_eur + tm_shape(nuts2_var, is.master = TRUE) +
  tm_polygons("agg_std", title = "SD mean priority rank",
              palette = sd_colors, labels = sd_labels, breaks = sd_breaks,
              border.col = "lightgrey", lwd = 0.3,
              auto.palette.mapping = FALSE) +
  tm_layout(title.size = title_size) +
  tm_format_Europe(title = "B", inner.margins = inner_margins,
                   legend.show = TRUE, legend.position = c("left", "top"),
                   legend.bg.color = "white",
                   title.position = c("right", "top"))

# Save map ----------------------------------------------------------------

file_mean_top <- "reports/figures/06_figure_03_A.png"
file_mean_sd <- "reports/figures/07_figure_03_B.png"
file_composite <- "reports/figures/08_figure_03.png"

save_tmap(tm_mean_top, file_mean_top, width = 1500, height = 1800)
save_tmap(tm_sd_top, file_mean_sd, width = 1500, height = 1800)

# Combine images using magick (couldn't figure a better way...)
img_mean_top <- magick::image_read(file_mean_top)
img_mean_sd <- magick::image_read(file_mean_sd)
img_composite <- magick::image_append(c(img_mean_top, img_mean_sd))
magick::image_write(img_composite, path = file_composite)

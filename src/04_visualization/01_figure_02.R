# Figure 2: Variation between rank patterns in NUTS2 units.
#
library(gridExtra)
library(magick)
library(maptools)
library(RColorBrewer)
library(sp)
library(tmap)

data(Europe)

# Helper functions --------------------------------------------------------

create_mean_map <- function(x, width, height, inner.margins = NULL) {
  mean_colors <- RColorBrewer::brewer.pal(10, "RdYlBu")
  mean_breaks <- seq(1, 0, -(1 / length(mean_colors)))
  mean_labels <- format((100 - mean_breaks * 100), nsmall = 0)
  mean_labels <- cbind(mean_labels[1:(length(mean_labels) - 1)],
                       gsub(" ", "", mean_labels[2:length(mean_labels)]))
  mean_labels[,2] <- paste(mean_labels[,2], "%")
  mean_labels <- apply(mean_labels, 1, paste, collapse = " - ")
  mean_labels[1] <- gsub(" 0 - ", "", mean_labels[1])

  # Manually categorize data to get an inverted legenf
  x$agg_mean_cat <-  cut(x$agg_mean, breaks = mean_breaks)
  x$agg_mean_cat <- factor(x$agg_mean_cat,
                           levels = rev(levels(x$agg_mean_cat)))

  tm_mean_top <- tm_shape(Europe) +
    tm_fill("lightgrey") +
    tm_format_Europe(inner.margins = inner.margins) +
    tm_shape(x, is.master = TRUE) +
    tm_polygons("agg_mean_cat", title = "Mean best X% of \nthe solutions", style = "fixed",
                palette = mean_colors, labels = mean_labels, breaks = mean_breaks,
                border.col = "lightgrey", lwd = 0.3,
                auto.palette.mapping = FALSE) +
    tm_layout(title.size = 1.5) +
    tm_format_Europe(title = "A", inner.margins = inner.margins,
                     legend.position = c("left", "top"),
                     legend.bg.color = "white",
                     title.position = c("right", "top"))
  return(tm_mean_top)
}

create_sd_map <-  function(x, upper.limit, step = 0.05, width, height,
                           inner.margins = NULL) {

  sd_colors <- rev(RColorBrewer::brewer.pal(upper.limit / step, "RdYlBu"))
  sd_breaks <- seq(0, upper.limit, step)
  sd_labels <- format(sd_breaks, nsmall = 0)
  sd_labels <- cbind(sd_labels[1:(length(sd_labels) - 1)],
                     gsub(" ", "", sd_labels[2:length(sd_labels)]))
  sd_labels[,2] <- sd_labels[,2]
  sd_labels <- apply(sd_labels, 1, paste, collapse = " - ")

  tm_sd_top <- tm_shape(Europe) +
    tm_fill("lightgrey") +
    tm_format_Europe(inner.margins = inner.margins) +
    tm_shape(x, is.master = TRUE) +
    tm_polygons("agg_std", title = "SD mean priority rank",
                palette = sd_colors, labels = sd_labels, breaks = sd_breaks,
                border.col = "lightgrey", lwd = 0.3,
                auto.palette.mapping = FALSE) +
    tm_layout(title.size = 1.5) +
    tm_format_Europe(title = "B", inner.margins = inner.margins,
                     legend.show = TRUE, legend.position = c("left", "top"),
                     legend.bg.color = "white",
                     title.position = c("right", "top"))
  return(tm_sd_top)
}

# Load data ---------------------------------------------------------------

nuts2_var_ds_nocosts <- "analyses/comparison/nuts2_rank_variation.shp"
nuts2_var_ds_costs <- "analyses/comparison/nuts2_rank_variation_costs.shp"

nuts2_var_nocosts <- maptools::readShapePoly(nuts2_var_ds_nocosts,
                                             proj4string = CRS("+init=epsg:3035"))
nuts2_var_costs <- maptools::readShapePoly(nuts2_var_ds_costs,
                                           proj4string = CRS("+init=epsg:3035"))

# Common parameters -------------------------------------------------------

img_width <- 3000
img_height <- 1800
inner_margins <- c(0.02, 0.02, 0.02, -0.05)

# Plot mean data ----------------------------------------------------------

tm_mean_top_nocosts <- create_mean_map(nuts2_var_nocosts, width = img_width,
                                       height = img_height,
                                       inner.margins = inner_margins)

tm_mean_top_costs <- create_mean_map(nuts2_var_costs, width = img_width,
                                     height = img_height,
                                     inner.margins = inner_margins)

# Plot SD data ------------------------------------------------------------

tm_sd_top_nocosts <- create_sd_map(nuts2_var_nocosts, upper.limit = 0.35,
                                   width = img_width, height = img_height,
                                   inner.margins = inner_margins)

tm_sd_top_costs <- create_sd_map(nuts2_var_costs, upper.limit = 0.20,
                                 step = 0.04,
                                 width = img_width,
                                 height = img_height,
                                 inner.margins = inner_margins)

# Save maps ---------------------------------------------------------------

file_mean_top_nocosts <- "reports/figures/figure02/01_figure_02_A_nocosts.png"
file_mean_sd_nocosts <- "reports/figures/figure02/02_figure_02_B_nocosts.png"
file_composite_nocosts <- "reports/figures/figure04/03_figure_02_nocosts.png"

file_mean_top_costs <- "reports/figures/figure02/04_figure_02_A_costs.png"
file_mean_sd_costs <- "reports/figures/figure02/05_figure_02_B_costs.png"
file_composite_costs <- "reports/figures/figure02/06_figure_02_costs.png"

save_tmap(tm_mean_top_nocosts, file_mean_top_nocosts, width = 1500, height = 1800)
save_tmap(tm_sd_top_nocosts, file_mean_sd_nocosts, width = 1500, height = 1800)
save_tmap(tm_mean_top_costs, file_mean_top_costs, width = 1500, height = 1800)
save_tmap(tm_sd_top_costs, file_mean_sd_costs, width = 1500, height = 1800)

# Combine images using magick (couldn't figure a better way...)
img_mean_top_nocosts <- magick::image_read(file_mean_top_nocosts)
img_mean_sd_nocosts <- magick::image_read(file_mean_sd_nocosts)
img_composite_nocosts <- magick::image_append(c(img_mean_top_nocosts,
                                                img_mean_sd_nocosts))
magick::image_write(img_composite_nocosts, path = file_composite_nocosts)

img_mean_top_costs <- magick::image_read(file_mean_top_costs)
img_mean_sd_costs <- magick::image_read(file_mean_sd_costs)
img_composite_costs <- magick::image_append(c(img_mean_top_costs,
                                                img_mean_sd_costs))
magick::image_write(img_composite_costs, path = file_composite_costs)

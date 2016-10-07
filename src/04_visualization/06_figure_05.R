library(maptools)
library(sp)
library(tmap)
library(viridis)

# Load data ---------------------------------------------------------------

nuts2_var_ds <- "analyses/comparison/nuts2_rank_variation.shp"
nuts2_var <- maptools::readShapePoly(nuts2_var_ds,
                                   proj4string = CRS("+init=epsg:3035"))

# Plot data ---------------------------------------------------------------

nuts2_var_maps <- tm_shape(nuts2_var) +
  tm_polygons(c("agg_mean", "agg_std"), title = c("", ""),
              palette = viridis::viridis(5),
              legend.is.portrait = TRUE) +
  tm_format_Europe(title = c("Mean priority rank", "Std or prriority ranks"),
                   inner.margins = c(0.02, 0.02, 0.02, -0.05),
                   legend.show = TRUE, legend.position = c("left", "top"),
                   legend.bg.color = "white") +
  tm_style_grey()

# Save map ----------------------------------------------------------------

save_tmap(nuts2_var_maps,"reports/figures/06_figure_05.png", width = 3000,
          height = 1800)

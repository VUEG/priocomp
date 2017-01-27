library(ggplot2)
library(dplyr)
library(raster)
library(rgdal)
library(zonator)

priocomp_project <- load_zproject(root = "analyses/zonation/priocomp/",
                                  debug = TRUE)

# We also need to original NUTS data to map which NUTS2 region is which
ppa_shp <- "data/processedeurostat/nuts_level2/NUTS_RG_01M_2013_level2.shp"
PPA_units_sp <- rgdal::readOGR(ppa_shp, ogrListLayers(ppa_shp))
nuts2_data <- as.data.frame(PPA_units_sp)

# Get the cost data
#cost_raster <- raster::raster("data/processed/features/eea/pop_density/pop_density_v5.tif")
#cost_raster_bin <- raster::raster("tests/tmp/pop_density_v5.tif")
#val_cr <- raster::getValues(cost_raster)
#unique(val_cr)
#hist(cost_raster)

# Helper functions --------------------------------------------------------

plot_features <- function(x) {
  plot(curves(results(x), groups = FALSE), main = x@name)
}

plot_groups <- function(x) {
  plot(curves(results(x), groups = TRUE), main = x@name, min = TRUE,
       mean = TRUE, max = TRUE)
}


# Globals -----------------------------------------------------------------

bdes_groups <- c("1" = "ES", "2" = "BD", "3" = "cost")
es_groups <- c("1" = "ES", "2" = "cost")
bd_groups <- c("1" = "amphibians", "2" = "birds",
               "3" = "mammals", "4" = "reptiles", "5" = "cost")

# Read in the PPA data
#variant1_ppa1 <- results1@ppa.lsm

# Combine with NUTS2-data
# variant1_ppa1 <- variant1_ppa1 %>%
#   left_join(., nuts2_data, by = c("Unit" = "ID"))

# p <- ggplot(variant1_ppa1, aes(x = Spp_distribution_sum, y = Mean_rank,
#                               size = Area, color = country))
# p + geom_point() + ylab("Mean rank\n") + xlab("\nDistribution sum")

# 04_abf_all_wgt --------------------------------------------------------------

variant4 <- zonator::get_variant(priocomp_project, 4)
zonator::groupnames(variant4) <- bdes_groups
plot_groups(variant4)

# 06_abf_all_wgt_cst ----------------------------------------------------------

variant6 <- zonator::get_variant(priocomp_project, 6)
zonator::groupnames(variant6) <- bdes_groups
plot_groups(variant6)

# 08_abf_es -------------------------------------------------------------------

variant8 <- zonator::get_variant(priocomp_project, 8)
zonator::groupnames(variant8) <- es_groups
plot_features(variant8)

# 10_abf_es_cst ---------------------------------------------------------------

variant10 <- zonator::get_variant(priocomp_project, 10)
zonator::groupnames(variant10) <- es_groups
plot_features(variant10)

# 12_abf_bd -------------------------------------------------------------------

variant12 <- zonator::get_variant(priocomp_project, 12)
zonator::groupnames(variant12) <- bd_groups
plot_groups(variant12)

# 14_abf_bd_cst ---------------------------------------------------------------

variant14 <- zonator::get_variant(priocomp_project, 14)
zonator::groupnames(variant14) <- bd_groups
plot_groups(variant14)

# 16_load_bd_all --------------------------------------------------------------

variant16 <- zonator::get_variant(priocomp_project, 16)
zonator::groupnames(variant16) <- bdes_groups
plot_groups(variant4)
plot_groups(variant16)

# 18_load_bd_all_cst ----------------------------------------------------------

variant18 <- zonator::get_variant(priocomp_project, 18)
zonator::groupnames(variant18) <- bdes_groups
plot_groups(variant6)
plot_groups(variant18)

# 20_load_bd_es ---------------------------------------------------------------

variant20 <- zonator::get_variant(priocomp_project, 20)
zonator::groupnames(variant20) <- es_groups
plot_features(variant8)
plot_features(variant20)

# 22_load_bd_es_cst -----------------------------------------------------------

variant22 <- zonator::get_variant(priocomp_project, 22)
zonator::groupnames(variant22) <- bdes_groups
plot_features(variant10)
plot_features(variant22)

# 23_load_rwr_all -------------------------------------------------------------

variant23 <- zonator::get_variant(priocomp_project, 23)
zonator::groupnames(variant23) <- bdes_groups
plot_groups(variant4)
plot_groups(variant23)

# 24_load_ilp_all -------------------------------------------------------------

variant24 <- zonator::get_variant(priocomp_project, 24)
zonator::groupnames(variant24) <- bdes_groups
plot_groups(variant4)
plot_groups(variant24)

# 26_load_ilp_all_cst ---------------------------------------------------------

variant26 <- zonator::get_variant(priocomp_project, 26)
zonator::groupnames(variant26) <- bdes_groups
plot_groups(variant6)
plot_groups(variant26)
3

library(ggplot2)
library(dplyr)
library(rgdal)
library(zonator)

priocomp_project <- load_zproject(root = "analyses/zonation/priocomp/",
                                  debug = TRUE)

# We also need to original NUTS data to map which NUTS2 region is which
ppa_shp <- "data/processed/eurostat/nuts/level2/NUTS_RG_01M_2013_level2.shp"
PPA_units_sp <- rgdal::readOGR(ppa_shp, ogrListLayers(ppa_shp))
nuts2_data <- as.data.frame(PPA_units_sp)


# Helper functions --------------------------------------------------------

plot_features <- function(x) {
  plot(curves(results(x), groups = FALSE), main = x@name)
}

plot_groups <- function(x) {
  plot(curves(results(x), groups = TRUE), main = x@name, min = TRUE,
       mean = TRUE, max = TRUE)
}


# Globals -----------------------------------------------------------------

bdes_groups <- c("1" = "ES", "2" = "BD")
bd_groups <- c("1" = "amphibians", "2" = "birds",
               "3" = "mammals", "4" = "reptiles")

# 01_caz -----------------------------------------------------------------------

variant1 <- zonator::get_variant(priocomp_project, 1)
groupnames(variant1) <- bdes_groups
plot_groups(variant1)

# Read in the PPA data
#variant1_ppa1 <- results1@ppa.lsm

# Combine with NUTS2-data
# variant1_ppa1 <- variant1_ppa1 %>%
#   left_join(., nuts2_data, by = c("Unit" = "ID"))

# p <- ggplot(variant1_ppa1, aes(x = Spp_distribution_sum, y = Mean_rank,
#                               size = Area, color = country))
# p + geom_point() + ylab("Mean rank\n") + xlab("\nDistribution sum")

# 02_abf_ --------------------------------------------------------------------

variant2 <- zonator::get_variant(priocomp_project, 2)
zonator::groupnames(variant2) <- bdes_groups
plot_groups(variant2)

# 03_caz_wgt ------------------------------------------------------------------

variant3 <- zonator::get_variant(priocomp_project, 3)
zonator::groupnames(variant3) <- bdes_groups
plot_groups(variant3)

# 04_abf_wgt ------------------------------------------------------------------

variant4 <- zonator::get_variant(priocomp_project, 4)
zonator::groupnames(variant4) <- bdes_groups
plot_groups(variant4)

# 05_caz_es -------------------------------------------------------------------

variant5 <- zonator::get_variant(priocomp_project, 5)
plot_features(variant5)

# 06_abf_es -------------------------------------------------------------------

variant6 <- zonator::get_variant(priocomp_project, 6)
plot_features(variant6)

# 07_caz_bd -------------------------------------------------------------------

variant7 <- zonator::get_variant(priocomp_project, 7)
zonator::groupnames(variant7) <- bd_groups
plot_groups(variant7)

# 08_caz_bd -------------------------------------------------------------------

variant8 <- zonator::get_variant(priocomp_project, 8)
zonator::groupnames(variant8) <- bd_groups
plot_groups(variant8)


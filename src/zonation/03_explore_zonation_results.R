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

# 01_caz -----------------------------------------------------------------------

variant1 <- zonator::get_variant(priocomp_project, 1)
groupnames(variant1) <- c("1" = "ES", "2" = "BD")
results1 <- results(variant1)
feature_curves1 <- curves(results1)
feature_grp_curves1 <- curves(results1, groups = TRUE)
plot(feature_grp_curves1)

# Read in the PPA data
variant1_ppa1 <- results@ppa.lsm

# Combine with NUTS2-data
variant1_ppa1 <- variant1_ppa1 %>%
  left_join(., nuts2_data, by = c("Unit" = "ID"))

p <- ggplot(variant1_ppa1, aes(x = Spp_distribution_sum, y = Mean_rank,
                               size = Area, color = country))
p + geom_point() + ylab("Mean rank\n") + xlab("\nDistribution sum")

# 02_abf_ --------------------------------------------------------------------

variant2 <- zonator::get_variant(priocomp_project, 2)
zonator::groupnames(variant2) <- c("1" = "ES", "2" = "BD")
results2 <- zonator::results(variant2)
feature_curves2 <- zonator::curves(results2)
feature_grp_curves2 <- zonator::curves(results2, groups = TRUE)
plot(feature_grp_curves2)

# 03_caz_wgt ------------------------------------------------------------------

variant3 <- zonator::get_variant(priocomp_project, 3)
zonator::groupnames(variant3) <- c("1" = "ES", "2" = "BD")
results3 <- zonator::results(variant3)
feature_curves3 <- zonator::curves(results3)
feature_grp_curves3 <- zonator::curves(results3, groups = TRUE)
plot(feature_grp_curves3)

# 04_abf_wgt ------------------------------------------------------------------

variant4 <- zonator::get_variant(priocomp_project, 4)
zonator::groupnames(variant4) <- c("1" = "ES", "2" = "BD")
results4 <- zonator::results(variant4)
feature_curves4 <- zonator::curves(results4)
feature_grp_curves4 <- zonator::curves(results4, groups = TRUE)
plot(feature_grp_curves4)

# 05_abf_es -------------------------------------------------------------------

variant5 <- zonator::get_variant(priocomp_project, 5)
results5 <- zonator::results(variant5)
feature_curves5 <- zonator::curves(results5)
plot(feature_curves5)

# 07_caz_bd -------------------------------------------------------------------

variant7 <- zonator::get_variant(priocomp_project, 7)
zonator::groupnames(variant7) <- c("1" = "amphibians", "2" = "birds",
                                   "3" = "mammals", "4" = "reptiles")
results7 <- zonator::results(variant7)
feature_curves7 <- zonator::curves(results7)
feature_grp_curves7 <- zonator::curves(results7, groups = TRUE)
plot(feature_grp_curves7)

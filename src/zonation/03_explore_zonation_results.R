library(ggplot2)
library(dplyr)
library(rgdal)
library(zonator)

priocomp_project <- load_zproject(root = "analyses/zonation/priocomp/",
                                  debug = TRUE)

# We also need to original NUTS data to map which NUTS2 region is which
ppa_shp <- "data/processed/nuts/NUTS_RG_01M_2013/level2/NUTS_RG_01M_2013_level2_subset.shp"
PPA_units_sp <- rgdal::readOGR(ppa_shp, ogrListLayers(ppa_shp))
nuts2_data <- as.data.frame(PPA_units_sp)

# 01_caz -----------------------------------------------------------------------

variant1 <- zonator::get_variant(priocomp_project, 1)
plot(zonator::curves(variant1))

# Read in the PPA data
variant1_results <- zonator::results(variant1)
variant1_ppa1 <- variant1_results@ppa.lsm

# Combine with NUTS2-data
variant1_ppa1 <- variant1_ppa1 %>%
  left_join(., nuts2_data, by = c("Unit" = "ID"))

p <- ggplot(variant1_ppa1, aes(x = Spp_distribution_sum, y = Mean_rank,
                               size = Area, color = country))
p + geom_point() + ylab("Mean rank\n") + xlab("\nDistribution sum")

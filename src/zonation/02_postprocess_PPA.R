library(rgdal)
library(zonator)

source("src/utils.R")

priocomp_root <- "analyses/zonation/priocomp/"
variants <- c("01_caz", "02_abf")
ppa_shp <- "data/processed/nuts/NUTS_RG_01M_2013/level2/NUTS_RG_01M_2013_level2_subset.shp"

postprocess_ppa(priocomp_root, variants, ppa_shp)

library(rgdal)
library(zonator)

source("src/utils.R")

priocomp_root <- "analyses/zonation/priocomp/"
variants <- c("01_caz", "02_abf", "03_caz_wgt", "04_abf_wgt")
ppa_shp <- "data/processed/eurostat/nuts/level2/NUTS_RG_01M_2013_level2.shp"

postprocess_ppa(priocomp_root, variants, ppa_shp)

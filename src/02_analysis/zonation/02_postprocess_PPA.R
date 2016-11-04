library(rgdal)
library(zonator)

source("src/00_lib/utils.R")

priocomp_root <- "analyses/zonation/priocomp/"
# variants <- c("01_caz", "02_abf", "03_caz_wgt", "04_abf_wgt", "05_caz_es",
#               "06_abf_es", "07_caz_bd", "08_abf_bd")

variants <- c("04_abf_wgt", "06_abf_es", "08_abf_bd")

ppa_shp <- "data/processed/eurostat/nuts_level2/NUTS_RG_01M_2013_level2.shp"

postprocess_ppa(priocomp_root, variants, ppa_shp)

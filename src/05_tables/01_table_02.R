library(dplyr)
library(DT)
library(tidyr)
library(zonator)

source("~/Dropbox/project.J-PriPA/japan-zsetup/R/00_lib/utils.R")

# Helper functions --------------------------------------------------------

get_perf_stats <- function(dat, level) {
  perf_dat <- dat %>%
    dplyr::group_by(variant) %>%
    dplyr::filter(pr_lost > level) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  return(perf_dat)
}

# Load variants and configure groups --------------------------------------

zproject <- zonator::load_zproject('analyses/zonation/priocomp/')

highlights <- c(0.1)

# RWR
rwr_nocosts <- get_stat_curves(zproject, variant_ids = 23, groups = FALSE) %>%
  dplyr::mutate(variant = gsub("23_load_rwr_all", "RWR_ALL", variant)) %>%
  dplyr::mutate(variant = factor(variant))

rwr_nocosts_grps <- get_stat_curves(zproject, variant_ids = 23, groups = TRUE) %>%
  dplyr::mutate(variant = gsub("23_load_rwr_all_1", "RWR_ALL_ES", variant)) %>%
  dplyr::mutate(variant = gsub("23_load_rwr_all_2", "RWR_ALL_BD", variant)) %>%
  dplyr::mutate(variant = factor(variant, levels = c("RWR_ALL_ES",
                                                     "RWR_ALL_BD"),
                                 ordered = TRUE))

rwr_nocosts <- dplyr::bind_rows(rwr_nocosts, rwr_nocosts_grps)

#rwr_costs <- get_stat_curves(zproject, variant_ids = 25, groups = TRUE) %>%
#  dplyr::mutate(variant = gsub("_1", "_es", variant)) %>%
#  dplyr::mutate(variant = gsub("_2", "_bd", variant)) %>%
#  dplyr::mutate(variant = factor(variant, levels = c("25_load_rwr_all_cst_es",
#                                                     "25_load_rwr_all_cst_bd"),
#                                 ordered = TRUE))

# Zonation
zon_nocosts <- get_stat_curves(zproject, variant_ids = 4, groups = FALSE) %>%
dplyr::mutate(variant = gsub("04_abf_all_wgt", "ZON_ALL", variant)) %>%
  dplyr::mutate(variant = factor(variant))

zon_nocosts_grps <- get_stat_curves(zproject, variant_ids = 4, groups = TRUE) %>%
  dplyr::mutate(variant = gsub("04_abf_all_wgt_1", "ZON_ALL_ES", variant)) %>%
  dplyr::mutate(variant = gsub("04_abf_all_wgt_2", "ZON_ALL_BD", variant)) %>%
  dplyr::mutate(variant = factor(variant, levels = c("ZON_ALL_ES",
                                                     "ZON_ALL_BD"),
                                 ordered = TRUE))

zon_nocosts <- dplyr::bind_rows(zon_nocosts, zon_nocosts_grps)

#zon_costs <- get_stat_curves(zproject, variant_ids = 6, groups = TRUE) %>%
#  dplyr::mutate(variant = gsub("_1", "_es", variant)) %>%
#  dplyr::mutate(variant = gsub("_2", "_bd", variant)) %>%
#  dplyr::mutate(variant = factor(variant, levels = c("06_abf_all_wgt_cst_es",
#                                                     "06_abf_all_wgt_cst_bd"),
#                                 ordered = TRUE))

# ILP
ilp_nocosts <- get_stat_curves(zproject, variant_ids = 24, groups = FALSE) %>%
  dplyr::mutate(variant = gsub("24_load_ilp_all", "ILP_ALL", variant)) %>%
  dplyr::mutate(variant = factor(variant))

ilp_nocosts_grps <- get_stat_curves(zproject, variant_ids = 24, groups = TRUE) %>%
  dplyr::mutate(variant = gsub("24_load_ilp_all_1", "ILP_ALL_ES", variant)) %>%
  dplyr::mutate(variant = gsub("24_load_ilp_all_2", "ILP_ALL_BD", variant)) %>%
  dplyr::mutate(variant = factor(variant, levels = c("ILP_ALL_ES",
                                                     "ILP_ALL_BD"),
                                 ordered = TRUE))

ilp_nocosts <- dplyr::bind_rows(ilp_nocosts, ilp_nocosts_grps)

#ilp_costs <- get_stat_curves(zproject, variant_ids = 26, groups = TRUE) %>%
#  dplyr::mutate(variant = gsub("_1", "_es", variant)) %>%
#  dplyr::mutate(variant = gsub("_2", "_bd", variant)) %>%
#  dplyr::mutate(variant = factor(variant, levels = c("26_load_ilp_all_cst_es",
#                                                     "26_load_ilp_all_cst_bd"),
#                                 ordered = TRUE))


# Combine data ------------------------------------------------------------

all_nocosts <- dplyr::bind_rows(rwr_nocosts, zon_nocosts, ilp_nocosts) %>%
  dplyr::mutate(variant = factor(variant, levels = c("RWR_ALL", "RWR_ALL_ES", "RWR_ALL_BD",
                                                     "ZON_ALL", "ZON_ALL_ES", "ZON_ALL_BD",
                                                     "ILP_ALL", "ILP_ALL_ES", "ILP_ALL_BD"),
                                 ordered = TRUE))

#all_costs <- dplyr::bind_rows(rwr_costs, zon_costs, ilp_costs) %>%
#  dplyr::mutate(variant = factor(variant, levels = c("25_load_rwr_all_cst_es",
#                                                     "25_load_rwr_all_cst_bd",
#                                                     "06_abf_all_wgt_cst_es",
#                                                     "06_abf_all_wgt_cst_bd",
#                                                     "26_load_ilp_all_cst_es",
#                                                     "26_load_ilp_all_cst_bd"),
#                                 ordered = TRUE))

all_perf_data <- get_perf_stats(all_nocosts, 0.9)

t2 <- all_perf_data %>%
  dplyr::select(variant, min, median, mean, max) %>%
  tidyr::gather(stat, value, -variant) %>%
  tidyr::spread(variant, value) %>%
  DT::datatable() %>%
  DT::formatPercentage(columns = 1:10, digits = 1)

DT::saveWidget(t2, '/home/jlehtoma/Dropbox/Projects/VU/OPERAs/SP2/priocomp/reports/tables/01_table_02.html')

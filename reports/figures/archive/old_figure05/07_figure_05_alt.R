library(dplyr)
library(ggplot2)
library(ggthemes)
library(tidyr)
library(zonator)
library(viridis)

source("~/Dropbox/project.J-PriPA/japan-zsetup/R/00_lib/utils.R")

# Load variants and configure groups --------------------------------------

zproject <- zonator::load_zproject('analyses/zonation/priocomp/')

highlights <- c(0.1)

# Separate ES and BD from ALL
v04_nocosts <- get_stat_curves(zproject, variant_ids = 4, groups = TRUE) %>%
  dplyr::mutate(variant = gsub("_1", "_es", variant)) %>%
  dplyr::mutate(variant = gsub("_2", "_bd", variant)) %>%
  dplyr::mutate(variant = factor(variant, levels = c("04_abf_all_wgt_es",
                                                     "04_abf_all_wgt_bd"),
                                 ordered = TRUE))

v06_costs <- get_stat_curves(zproject, variant_ids = 6, groups = TRUE) %>%
  dplyr::mutate(variant = gsub("_1", "_es", variant)) %>%
  dplyr::mutate(variant = gsub("_2", "_bd", variant)) %>%
  dplyr::mutate(variant = factor(variant, levels = c("06_abf_all_wgt_cst_es",
                                                     "06_abf_all_wgt_cst_bd"),
                                 ordered = TRUE))

# No costs
all_nocosts <- get_stat_curves(zproject, variant_ids = c(8, 12, 19, 20)) %>%
  dplyr::bind_rows(v04_nocosts, .) %>%
  dplyr::mutate(variant = factor(variant,
                                 levels = c("08_abf_es", "12_abf_bd",
                                            "04_abf_all_wgt_es", "04_abf_all_wgt_bd",
                                            "20_load_bd_es", "19_load_es_bd"),
                                 ordered = TRUE))
# Costs
all_costs <- get_stat_curves(zproject, variant_ids = c(10, 14, 21, 22)) %>%
  dplyr::bind_rows(v06_costs, .) %>%
  dplyr::mutate(variant = factor(variant,
                                 levels = c("10_abf_es_cst", "14_abf_bd_cst",
                                            "06_abf_all_wgt_cst_es", "06_abf_all_wgt_cst_bd",
                                            "22_load_bd_es_cst", "21_load_es_bd_cst"),
                                 ordered = TRUE))

p1 <- all_nocosts %>%
  plot_curves(title = "", non_param = TRUE, invert_x = TRUE,
              nrow = 3, ncol = 2, highlights = highlights,
              plot_min = TRUE, plot_max = TRUE,
              labels = c("ES", "BD",
                         "ES (rank ALL)", "BD (rank ALL)",
                         "ES (rank BD)", "BD (rank ES)"),
              ylab = "Fraction of feature occurrence level covered\n",
              xlab = "\nFraction of the landscape")

p2 <- all_costs %>%
  plot_curves(title = "", non_param = TRUE, invert_x = TRUE,
              nrow = 3, ncol = 2, highlights = highlights,
              plot_min = TRUE, plot_max = TRUE,
              labels = c("ES", "BD",
                         "ES (rank ALL)", "BD (rank ALL)",
                         "ES (rank BD)", "BD (rank ES)"),
              ylab = "Fraction of feature occurrence level covered\n",
              xlab = "\nFraction of the landscape")

# Save Figures ------------------------------------------------------------

ggsave("reports/figures/figure05/01_figure_05_nocosts_alt.png", p1,
       width = 6, height = 9)
ggsave("reports/figures/figure05/02_figure_05_costs_alt.png", p2,
       width = 6, height = 9)

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

# RWR
rwr_nocosts <- get_stat_curves(zproject, variant_ids = 23, groups = TRUE) %>%
  dplyr::mutate(variant = gsub("_1", "_es", variant)) %>%
  dplyr::mutate(variant = gsub("_2", "_bd", variant)) %>%
  dplyr::mutate(variant = factor(variant, levels = c("23_load_rwr_all_es",
                                                     "23_load_rwr_all_bd"),
                                 ordered = TRUE))

rwr_costs <- get_stat_curves(zproject, variant_ids = 25, groups = TRUE) %>%
  dplyr::mutate(variant = gsub("_1", "_es", variant)) %>%
  dplyr::mutate(variant = gsub("_2", "_bd", variant)) %>%
  dplyr::mutate(variant = factor(variant, levels = c("25_load_rwr_all_cst_es",
                                                     "25_load_rwr_all_cst_bd"),
                                 ordered = TRUE))

# Zonation
zon_nocosts <- get_stat_curves(zproject, variant_ids = 4, groups = TRUE) %>%
  dplyr::mutate(variant = gsub("_1", "_es", variant)) %>%
  dplyr::mutate(variant = gsub("_2", "_bd", variant)) %>%
  dplyr::mutate(variant = factor(variant, levels = c("04_abf_all_wgt_es",
                                                     "04_abf_all_wgt_bd"),
                                 ordered = TRUE))

zon_costs <- get_stat_curves(zproject, variant_ids = 6, groups = TRUE) %>%
  dplyr::mutate(variant = gsub("_1", "_es", variant)) %>%
  dplyr::mutate(variant = gsub("_2", "_bd", variant)) %>%
  dplyr::mutate(variant = factor(variant, levels = c("06_abf_all_wgt_cst_es",
                                                     "06_abf_all_wgt_cst_bd"),
                                 ordered = TRUE))

# ILP
ilp_nocosts <- get_stat_curves(zproject, variant_ids = 24, groups = TRUE) %>%
  dplyr::mutate(variant = gsub("_1", "_es", variant)) %>%
  dplyr::mutate(variant = gsub("_2", "_bd", variant)) %>%
  dplyr::mutate(variant = factor(variant, levels = c("24_load_ilp_all_es",
                                                     "24_load_ilp_all_bd"),
                                 ordered = TRUE))

ilp_costs <- get_stat_curves(zproject, variant_ids = 26, groups = TRUE) %>%
  dplyr::mutate(variant = gsub("_1", "_es", variant)) %>%
  dplyr::mutate(variant = gsub("_2", "_bd", variant)) %>%
  dplyr::mutate(variant = factor(variant, levels = c("26_load_ilp_all_cst_es",
                                                     "26_load_ilp_all_cst_bd"),
                                 ordered = TRUE))


# Combine data ------------------------------------------------------------

all_nocosts <- dplyr::bind_rows(rwr_nocosts, zon_nocosts, ilp_nocosts) %>%
  dplyr::mutate(variant = factor(variant, levels = c("23_load_rwr_all_es",
                                                     "23_load_rwr_all_bd",
                                                     "04_abf_all_wgt_es",
                                                     "04_abf_all_wgt_bd",
                                                     "24_load_ilp_all_es",
                                                     "24_load_ilp_all_bd"),
                                 ordered = TRUE))

all_costs <- dplyr::bind_rows(rwr_costs, zon_costs, ilp_costs) %>%
  dplyr::mutate(variant = factor(variant, levels = c("25_load_rwr_all_cst_es",
                                                     "25_load_rwr_all_cst_bd",
                                                     "06_abf_all_wgt_cst_es",
                                                     "06_abf_all_wgt_cst_bd",
                                                     "26_load_ilp_all_cst_es",
                                                     "26_load_ilp_all_cst_bd"),
                                 ordered = TRUE))

# Plot --------------------------------------------------------------------

labels <- c("RWR_ALL ES", "RWR_ALL BD",
            "ZON_ALL ES", "ZON_ALL BD",
            "ILP_ALL ES", "ILP_ALL BD")

p1 <- all_nocosts %>%
  plot_curves(title = "No costs", non_param = TRUE, invert_x = TRUE,
              nrow = 3, ncol = 2, highlights = highlights,
              plot_min = TRUE, plot_max = TRUE,
              labels = labels,
              ylab = "Fraction of feature occurrence level covered\n",
              xlab = "\nFraction of the landscape")

p2 <- all_costs %>%
  plot_curves(title = "Costs", non_param = TRUE, invert_x = TRUE,
              nrow = 3, ncol = 2, highlights = highlights,
              plot_min = TRUE, plot_max = TRUE,
              labels = labels,
              ylab = "Fraction of feature occurrence level covered\n",
              xlab = "\nFraction of the landscape")

# Save Figures ------------------------------------------------------------

ggsave("reports/figures/figure05/02_figure_05_nocosts_methods.png", p1,
       width = 6, height = 9)
ggsave("reports/figures/figure05/03_figure_05_costs_methods.png", p2,
       width = 6, height = 9)


# Extra - get exact figures -----------------------------------------------

get_perf_stats <- function(dat, level) {
  perf_dat <- dat %>%
    dplyr::group_by(variant) %>%
    dplyr::filter(pr_lost > level) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  return(perf_dat)
}

make_long <- function(dat) {
  dat_long <- dat %>%
    tidyr::gather(stat, value, -pr_lost, -variant) %>%
    dplyr::mutate(variant = gsub("24_load_ilp_all", "ILP", variant)) %>%
    dplyr::mutate(variant = gsub("04_abf_all_wgt", "ZON", variant)) %>%
    dplyr::mutate(variant = gsub("23_load_rwr_all", "RWR", variant)) %>%
    dplyr::mutate(variant = toupper(variant)) %>%
    tidyr::separate(variant, into = c("method", "group")) %>%
    dplyr::mutate(method = factor(method, levels = c("ILP", "ZON", "RWR"),
                                  ordered = TRUE)) %>%
    dplyr::mutate(group = factor(group, levels = c("ES", "BD"),
                                 ordered = TRUE))
  return(dat_long)
}

plot_perf_stats <- function(dat, title) {
  p <- ggplot(dat, aes(x = method, y = value, fill = group)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ stat, scales = "free_x") + coord_flip() +
    ggtitle(title)
  return(p)
}

# Top 2%
perf_2 <- get_perf_stats(all_nocosts, 0.98)
plot_perf_stats(make_long(perf_2), "Top 2%")

# Top 10%
perf_10 <- get_perf_stats(all_nocosts, 0.90)
plot_perf_stats(make_long(perf_10), "Top 10%")

# Top 25%
perf_25 <- get_perf_stats(all_nocosts, 0.75)
plot_perf_stats(make_long(perf_25), "Top 25%")


# Statistical significance ------------------------------------------------

rwr_all_nocosts <- zonator::get_variant(zproject, 23)
rwr_all_nocosts_curves <- zonator::curves(rwr_all_nocosts)
rwr_es_nocosts <- rwr_all_nocosts_curves %>%
  dplyr::select(pr_lost, woodprod_average:species_richness_vascular_plants)
rwr_bd_nocosts <- rwr_all_nocosts_curves %>%
  dplyr::select(pr_lost, alytes_cisternasii:zamenis_situla)

rwr_es_nocosts_top10 <- rwr_es_nocosts %>%
  dplyr::filter(pr_lost > 0.90) %>%
  dplyr::slice(1) %>%
  dplyr::select(-pr_lost) %>%
  tidyr::gather(feature, value) %>%
  dplyr::mutate(group = "ES")

rwr_bd_nocosts_top10 <- rwr_bd_nocosts %>%
  dplyr::filter(pr_lost > 0.90) %>%
  dplyr::slice(1) %>%
  dplyr::select(-pr_lost) %>%
  tidyr::gather(feature, value) %>%
  dplyr::mutate(group = "BD")

rwr_all_nocosts_data <- dplyr::bind_rows(rwr_es_nocosts_top10, rwr_bd_nocosts_top10)

wilcox.test(value ~ group, data = rwr_all_nocosts_data)

library(dplyr)
library(ggplot2)
library(ggthemes)
library(tidyr)
library(zonator)
library(viridis)


# Load variants and configure groups --------------------------------------

zproject <- zonator::load_zproject('analyses/zonation/priocomp/')

all_groups <- c("1" = "ZON_ES", "2" = "ZON_BD", "3" = "ZON_CST")
es_groups <- c("1" = "ZON_ES", "2" = "ZON_CST")
bd_taxon_groups <- c("1" = "amphibians", "2" = "birds", "3" = "mammals",
                     "4" = "reptiles", "5" = "ZON_CST")
bd_groups <- c("1" = "ZON_BD", "2" = "ZON_CST")

# ALL
v04_abf_all_wgt <- zonator::get_variant(zproject, 4)
v06_abf_all_wgt_cst <- zonator::get_variant(zproject, 6)
zonator::groupnames(v04_abf_all_wgt) <- all_groups
zonator::groupnames(v06_abf_all_wgt_cst) <- all_groups

# ES
v08_abf_es <- zonator::get_variant(zproject, 8)
v10_abf_es_cst <- zonator::get_variant(zproject, 10)
zonator::groupnames(v08_abf_es) <- es_groups
zonator::groupnames(v10_abf_es_cst) <- es_groups

# BD
v12_abf_bd <- zonator::get_variant(zproject, 12)
v14_abf_bd_cst <- zonator::get_variant(zproject, 14)
# NOTE: variants 12 and 14 have been run with taxon groups. Convert
# these to more generic groups.
generic_groups <- c(rep(1, zonator::nfeatures(v12_abf_bd) - 1), 2)
zonator::groups(v12_abf_bd) <- generic_groups
zonator::groupnames(v12_abf_bd) <- bd_groups
zonator::groups(v14_abf_bd_cst) <- generic_groups
zonator::groupnames(v14_abf_bd_cst) <- bd_groups

# Load ranking from ES, features from BD
v19_load_es_bd <- zonator::get_variant(zproject, 19)
zonator::groups(v19_load_es_bd) <- generic_groups
zonator::groupnames(v19_load_es_bd) <- bd_groups

# Load ranking from BD, features from ES
v20_load_bd_es <- zonator::get_variant(zproject, 20)
zonator::groupnames(v20_load_bd_es) <- es_groups

# Load ranking from ES, features from BD, COSTS
v21_load_es_bd_cst <- zonator::get_variant(zproject, 21)
zonator::groups(v21_load_es_bd_cst) <- generic_groups
zonator::groupnames(v21_load_es_bd_cst) <- bd_groups

# Load ranking from BD, features from ES, COSTS
v22_load_bd_es_cst <- zonator::get_variant(zproject, 22)
zonator::groupnames(v22_load_bd_es_cst) <- es_groups

# Get group curves data ---------------------------------------------------

v04_grp_crvs <- zonator::curves(v04_abf_all_wgt, groups = TRUE)
v06_grp_crvs <- zonator::curves(v06_abf_all_wgt_cst, groups = TRUE)
v08_grp_crvs <- zonator::curves(v08_abf_es, groups = TRUE)
v10_grp_crvs <- zonator::curves(v10_abf_es_cst, groups = TRUE)
v12_grp_crvs <- zonator::curves(v12_abf_bd, groups = TRUE)
v14_grp_crvs <- zonator::curves(v14_abf_bd_cst, groups = TRUE)
v19_grp_crvs <- zonator::curves(v19_load_es_bd, groups = TRUE)
v20_grp_crvs <- zonator::curves(v20_load_bd_es, groups = TRUE)
v21_grp_crvs <- zonator::curves(v21_load_es_bd_cst, groups = TRUE)
v22_grp_crvs <- zonator::curves(v22_load_bd_es_cst, groups = TRUE)

# Re-arrange data ---------------------------------------------------------

v04_grp_mean <- v04_grp_crvs %>%
  dplyr::select(pr_lost, ZON_ES = mean.ZON_ES, ZON_BD = mean.ZON_BD) %>%
  tidyr::gather(variant, ave_pr, -pr_lost)
v06_grp_mean <- v06_grp_crvs %>%
  dplyr::select(pr_lost, ZON_ES = mean.ZON_ES, ZON_BD = mean.ZON_BD) %>%
  tidyr::gather(variant, ave_pr, -pr_lost)

v08_mean <- v08_grp_crvs %>%
  dplyr::select(pr_lost, ave_pr = mean.ZON_ES) %>%
  dplyr::mutate(variant = "ES")
v10_mean <- v10_grp_crvs %>%
  dplyr::select(pr_lost, ave_pr = mean.ZON_ES) %>%
  dplyr::mutate(variant = "ES")

v12_mean <- v12_grp_crvs %>%
  dplyr::select(pr_lost, ave_pr = mean.ZON_BD) %>%
  dplyr::mutate(variant = "BD")
v14_mean <- v14_grp_crvs %>%
  dplyr::select(pr_lost, ave_pr = mean.ZON_BD) %>%
  dplyr::mutate(variant = "BD")

v19_mean <- v19_grp_crvs %>%
  dplyr::select(pr_lost, ave_pr = mean.ZON_BD) %>%
  dplyr::mutate(variant = "BD (rank ES)")
v20_mean <- v20_grp_crvs %>%
  dplyr::select(pr_lost, ave_pr = mean.ZON_ES) %>%
  dplyr::mutate(variant = "ES (rank BD)")
v21_mean <- v21_grp_crvs %>%
  dplyr::select(pr_lost, ave_pr = mean.ZON_BD) %>%
  dplyr::mutate(variant = "BD (rank ES)")
v22_mean <- v22_grp_crvs %>%
  dplyr::select(pr_lost, ave_pr = mean.ZON_ES) %>%
  dplyr::mutate(variant = "ES (rank BD)")

# Combine data ------------------------------------------------------------

# No costs
v04_grp_mean$variant <- gsub("ZON_BD", "BD (rank ALL)", v04_grp_mean$variant)
v04_grp_mean$variant <- gsub("ZON_ES", "ES (rank ALL)", v04_grp_mean$variant)

datag_perf_nocosts <- dplyr::bind_rows(list(v04_grp_mean, v08_mean, v12_mean,
                                       v19_mean, v20_mean))

# Costs
v06_grp_mean$variant <- gsub("ZON_BD", "BD (rank ALL)", v06_grp_mean$variant)
v06_grp_mean$variant <- gsub("ZON_ES", "ES (rank ALL)", v06_grp_mean$variant)

datag_perf_costs <- dplyr::bind_rows(list(v06_grp_mean, v10_mean, v14_mean,
                                          v21_mean, v22_mean))


# Plot mean curves --------------------------------------------------------

x_lab <- "\nFraction of the landscape"
x_scale <- scale_x_continuous(breaks = seq(0, 1, 0.2),
                              labels = paste(100 * seq(1, 0, -0.2), "%"))
y_scale <- scale_y_continuous(breaks = seq(0, 1, 0.2),
                              labels = paste(100 * seq(0, 1, 0.2), "%"))

# No costs
p1 <- ggplot2::ggplot(datag_perf_nocosts, aes(x = pr_lost, y = ave_pr,
                                          color = variant,
                                          linetype = variant)) +
  geom_vline(xintercept = 0.9, alpha = 0.5, linetype = 3) +
  annotate("text", x = 0.93, y = 1, label = "10%", size = 3) +
  geom_line(size = 0.9) + xlab(x_lab) +
  scale_linetype_manual("", values = rep(1:3, 2)) +
  scale_color_manual("", values =  rev(rep(viridis(2, end = 0.7), 3, each = 3))) +
  x_scale + y_scale + ylab("Average feature distribution covered\n") +
  theme_minimal() +
  theme(legend.position = c(0.25, 0.075),
        legend.justification = c(0.5, 0),
        legend.key.width = unit(1,"cm"))

# Costs
p2 <- ggplot2::ggplot(datag_perf_costs, aes(x = pr_lost, y = ave_pr,
                                      color = variant,
                                      linetype = variant)) +
  geom_vline(xintercept = 0.9, alpha = 0.5, linetype = 3) +
  annotate("text", x = 0.93, y = 1, label = "10%", size = 3) +
  geom_line(size = 0.9) + xlab(x_lab) +
  scale_linetype_manual("", values = rep(1:3, 2)) +
  scale_color_manual("", values =  rev(rep(viridis(2, end = 0.7), 3, each = 3))) +
  x_scale + y_scale + ylab("Average feature distribution covered\n") +
  theme_minimal() +
  theme(legend.position = c(0.25, 0.075),
        legend.justification = c(0.5, 0),
        legend.key.width = unit(1,"cm"))

# Save Figure -------------------------------------------------------------

ggsave("reports/figures/figure05/01_figure_05_nocosts.png", p1,
       width = 6, height = 6)
ggsave("reports/figures/figure05/02_figure_05_costs.png", p2,
       width = 6, height = 6)

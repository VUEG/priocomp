library(dplyr)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(tidyr)
library(zonator)
library(viridis)

# Load variants and get curves data ---------------------------------------

zproject <- zonator::load_zproject('analyses/zonation/priocomp/')

v04_abf_all_wgt <- zonator::get_variant(zproject, 4)
v06_abf_es <- zonator::get_variant(zproject, 6)
v08_abf_bd <- zonator::get_variant(zproject, 8)

groupnames(v04_abf_all_wgt) <- c("1" = "ES", "2" = "BD")
v04_crvs <- zonator::curves(v04_abf_all_wgt, groups = FALSE)
v04_grp_crvs <- zonator::curves(v04_abf_all_wgt, groups = TRUE)

v06_crvs <- zonator::curves(v06_abf_es, groups = FALSE)

v08_crvs <- zonator::curves(v08_abf_bd, groups = FALSE)

# Re-arrange data ---------------------------------------------------------

v04_mean <- v04_crvs %>%
  dplyr::select(pr_lost, ave_pr, w_pr) %>%
  dplyr::mutate(variant = "ZON")

# Construct a version where groups {"ES", "BD"} are separated
v04_grp_mean <- v04_grp_crvs %>%
  dplyr::select(pr_lost, mean.ES, mean.BD) %>%
  tidyr::gather(variant, ave_pr, -pr_lost) %>%
  dplyr::mutate(variant = gsub("mean\\.", "Full ", variant))

# Get the average proportion remaining over all features in respective
# variants
v06_mean <- v06_crvs %>%
  dplyr::select(pr_lost, ave_pr) %>%
  dplyr::mutate(variant = "ES only")

v08_mean <- v08_crvs %>%
  dplyr::select(pr_lost, ave_pr) %>%
  dplyr::mutate(variant = "BD only")

# Pre-load ----------------------------------------------------------------

# zonator can't deal with pre-loaded variants, so manually load up the
# curves data

v11_crvs <- zonator::read_curves("analyses/zonation/priocomp/11_load_abf_es_bd/11_load_abf_es_bd_out/11_load_abf_es_bd.curves.txt")
v11_mean <- v11_crvs %>%
  dplyr::select(pr_lost, ave_pr) %>%
  dplyr::mutate(variant = "BD only (rank ES only)")

v12_crvs <- zonator::read_curves("analyses/zonation/priocomp/12_load_abf_bd_es/12_load_abf_bd_es_out/12_load_abf_bd_es.curves.txt")
v12_mean <- v12_crvs %>%
  dplyr::select(pr_lost, ave_pr) %>%
  dplyr::mutate(variant = "ES only (rank BD only)")

v13_crvs <- zonator::read_curves("analyses/zonation/priocomp/13_load_abf_wgt_rwr_all/13_load_abf_wgt_rwr_all_out/13_load_abf_wgt_rwr_all.curves.txt")
v13_mean <- v13_crvs %>%
  dplyr::select(pr_lost, ave_pr, w_pr) %>%
  dplyr::mutate(variant = "RWR")

v14_crvs <- zonator::read_curves("analyses/zonation/priocomp/14_load_abf_wgt_ilp_all/14_load_abf_wgt_ilp_all_out/14_load_abf_wgt_ilp_all.curves.txt")
v14_mean <- v14_crvs %>%
  dplyr::select(pr_lost, ave_pr, w_pr) %>%
  dplyr::mutate(variant = "ILP")


# Combine data ------------------------------------------------------------

# Between data groups

v04_grp_mean$variant <- gsub("Full BD", "BD (rank ES+BD)",
                             v04_grp_mean$variant)
v04_grp_mean$variant <- gsub("Full ES", "ES (rank BD+ES)",
                             v04_grp_mean$variant)
datag_perf <- dplyr::bind_rows(list(v04_grp_mean, v06_mean, v08_mean,
                                        v11_mean, v12_mean))

# Between methods
metho_perf <- dplyr::bind_rows(list(dplyr::filter(v04_mean,
                                                      variant == "ZON"),
                                        v13_mean, v14_mean))

# Plot mean curves --------------------------------------------------------

x_lab <- "\nFraction of landscape"
x_scale <- scale_x_continuous(breaks = seq(0, 1, 0.2),
                              labels = paste(100 * seq(1, 0, -0.2), "%"))
y_scale <- scale_y_continuous(breaks = seq(0, 1, 0.2),
                              labels = paste(100 * seq(0, 1, 0.2), "%"))

p1 <- ggplot2::ggplot(datag_perf, aes(x = pr_lost, y = ave_pr,
                                          color = variant,
                                          linetype = variant)) +
  geom_line(size = 0.9) + xlab(x_lab) +
  scale_linetype_manual("", values = rep(1:3, 2)) +
  scale_color_manual("", values =  rev(rep(viridis(2, end = 0.7), 3, each = 3))) +
  x_scale + y_scale + ylab("Average proportion of features covered\n") +
  ggtitle("A") + theme_minimal() +
  theme(legend.position = c(0.25, 0.075),
        legend.justification = c(0.5, 0),
        legend.key.width = unit(1,"cm"))

p2 <- ggplot2::ggplot(metho_perf, aes(x = pr_lost, y = w_pr,
                                        color = variant)) +
  geom_line(size = 1) + x_scale + y_scale + xlab(x_lab) + ylab("") +
  ggtitle("B") +  scale_color_viridis("", discrete = TRUE) +
  theme_minimal() +
  theme(legend.position = c(0.25, 0.21),
        legend.justification = c(0.5, 0),
        legend.key.width = unit(1,"cm"))

fig3 <- gridExtra::grid.arrange(p1, p2, ncol = 2, nrow = 1)

# Save Figure -------------------------------------------------------------

ggsave("reports/figures/04_figure_03.png", fig3)

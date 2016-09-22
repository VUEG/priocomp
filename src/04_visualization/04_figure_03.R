library(dplyr)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(zonator)
library(viridis)

zproject <- zonator::load_zproject('analyses/zonation/priocomp/')


# Load variants and get curves data ---------------------------------------

v04_abf_all_wgt <- zonator::get_variant(zproject, 4)
v06_abf_es <- zonator::get_variant(zproject, 6)
v08_abf_bd <- zonator::get_variant(zproject, 8)

groupnames(v04_abf_all_wgt) <- c("1" = "ES", "2" = "BD")
v04_crvs <- zonator::curves(v04_abf_all_wgt, groups = FALSE)
v04_grp_crvs <- zonator::curves(v04_abf_all_wgt, groups = TRUE)

v06_crvs <- zonator::curves(v06_abf_es, groups = FALSE)

v08_crvs <- zonator::curves(v08_abf_bd, groups = FALSE)


# Re-arrange data ---------------------------------------------------------

# Get the average proportion remaining over all features in respective variants
v04_mean <- v04_crvs %>%
  dplyr::select(pr_lost, ave_pr, w_pr) %>%
  dplyr::mutate(variant = "ALL")

v06_mean <- v06_crvs %>%
  dplyr::select(pr_lost, ave_pr, w_pr) %>%
  dplyr::mutate(variant = "ES")

v08_mean <- v08_crvs %>%
  dplyr::select(pr_lost, ave_pr, w_pr) %>%
  dplyr::mutate(variant = "BD")

# Combine all the average data
mean_perf <- dplyr::bind_rows(list(v04_mean, v06_mean, v08_mean))

# Pre-load ----------------------------------------------------------------

v09_crvs <- zonator::read_curves("analyses/zonation/priocomp/09_load_es/09_load_es_out/09_load_es.curves.txt")
v09_mean <- v09_crvs %>%
  dplyr::select(pr_lost, ave_pr, w_pr) %>%
  dplyr::mutate(variant = "PES")

v10_crvs <- zonator::read_curves("analyses/zonation/priocomp/10_load_bd/10_load_bd_out/10_load_bd.curves.txt")
v10_mean <- v10_crvs %>%
  dplyr::select(pr_lost, ave_pr, w_pr) %>%
  dplyr::mutate(variant = "PBD")

mean_pl_perf <- dplyr::bind_rows(list(v04_mean, v09_mean, v10_mean))

# Plot mean curves --------------------------------------------------------

p1 <- ggplot2::ggplot(mean_perf, aes(x = pr_lost, y = ave_pr,
                                    color = variant)) +
  geom_line(size = 1) + xlab("") +
  ylab("Average feature remaining\n") + ggtitle("A") +
  scale_color_viridis("", discrete = TRUE) +
  theme_minimal()

p2 <- ggplot2::ggplot(mean_pl_perf, aes(x = pr_lost, y = ave_pr,
                                        color = variant)) +
  geom_line(size = 1) + xlab("\nLandscape lost") +
  ylab("Average feature remaining\n") + ggtitle("B") +
  scale_color_viridis("", discrete = TRUE) + theme_minimal()

p3 <- ggplot2::ggplot(mean_pl_perf, aes(x = pr_lost, y = w_pr,
                                        color = variant)) +
  geom_line(size = 1) + xlab("\nLandscape lost") +
  ylab("Average feature remaining\n") + ggtitle("C") +
  scale_color_viridis("", discrete = TRUE) + theme_minimal()

fig3 <- gridExtra::grid.arrange(p1, p2, p3, ncol = 1, nrow = 3)

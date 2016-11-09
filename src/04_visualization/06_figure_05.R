library(dplyr)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(tidyr)
library(zonator)
library(viridis)

# Helper functions --------------------------------------------------------

get_perf_level <- function(data, x) {
  perf_data <- data %>%
    dplyr::group_by(method, group, feature) %>%
    dplyr::filter(pr_lost >= x) %>%
    dplyr::arrange(pr_lost) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(level = x)
  return(perf_data)
}

# Load variants and get curves data ---------------------------------------

zproject <- zonator::load_zproject('analyses/zonation/priocomp/')
v04_abf_all_wgt <- zonator::get_variant(zproject, 4)

v04_crvs <- zonator::curves(v04_abf_all_wgt, groups = FALSE)
v13_crvs <- zonator::read_curves("analyses/zonation/priocomp/13_load_abf_wgt_rwr_all/13_load_abf_wgt_rwr_all_out/13_load_abf_wgt_rwr_all.curves.txt")
v14_crvs <- zonator::read_curves("analyses/zonation/priocomp/14_load_abf_wgt_ilp_all/14_load_abf_wgt_ilp_all_out/14_load_abf_wgt_ilp_all.curves.txt")

# Mean curves data --------------------------------------------------------

v04_mean <- v04_crvs %>%
  dplyr::select(pr_lost, ave_pr, w_pr) %>%
  dplyr::mutate(variant = "ZON_ALL")

# FIXME! For now (2016-11-08), calculate means manually. Groups file
# definition were not correct in the last run.
v04_es_mean <- v04_crvs %>%
  dplyr::select(woodprod_average:species_richness_vascular_plants) %>%
  dplyr::mutate(ave_pr = rowMeans(., na.rm = TRUE)) %>%
  dplyr::select(ave_pr) %>%
  dplyr::bind_cols(dplyr::select(v04_crvs, pr_lost), .) %>%
  dplyr::mutate(variant = "ZON_ES")

v04_bd_mean <- v04_crvs %>%
  dplyr::select(-(pr_lost:species_richness_vascular_plants)) %>%
  dplyr::mutate(ave_pr = rowMeans(., na.rm = TRUE)) %>%
  dplyr::select(ave_pr) %>%
  dplyr::bind_cols(dplyr::select(v04_crvs, pr_lost), .) %>%
  dplyr::mutate(variant = "ZON_BD")

v04_grp_mean <- dplyr::bind_rows(v04_es_mean, v04_bd_mean)

v13_crvs <- zonator::read_curves("analyses/zonation/priocomp/13_load_abf_wgt_rwr_all/13_load_abf_wgt_rwr_all_out/13_load_abf_wgt_rwr_all.curves.txt")
v13_mean <- v13_crvs %>%
  dplyr::select(pr_lost, ave_pr, w_pr) %>%
  dplyr::mutate(variant = "RWR_ALL")

# FIXME! For now (2016-11-08), calculate means manually. Groups file
# definition were not correct in the last run.
v13_es_mean <- v13_crvs %>%
  dplyr::select(f1:f9) %>%
  dplyr::mutate(ave_pr = rowMeans(., na.rm = TRUE)) %>%
  dplyr::select(ave_pr) %>%
  dplyr::bind_cols(dplyr::select(v13_crvs, pr_lost), .) %>%
  dplyr::mutate(variant = "RWR_ES")

v13_bd_mean <- v13_crvs %>%
  dplyr::select(-(pr_lost:f9)) %>%
  dplyr::mutate(ave_pr = rowMeans(., na.rm = TRUE)) %>%
  dplyr::select(ave_pr) %>%
  dplyr::bind_cols(dplyr::select(v13_crvs, pr_lost), .) %>%
  dplyr::mutate(variant = "RWR_BD")

v13_grp_mean <- dplyr::bind_rows(v13_es_mean, v13_bd_mean)

v14_crvs <- zonator::read_curves("analyses/zonation/priocomp/14_load_abf_wgt_ilp_all/14_load_abf_wgt_ilp_all_out/14_load_abf_wgt_ilp_all.curves.txt")
v14_mean <- v14_crvs %>%
  dplyr::select(pr_lost, ave_pr, w_pr) %>%
  dplyr::mutate(variant = "ILP_ALL")

# FIXME! For now (2016-11-08), calculate means manually. Groups file
# definition were not correct in the last run.
v14_es_mean <- v14_crvs %>%
  dplyr::select(f1:f9) %>%
  dplyr::mutate(ave_pr = rowMeans(., na.rm = TRUE)) %>%
  dplyr::select(ave_pr) %>%
  dplyr::bind_cols(dplyr::select(v14_crvs, pr_lost), .) %>%
  dplyr::mutate(variant = "ILP_ES")

v14_bd_mean <- v14_crvs %>%
  dplyr::select(-(pr_lost:f9)) %>%
  dplyr::mutate(ave_pr = rowMeans(., na.rm = TRUE)) %>%
  dplyr::select(ave_pr) %>%
  dplyr::bind_cols(dplyr::select(v14_crvs, pr_lost), .) %>%
  dplyr::mutate(variant = "ILP_BD")

v14_grp_mean <- dplyr::bind_rows(v14_es_mean, v14_bd_mean)

# Between methods
method_perf <- dplyr::bind_rows(v04_mean, v04_grp_mean, v13_mean,
                                v13_grp_mean, v14_mean, v14_grp_mean)

x_lab <- "\nFraction of the landscape"
x_scale <- scale_x_continuous(breaks = seq(0, 1, 0.2),
                              labels = paste(100 * seq(1, 0, -0.2), "%"))
y_scale <- scale_y_continuous(breaks = seq(0, 1, 0.2),
                              labels = paste(100 * seq(0, 1, 0.2), "%"))
# Define x-axis vlines that are used to link to Fig 5
vlines_x <- c(0.75, 0.9, 0.98)
vlines_labs <- c("25%", "10%", "2%")
vlines_labs_x <- vlines_x + 0.03
vlines_labs_y <- 1.0

p1 <- ggplot2::ggplot(method_perf, aes(x = pr_lost, y = ave_pr,
                                       color = variant,
                                       linetype = variant)) +
  geom_vline(xintercept = vlines_x, alpha = 0.5, linetype = 3) +
  annotate("text", x = vlines_labs_x, y = vlines_labs_y,
           label = vlines_labs, size = 3) +
  geom_line(size = 1) + x_scale + y_scale + xlab(x_lab) + ylab("") +
  scale_linetype_manual("", values = rep(1:3, 3)) +
  scale_color_manual("", values =  rev(rep(viridis(3, end = 0.9), 1, each = 3))) +
  ylab("Average feature distribution covered\n") + ggtitle("A") + theme_minimal() +
  theme(legend.position = c(0.15, 0.05),
        legend.justification = c(0.5, 0),
        legend.key.width = unit(1,"cm"))

# All curves data ---------------------------------------------------------

v04_crvs$method <- "ZON"
v13_crvs$method <- "RWR"
v14_crvs$method <- "ILP"

v04_mean <- v04_crvs %>%
  dplyr::select(pr_lost, ave_pr, w_pr) %>%
  dplyr::mutate(variant = "ZON_ALL")

# zonator::read_curves() doesn't produce feature names, but they are the same
# as in 04.
names(v13_crvs) <- names(v04_crvs)
names(v14_crvs) <- names(v04_crvs)

# ES features
all_crvs_es <- dplyr::bind_rows(v04_crvs, v13_crvs, v14_crvs) %>%
  dplyr::select(-(cost:ext2)) %>%
  dplyr::select(pr_lost:species_richness_vascular_plants, method) %>%
  tidyr::gather(feature, pr_rem, -pr_lost, -method) %>%
  dplyr::mutate(group = "ES") %>%
  dplyr::select(pr_lost, method, group, feature, pr_rem)

# BD features
all_crvs_bd <- dplyr::bind_rows(v04_crvs, v13_crvs, v14_crvs) %>%
  dplyr::select(-(cost:species_richness_vascular_plants)) %>%
  tidyr::gather(feature, pr_rem, -pr_lost, -method) %>%
  dplyr::mutate(group = "BD") %>%
  dplyr::select(pr_lost, method, group, feature, pr_rem)

# ALL is effectively ES + BD
all_crvs <- dplyr::bind_rows(all_crvs_es, all_crvs_bd)
all_crvs_all <- all_crvs
all_crvs_all$group <- "ALL"
all_crvs <- dplyr::bind_rows(all_crvs, all_crvs_all)

perf_levels <- c(0.75, 0.9, 0.98)
perf_details <- dplyr::bind_rows(lapply(perf_levels,
                                        function(y) {get_perf_level(all_crvs, y)}))

perf_details$level <- factor(perf_details$level,
                             levels = perf_levels,
                             labels = paste0(100 - perf_levels * 100, "%"))

p2 <- ggplot(perf_details, aes(x = factor(group), y = pr_rem,
                                 fill = factor(method))) +
  scale_fill_manual("Method", values =  rev(rep(viridis(3, end = 0.9), 3, each = 1))) +
  geom_boxplot(outlier.colour = "darkgrey", outlier.size = 0.2) + facet_wrap(~level) +
  ylab("Feature distribution covered\n") + xlab("\nData group") +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = paste0(seq(0, 100, 25), "%")) +
  ggtitle("B") + theme_minimal()

# Combine panels
fig5 <- gridExtra::grid.arrange(p1, p2, ncol = 2, nrow = 1)

# Save Figure -------------------------------------------------------------

ggsave("reports/figures/12_figure_05.png", fig5, width = 12, height = 6)

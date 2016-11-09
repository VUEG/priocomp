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

v04_crvs <- zonator::curves(v04_abf_all_wgt, groups = FALSE)
v06_crvs <- zonator::curves(v06_abf_es, groups = FALSE)
v08_crvs <- zonator::curves(v08_abf_bd, groups = FALSE)

# Re-arrange data ---------------------------------------------------------

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

# Combine data ------------------------------------------------------------

# Between data groups

v04_grp_mean$variant <- gsub("ZON_BD", "BD (rank ES+BD)", v04_grp_mean$variant)
v04_grp_mean$variant <- gsub("ZON_ES", "ES (rank BD+ES)", v04_grp_mean$variant)

datag_perf <- dplyr::bind_rows(list(v04_grp_mean, v06_mean, v08_mean,
                                    v11_mean, v12_mean))

# Between methods
method_perf <- dplyr::bind_rows(list(dplyr::filter(v04_mean, variant == "ZON_ALL"),
                                     dplyr::filter(v04_grp_mean,
                                                   variant %in% c("BD (rank ES+BD)", "ES (rank BD+ES)")),
                                         v13_mean, v13_grp_mean, v14_mean, v14_grp_mean))
method_perf$variant <- gsub("BD\\ \\(rank\\ ES\\+BD\\)", "ZON_BD", method_perf$variant)
method_perf$variant <- gsub("ES\\ \\(rank\\ BD\\+ES\\)", "ZON_ES", method_perf$variant)

# Plot mean curves --------------------------------------------------------

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

p1 <- ggplot2::ggplot(datag_perf, aes(x = pr_lost, y = ave_pr,
                                          color = variant,
                                          linetype = variant)) +
  geom_vline(xintercept = vlines_x[2], alpha = 0.5, linetype = 3) +
  annotate("text", x = vlines_labs_x[2], y = vlines_labs_y,
           label = vlines_labs[2], size = 3) +
  geom_line(size = 0.9) + xlab(x_lab) +
  scale_linetype_manual("", values = rep(1:3, 2)) +
  scale_color_manual("", values =  rev(rep(viridis(2, end = 0.7), 3, each = 3))) +
  x_scale + y_scale + ylab("Average feature distribution covered\n") +
  ggtitle("A") + theme_minimal() +
  theme(legend.position = c(0.25, 0.075),
        legend.justification = c(0.5, 0),
        legend.key.width = unit(1,"cm"))

p2 <- ggplot2::ggplot(method_perf, aes(x = pr_lost, y = ave_pr,
                                      color = variant,
                                      linetype = variant)) +
  geom_vline(xintercept = vlines_x, alpha = 0.5, linetype = 3) +
  annotate("text", x = vlines_labs_x, y = vlines_labs_y,
           label = vlines_labs, size = 3) +
  geom_line(size = 1) + x_scale + y_scale + xlab(x_lab) + ylab("") +
  scale_linetype_manual("", values = rep(1:3, 3)) +
  scale_color_manual("", values =  rev(rep(viridis(3, end = 0.9), 1, each = 3))) +
  ggtitle("B") + theme_minimal() +
  theme(legend.position = c(0.15, 0.05),
        legend.justification = c(0.5, 0),
        legend.key.width = unit(1,"cm"))

fig3 <- gridExtra::grid.arrange(p1, p2, ncol = 2, nrow = 1)

# Save Figure -------------------------------------------------------------

ggsave("reports/figures/11_figure_04.png", fig3, width = 12, height = 6)

# Extra -------------------------------------------------------------------

# Performance levels for the top 25%, 10% and 2% of the solution

get_perf_level <- function(data, variant_str, x) {
  if (!variant_str %in% unique(data$variant)) {
    stop("Bad variant name")
  }
  perf_data <- data %>%
    dplyr::filter(variant == variant_str) %>%
    dplyr::filter(pr_lost >= x) %>%
    dplyr::arrange(pr_lost) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(level = x)
  return(perf_data)
}

bd_only_top10 <- get_perf_level(v08_mean, "BD only", 0.9)
bdes_top10 <- get_perf_level(v04_grp_mean, "BD (rank ES+BD)", 0.9)
bd_esonly_top10 <- get_perf_level(v11_mean, "BD only (rank ES only)", 0.9)

es_only_top10 <- get_perf_level(v06_mean, "ES only", 0.9)
esbd_top10 <- get_perf_level(v04_grp_mean, "ES (rank BD+ES)", 0.9)
es_bdonly_top10 <- get_perf_level(v12_mean, "ES only (rank BD only)", 0.9)

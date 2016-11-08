library(dplyr)
library(ggplot2)
library(ggthemes)
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

# Re-arrange data ---------------------------------------------------------

v04_crvs$method <- "ZON"
v13_crvs$method <- "RWR"
v14_crvs$method <- "ILP"

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
                             labels = paste0("Top ", 100 - perf_levels * 100, "%"))

p1 <- ggplot(perf_details, aes(x = factor(group), y = pr_rem,
                                 fill = factor(method))) +
  scale_fill_manual("Method", values =  rev(rep(viridis(3, end = 0.9), 3, each = 1))) +
  geom_boxplot(outlier.colour = "darkgrey", outlier.size = 0.2) + facet_wrap(~level) +
  ylab("Average feature distribution covered\n") + xlab("\nData group") +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = paste0(seq(0, 100, 25), "%")) +
  theme_minimal()

# Save Figure -------------------------------------------------------------

ggsave("reports/figures/12_figure_05.png", p1, width = 12, height = 6)

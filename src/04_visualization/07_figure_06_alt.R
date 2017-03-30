library(dplyr)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(hrbrthemes)
library(tidyr)
library(raster)
library(zonator)
library(viridis)

# Helper functions --------------------------------------------------------

create_cost_plot <- function(x, title = NULL) {

  x_lab <- "\nFraction of the landscape"
  x_scale <- scale_x_continuous(breaks = seq(0, 1, 0.1),
                                labels = paste(100 * seq(0, 1, 0.1), "%"))
  y_scale <- scale_y_continuous(breaks = seq(0, 1, 0.2))
  # Define x-axis vlines that are used to link to Fig 5
  vlines_x <- c(0.9)
  vlines_labs <- c("Top 10%")
  vlines_labs_x <- vlines_x + 0.04
  vlines_labs_y <- 1.03

  p1 <- ggplot2::ggplot(x, aes(x = pr_lost, y = cost, color = variant,
                               linetype = cost_type)) +
    geom_vline(xintercept = vlines_x, alpha = 0.5, linetype = 3) +
    annotate("text", x = vlines_labs_x, y = vlines_labs_y,
             label = vlines_labs, size = 3) +
    geom_line(size = 0.8) + x_scale + y_scale + xlab(x_lab) + ylab("") +
    scale_color_manual("", values =  rev(viridis(3, end = 0.9))) +
    scale_linetype_manual("", values = c("solid", "dotted")) +
    ylab("Relateive cost of the solution\n") + ggtitle(title) +
    theme_ipsum_rc() +
    theme(legend.position = c(0.20, 0.35),
          legend.justification = c(0.5, 0),
          legend.key.width = unit(1,"cm"))
  return(p1)
}

create_perf_plot <- function(x, title = NULL) {

  x_lab <- "\nFraction of the landscape"
  x_scale <- scale_x_continuous(breaks = seq(0, 1, 0.2),
                                labels = paste(100 * seq(1, 0, -0.2), "%"))
  y_scale <- scale_y_continuous(breaks = seq(0, 1, 0.2),
                                labels = paste(100 * seq(0, 1, 0.2), "%"))
  # Define x-axis vlines that are used to link to Fig 5
  vlines_x <- c(0.9)
  vlines_labs <- c("10%")
  vlines_labs_x <- vlines_x + 0.02
  vlines_labs_y <- 1.03

  p1 <- ggplot2::ggplot(x, aes(x = pr_lost, y = ave_pr, color = method)) +
    geom_vline(xintercept = vlines_x, alpha = 0.5, linetype = 3) +
    facet_wrap(~ group) +
    annotate("text", x = vlines_labs_x, y = vlines_labs_y,
             label = vlines_labs, size = 3) +
    geom_line(size = 1) + x_scale + y_scale + xlab(x_lab) + ylab("") +
    scale_color_manual("", values =  rev(rep(viridis(3, end = 0.9), 3))) +
    ylab("Average feature distribution covered\n") + ggtitle(title) +
    theme_minimal() +
    theme(legend.position = c(0.20, 0.05),
          legend.justification = c(0.5, 0),
          legend.key.width = unit(1,"cm"))
  return(p1)
}

get_costs <- function(cost_raster, rank_raster, n = 1002,
                      rescale_cost = TRUE) {
  cost_raster_src <- raster::getValues(cost_raster)
  rank_raster_src <- raster::getValues(rank_raster)

  # Rank raster will always have <= informative cells
  mask <- is.na(rank_raster_src)
  # Get only non-NoData data
  cost_raster_src <- cost_raster_src[!mask]
  rank_raster_src <- rank_raster_src[!mask]

  # Make matrix
  cost_matrix <- cbind(rank_raster_src, cost_raster_src)
  # Order by rank
  cost_matrix <- cost_matrix[sort.list(cost_matrix[,1], decreasing = TRUE), ]
  # Instead of proportion remaining, use proportion lost
  cost_matrix[,1] <- 1.0 - cost_matrix[,1]
  # Some cost data can be NA, replace with 0
  cost_matrix[is.na(cost_matrix[,2]),2] <- 0
  # Calculate cumulative sum
  cost_matrix <- cbind(cost_matrix, cumsum(cost_matrix[,2]))
  # Select every nth row
  every <- nrow(cost_matrix) / n
  cost_matrix <- cost_matrix[seq(1, nrow(cost_matrix), every), ]
  # Coerce into a dataframe
  cost_df <- as.data.frame(cost_matrix)
  names(cost_df) <- c("pr_lost", "cost_per_cell", "cost")

  if (rescale_cost) {
    cost_df$cost <- cost_df$cost / max(cost_df$cost)
  }

  return(cost_df)
}

get_grp_means <- function(variant, groups) {
  zonator::groupnames(variant) <- groups
  grp_crvs <- zonator::curves(variant, groups = TRUE)

  # Use the first part of the first group label as a method label
  method_label <- unlist(strsplit(groups, "_"))[1]
  # Generate current group labels
  group_labels <- paste0("mean.", groups[1:2])

  # Get group means, must use SE with dplyr.
  # NOTE: mean for ALL needs to calculated from ES + BD because CST is
  # included in the non-grouped curves
  lng_grp_1 <- group_labels[1]
  lng_grp_2 <- group_labels[2]
  sht_grp_1 <- groups[[1]]
  sht_grp_2 <- groups[[2]]
  sht_grp_3 <- paste0(method_label, "_ALL")
  mutate_call <- lazyeval::interp(~ (a + b) / 2,
                                  a = as.name(sht_grp_1),
                                  b = as.name(sht_grp_2))
  grp_means <- grp_crvs %>%
    dplyr::select_("pr_lost", lng_grp_1, lng_grp_2) %>%
    dplyr::rename_(.dots = setNames(list(lng_grp_1, lng_grp_2),
                                    c(sht_grp_1, sht_grp_2))) %>%
    dplyr::mutate_(.dots = setNames(list(mutate_call), sht_grp_3)) %>%
    tidyr::gather_("variant", "ave_pr", c(sht_grp_1, sht_grp_2, sht_grp_3))

  return(grp_means)
}

get_perf_level <- function(data, x) {
  perf_data <- data %>%
    dplyr::group_by(variant, group, feature) %>%
    dplyr::filter(pr_lost >= x) %>%
    dplyr::arrange(pr_lost) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(level = x)
  return(perf_data)
}

# Load variants -----------------------------------------------------------

zproject <- zonator::load_zproject('analyses/zonation/priocomp/')

zon_groups <- c("1" = "ZON_ES", "2" = "ZON_BD", "3" = "ZON_CST")
rwr_groups <- c("1" = "RWR_ES", "2" = "RWR_BD", "3" = "RWR_CST")
ilp_groups <- c("1" = "ILP_ES", "2" = "ILP_BD", "3" = "ILP_CST")

v04_abf_all_wgt <- zonator::get_variant(zproject, 4)
v06_abf_all_wgt_cst <- zonator::get_variant(zproject, 6)

v23_load_rwr_all <- zonator::get_variant(zproject, 23)
v25_load_rwr_all_cst <- zonator::get_variant(zproject, 25)

v24_load_ilp_all <- zonator::get_variant(zproject, 24)
v26_load_ilp_all_cst <- zonator::get_variant(zproject, 26)

# Mean group curves data --------------------------------------------------

v04_grp_mean <- get_grp_means(v04_abf_all_wgt, zon_groups)
v06_grp_mean <- get_grp_means(v06_abf_all_wgt_cst, zon_groups)

v23_grp_mean <- get_grp_means(v23_load_rwr_all, rwr_groups)
v25_grp_mean <- get_grp_means(v25_load_rwr_all_cst, rwr_groups)

v24_grp_mean <- get_grp_means(v24_load_ilp_all, ilp_groups)
v26_grp_mean <- get_grp_means(v26_load_ilp_all_cst, ilp_groups)

# Combine performance data from different methods, for no-costs and costs
method_perf_nocost <- dplyr::bind_rows(v04_grp_mean, v23_grp_mean,
                                       v24_grp_mean)
method_perf_cost <- dplyr::bind_rows(v06_grp_mean, v25_grp_mean,
                                     v26_grp_mean)

# Cost data ---------------------------------------------------------------

# Read in the cost data rasters
cost_raster <- raster::raster("data/processed/features/eea/pop_density/pop_density_v5.tif")
# No costs
v04_rank_raster <- zonator::rank_raster(v04_abf_all_wgt)
v23_rank_raster <- zonator::rank_raster(v23_load_rwr_all)
v24_rank_raster <- zonator::rank_raster(v24_load_ilp_all)

v04_cost <- get_costs(cost_raster = cost_raster,
                      rank_raster = v04_rank_raster)
v04_cost$variant <- "ZON_ALL"
v23_cost <- get_costs(cost_raster = cost_raster,
                      rank_raster = v23_rank_raster)
v23_cost$variant <- "RWR_ALL"
v24_cost <- get_costs(cost_raster = cost_raster,
                      rank_raster = v24_rank_raster)
v24_cost$variant <- "ILP_ALL"

nocost_perf_cost <- dplyr::bind_rows(v04_cost, v23_cost, v24_cost)

# With costs
v06_rank_raster <- zonator::rank_raster(v06_abf_all_wgt_cst)
v25_rank_raster <- zonator::rank_raster(v25_load_rwr_all_cst)
v26_rank_raster <- zonator::rank_raster(v26_load_ilp_all_cst)

v06_cost <- get_costs(cost_raster = cost_raster,
                      rank_raster = v06_rank_raster)
v06_cost$variant <- "ZON_ALL"
v25_cost <- get_costs(cost_raster = cost_raster,
                      rank_raster = v25_rank_raster)
v25_cost$variant <- "RWR_ALL"
v26_cost <- get_costs(cost_raster = cost_raster,
                      rank_raster = v26_rank_raster)
v26_cost$variant <- "ILP_ALL"

cost_perf_cost <- dplyr::bind_rows(v06_cost, v25_cost, v26_cost)

# Combine no-cost and cost data
nocost_perf_cost$cost_type <- "Cost not included"
cost_perf_cost$cost_type <- "Cost included"
all_cost_perf_cost <- dplyr::bind_rows(nocost_perf_cost, cost_perf_cost)

# HACK
method_perf_cost <- method_perf_cost %>%
  tidyr::separate(variant, sep = "_", into = c("method", "group"))
method_perf_nocost <- method_perf_nocost %>%
  tidyr::separate(variant, sep = "_", into = c("method", "group"))

# Plot data ---------------------------------------------------------------

p1 <- create_perf_plot(method_perf_nocost, title = "A")
p2 <- create_perf_plot(method_perf_cost, title = "A")

#p3 <- create_cost_plot(cost_perf_cost, title = "B")
#p3_2 <- create_cost_plot(nocost_perf_cost, title = "B")
p3 <- create_cost_plot(all_cost_perf_cost)

ggsave("reports/figures/figure06/04_figure_06_just_costs.png",
       p3, width = 6, height = 6)

# All curves data ---------------------------------------------------------

v06_crvs <- zonator::curves(v06_abf_all_wgt_cst)
v25_crvs <- zonator::curves(v25_load_rwr_all_cst)
v26_crvs <- zonator::curves(v26_load_ilp_all_cst)

v06_crvs$variant <- "ZON"
v25_crvs$variant <- "RWR"
v26_crvs$variant <- "ILP"

# ES features
all_crvs_es <- dplyr::bind_rows(v06_crvs, v25_crvs, v26_crvs) %>%
  dplyr::select(-(cost:ext2)) %>%
  dplyr::select(pr_lost:species_richness_vascular_plants, variant) %>%
  tidyr::gather(feature, pr_rem, -pr_lost, -variant) %>%
  dplyr::mutate(group = "ES") %>%
  dplyr::select(pr_lost, variant, group, feature, pr_rem)

# BD features
all_crvs_bd <- dplyr::bind_rows(v06_crvs, v25_crvs, v26_crvs) %>%
  dplyr::select(-(cost:species_richness_vascular_plants)) %>%
  tidyr::gather(feature, pr_rem, -pr_lost, -variant) %>%
  dplyr::mutate(group = "BD") %>%
  dplyr::select(pr_lost, variant, group, feature, pr_rem)

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

perf_details <- perf_details %>%
  dplyr::filter(level == "10%")

p4 <- ggplot(perf_details, aes(x = factor(group), y = pr_rem,
                                 fill = factor(variant))) +
  scale_fill_manual("Method", values =  rev(rep(viridis(3, end = 0.9), 3, each = 1))) +
  geom_boxplot(outlier.colour = "darkgrey", outlier.size = 0.2) +
  ylab("Feature distribution covered\n") + xlab("\nData group") +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = paste0(seq(0, 100, 25), "%")) +
  theme_minimal()

ggsave("reports/figures/figure06/01_figure_06_boxplot.png", p4, width = 5, height = 5)

# Combine panels
fig5 <- gridExtra::grid.arrange(p2, p3, p4, ncol = 2, nrow = 2)

# Save Figure -------------------------------------------------------------

ggsave("reports/figures/figure06/01_figure_06.png", fig5, width = 9, height = 9)


#  Extra: cost distribution -----------------------------------------------

cost_raster_src <- raster::getValues(cost_raster)
cost_raster_src <- cost_raster_src[!is.na(cost_raster_src)]
# Rescale into [0, 1]
cost_raster_src <- cost_raster_src / max(cost_raster_src)

costs <- data.frame(value = cost_raster_src,
                    value_cat = cut(cost_raster_src, breaks = seq(0, 1, 0.1),
                                    include.lowest = TRUE))
costs$value_cat <-  factor(costs$value_cat,
                           levels = levels(costs$value_cat))

agg_costs <- costs %>%
  dplyr::group_by(value_cat) %>%
  dplyr::summarise(
    count = n(),
    total = sum(value)
  )

p5 <- ggplot2::ggplot(agg_costs, aes(x = value_cat, y = count)) +
  ggplot2::geom_bar(stat = "identity")

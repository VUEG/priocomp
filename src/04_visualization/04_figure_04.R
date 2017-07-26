# Figure 4:
#
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

create_boxplot <- function(x, title = "", draw_legend = TRUE,
                           draw_y_label = TRUE) {
  # ES features
  all_crvs_es <- x %>%
    dplyr::select(-(cost:ext2)) %>%
    dplyr::select(pr_lost:species_richness_vascular_plants, variant) %>%
    tidyr::gather(feature, pr_rem, -pr_lost, -variant) %>%
    dplyr::mutate(group = "ES") %>%
    dplyr::select(pr_lost, variant, group, feature, pr_rem)

  # BD features
  all_crvs_bd <- x %>%
    dplyr::select(-(cost:species_richness_vascular_plants)) %>%
    tidyr::gather(feature, pr_rem, -pr_lost, -variant) %>%
    dplyr::mutate(group = "BD") %>%
    dplyr::select(pr_lost, variant, group, feature, pr_rem)

  # ALL is effectively ES + BD
  all_crvs <- dplyr::bind_rows(all_crvs_es, all_crvs_bd)
  all_crvs_all <- all_crvs
  all_crvs_all$group <- "ALL"
  all_crvs <- dplyr::bind_rows(all_crvs, all_crvs_all)

  perf_levels <- c(0.9, 0.98)
  perf_details <- dplyr::bind_rows(lapply(perf_levels,
                                          function(y) {get_perf_level(all_crvs, y)}))

  perf_details$level <- factor(perf_details$level,
                               levels = perf_levels,
                               labels = paste0(100 - perf_levels * 100, "%"))

  p <- ggplot(perf_details, aes(x = factor(group), y = pr_rem,
                                 fill = factor(variant))) +
    scale_fill_manual("Method", values =  rev(rep(viridis(3, end = 0.9), 3, each = 1))) +
    geom_boxplot(outlier.colour = "darkgrey", outlier.size = 0.2, color = "darkgrey") +
    facet_wrap(~level) +
    xlab("\nData group") +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = paste0(seq(0, 100, 25), "%")) +
    ggtitle(title) + theme_minimal()
  if (draw_y_label) {
    p <- p + ylab("Feature distribution covered\n")
  } else {
    p <- p + ylab("")
  }
  if (!draw_legend) {
    p <- p + theme(legend.position = "none")
  }
  return(p)
}

create_cost_plot <- function(x, title = NULL, draw_legend = TRUE,
                             draw_y_label = TRUE) {

  x_lab <- "\nFraction of the landscape selected"
  x_scale <- scale_x_continuous(breaks = seq(0, 1, 0.2),
                                labels = paste(100 * seq(0, 1, 0.2), "%"))
  y_scale <- scale_y_continuous(breaks = seq(0, 1, 0.2))

  p1 <- ggplot2::ggplot(x, aes(x = pr_lost, y = cost, color = variant)) +
    geom_line(size = 1) + x_scale + y_scale + xlab(x_lab) + ylab("") +
    scale_color_manual("", values =  rev(viridis(3, end = 0.9))) +
    ggtitle(title) + theme_minimal() +
    theme(legend.key.width = unit(1,"cm"))
  if (draw_y_label) {
    p1 <- p1 + ylab("Cumulative cost of the solution\n")
  }
  if (!draw_legend) {
    p1 <- p1 + theme(legend.position = "none")
  }
  return(p1)
}

create_perf_plot <- function(x, title = NULL, draw_legend = TRUE,
                             draw_y_label = TRUE) {

  x_lab <- "\nFraction of the landscape selected"
  x_scale <- scale_x_continuous(breaks = seq(0, 1, 0.2),
                                labels = paste(100 * seq(1, 0, -0.2), "%"))
  y_scale <- scale_y_continuous(breaks = seq(0, 1, 0.2),
                                labels = paste(100 * seq(0, 1, 0.2), "%"))
  # Define x-axis vlines that are used to link to Fig 5
  vlines_x <- c(0.9, 0.98)
  vlines_labs <- c("10%", "2%")
  vlines_labs_x <- vlines_x + 0.02
  vlines_labs_y <- 1.03

  p1 <- ggplot2::ggplot(x, aes(x = pr_lost, y = ave_pr, color = variant,
                               linetype = variant)) +
    geom_vline(xintercept = vlines_x, alpha = 0.5, linetype = 3) +
    annotate("text", x = vlines_labs_x, y = vlines_labs_y,
             label = vlines_labs, size = 3) +
    geom_line(size = 1) + x_scale + y_scale + xlab(x_lab) +
    scale_linetype_manual("", values = rep(1:3, 3)) +
    scale_color_manual("", values =  rev(rep(viridis(3, end = 0.9), 1, each = 3))) +
    ggtitle(title) + theme_minimal() +
    theme(legend.key.width = unit(1,"cm"))
  if (draw_y_label) {
    p1 <- p1 + ylab("Average feature distribution covered\n")
  } else {
    p1 <- p1 + ylab("")
  }
  if (!draw_legend) {
    p1 <- p1 + theme(legend.position = "none")
  }
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

match_grob_width <- function(p1, p2) {
  # Get the ggplot grobs
  g1 <- ggplotGrob(p1)
  g2 <- ggplotGrob(p2)
  # Add a column to g1 and g2 (any width will do)
  #g1 <- gtable_add_cols(g1, unit(1, "mm"))
  # Set g5 width to be equal to g6 widths
  g1$widths = g2$widths
  return(g1)
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
v04_rank_raster <- zonator::rank_raster(v04_abf_all_wgt)
v06_rank_raster <- zonator::rank_raster(v06_abf_all_wgt_cst)
v23_rank_raster <- zonator::rank_raster(v23_load_rwr_all)
v24_rank_raster <- zonator::rank_raster(v24_load_ilp_all)
v25_rank_raster <- zonator::rank_raster(v25_load_rwr_all_cst)
v26_rank_raster <- zonator::rank_raster(v26_load_ilp_all_cst)

v04_cost <- get_costs(cost_raster = cost_raster,
                      rank_raster = v04_rank_raster)
v04_cost$variant <- "ZON_ALL"
v06_cost <- get_costs(cost_raster = cost_raster,
                      rank_raster = v06_rank_raster)
v06_cost$variant <- "ZON_ALL"

v23_cost <- get_costs(cost_raster = cost_raster,
                      rank_raster = v23_rank_raster)
v23_cost$variant <- "RWR_ALL"
v24_cost <- get_costs(cost_raster = cost_raster,
                      rank_raster = v24_rank_raster)
v24_cost$variant <- "ILP_ALL"
v25_cost <- get_costs(cost_raster = cost_raster,
                      rank_raster = v25_rank_raster)
v25_cost$variant <- "RWR_ALL"
v26_cost <- get_costs(cost_raster = cost_raster,
                      rank_raster = v26_rank_raster)
v26_cost$variant <- "ILP_ALL"

cost_perf_nocost <- dplyr::bind_rows(v04_cost, v23_cost, v24_cost)
cost_perf_cost <- dplyr::bind_rows(v06_cost, v25_cost, v26_cost)


# Tests -------------------------------------------------------------------


# cost_perf_cost %>%
#   dplyr::group_by(variant) %>%
#   dplyr::summarise(
#     min = min(cost),
#     mean = mean(cost),
#     max = max(cost)
#   )
#
# ggplot(cost_perf_cost, aes(x = pr_lost, y = log(cost), color = variant)) +
#   geom_line()

# Tests -------------------------------------------------------------------

# All curves data ---------------------------------------------------------

v04_crvs <- data.frame(zonator::curves(v04_abf_all_wgt))
v23_crvs <- data.frame(zonator::curves(v23_load_rwr_all))
v24_crvs <- data.frame(zonator::curves(v24_load_ilp_all))

v04_crvs$variant <- "ZON"
v23_crvs$variant <- "RWR"
v24_crvs$variant <- "ILP"

v06_crvs <- data.frame(zonator::curves(v06_abf_all_wgt_cst))
v25_crvs <- data.frame(zonator::curves(v25_load_rwr_all_cst))
v26_crvs <- data.frame(zonator::curves(v26_load_ilp_all_cst))

v06_crvs$variant <- "ZON"
v25_crvs$variant <- "RWR"
v26_crvs$variant <- "ILP"

all_dat_nocosts <- dplyr::bind_rows(v04_crvs, v23_crvs, v24_crvs)
all_dat_costs <- dplyr::bind_rows(v06_crvs, v25_crvs, v26_crvs)

# Plot data ---------------------------------------------------------------

p1 <- create_perf_plot(method_perf_nocost, title = "A", draw_legend = FALSE)
p2 <- create_perf_plot(method_perf_cost, title = "B", draw_y_label = FALSE)

p3 <- create_cost_plot(cost_perf_nocost, title = "E", draw_legend = FALSE)
p4 <- create_cost_plot(cost_perf_cost, title = "F", draw_y_label = FALSE)

p5 <- create_boxplot(all_dat_nocosts, "C", draw_legend = FALSE)
p6 <- create_boxplot(all_dat_costs, "D", draw_y_label = FALSE)

# Match subplot widths
p1 <- match_grob_width(p1, p2)
p3 <- match_grob_width(p3, p4)
p5 <- match_grob_width(p5, p6)

# Combine panels
fig5 <- gridExtra::grid.arrange(p1, p2, p5, p6, p3, p4, ncol = 2, nrow = 3)

# Save Figure -------------------------------------------------------------

ggsave("reports/figures/figure04/01_figure_04.png", fig5, width = 10, height = 8.5)


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

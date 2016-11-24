library(cowplot)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(readr)
library(scales)
library(tidyr)
library(viridis)

# Utility functions -------------------------------------------------------

# Extract given statistic. NOTE: no reality checking is done, use with care.
extract_stat <- function(x, stat) {
  requested_stat <- x %>%
    dplyr::select(f1_type, f1_method, f2_type, f2_method, tau, cmcs, jac_01, jac_09) %>%
    tidyr::gather(statistic, value, -f1_type, -f1_method, -f2_type, -f2_method) %>%
    dplyr::filter(statistic == stat) %>%
    dplyr::select(-statistic)
  return(requested_stat)
}

# Generate self-crossing comparison stats (i.e. 1.0).
generate_self_cross <- function(stat_name, stat_value = 1.0) {
  method_codes <- c("RWR", "ZON", "ILP")
  type_codes <- c("ALL", "ALL_WGT", "ES", "BD")

  keys <- c()
  f_methods <- c()
  f_types <- c()

  for (i in method_codes) {
    for (j in type_codes) {
      keys <- c(keys, paste(i, j, i, j, sep = "_"))
      f_methods <- c(f_methods, i)
      f_types <- c(f_types, j)
    }
  }
  df <- data_frame(key = keys, f1_method = f_methods, f1_type = f_types,
                   f2_method = f_methods, f2_type = f_types)
  for (name in stat_name) {
    df[,stat_name] <- stat_value
  }
  return(df)
}

# Custom grid arrangement with only one legend
grid_arrange_shared_legend <- function(..., ncol = length(list(...)),
                                       nrow = 1, position = c("bottom", "right", "left")) {

  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  # Insert legend into the list of grobs
  gll <- gl[1:5]
  gll[[6]] <- legend
  gll[7:9] <- gl[7:9]
  gl <- c(gl, ncol = ncol, nrow = nrow)
  combined <- do.call(arrangeGrob, gll)
  return(ggdraw(combined))
}


# Simplify a character string name of a method used so that a 3-character
# code is returned.
match_method <- Vectorize(
  function(x) {
    if (grepl(".+\\/zonation\\/.+", x)) {
      method <- "ZON"
    } else if (grepl(".+\\/RWR\\/.+", x)) {
      method <- "RWR"
    } else if (grepl(".+\\/ILP\\/.+", x)) {
      method <- "ILP"
    } else {
      stop("Method", x, " not matched")
    }
    return(method)
  }, c("x"), USE.NAMES = FALSE, SIMPLIFY = TRUE)

# Simplify a character string name of a data collection type used so that
# a 3-character code is returned.
match_type <- Vectorize(
  function(x) {
    if (grepl("_all.tif", x) | grepl("_all_stats", x) | grepl("02_abf", x)) {
      mtype <- "ALL"
    } else if (grepl("_all_weights.tif", x) | grepl("_all_weights_stats", x) | grepl("wgt", x)) {
      mtype <- "ALL_WGT"
    } else if (grepl("_es", x)) {
      mtype <- "ES"
    } else if (grepl("_bd", x)) {
      mtype <- "BD"
    } else {
      stop("Type ", x, " not matched")
    }
    return(mtype)
  }, c("x"), USE.NAMES = FALSE, SIMPLIFY = TRUE)


# Plot a cross-comparison stat heatmap
plot_stat <- function(x, title, ...) {

  rwr_rwr_stat <- x %>%
    filter(f1_method == "RWR" & f2_method == "RWR")
  rwr_zon_stat <- x %>%
    filter(f1_method == "RWR" & f2_method == "ZON")
  rwr_ilp_stat <- x %>%
    filter(f1_method == "RWR" & f2_method == "ILP")
  zon_zon_stat <- x %>%
    filter(f1_method == "ZON" & f2_method == "ZON")
  zon_ilp_stat <- x %>%
    filter(f1_method == "ZON" & f2_method == "ILP")
  ilp_ilp_stat <- x %>%
    filter(f1_method == "ILP" & f2_method == "ILP")

  create_empty_subplot <- function(xx) {
    sub_p <- ggplot(xx , aes(x = f1_type, y = f2_type, fill = value)) +
      geom_blank() +
      coord_equal() + coord_flip() +
      labs(x = NULL, y = NULL, title = "") +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank())
    return(sub_p)
  }

  create_subplot <- function(xx, min_lim = 0.0, max_lim = 1.0, step = 0.25,
                             axis_titles = TRUE, axis_text = TRUE,
                             margins = unit(c(1, 1, 1, 1), "mm")) {

    sub_p <- ggplot(xx , aes(x = f1_type, y = f2_type, fill = value)) +
      geom_tile(color = "white", size = 0.1) +
      geom_text(aes(label = sprintf("%0.2f", round(value, digits = 2))),
                color = "black", size = 4) +
      scale_fill_viridis(name = title, label = comma,
                         limits = c(min_lim, max_lim),
                         breaks = seq(min_lim, max_lim, by = step)) +
      coord_equal() +
      labs(x = unique(xx$f1_method), y = unique(xx$f2_method)) +
      theme_tufte(base_family = "Helvetica") +
      theme(title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_text(size = 7),
            legend.text = element_text(size = 12),
            plot.margin = margins)

    if (!axis_titles) {
      sub_p <- sub_p + theme(axis.title = element_blank())
    }
    if (!axis_text) {
      sub_p <- sub_p + theme(axis.text = element_blank())
    }

    return(sub_p)
  }

  rwr_rwr_p <- create_subplot(rwr_rwr_stat, axis_titles = FALSE,
                              axis_text = FALSE, ...)
  rwr_zon_p <- create_subplot(rwr_zon_stat, axis_titles = FALSE,
                              axis_text = FALSE, ...)
  rwr_ilp_p <- create_subplot(rwr_ilp_stat, axis_titles = FALSE,
                              axis_text = FALSE, ...)
  # Make a blank panel to set the layout correctly
  zon_rwr_p <- create_empty_subplot(zon_zon_stat)
  zon_zon_p <- create_subplot(zon_zon_stat, axis_titles = FALSE,
                              axis_text = FALSE, ...)
  zon_ilp_p <- create_subplot(zon_ilp_stat, axis_titles = FALSE,
                              axis_text = FALSE, ...)
  ilp_rwr_p <- create_empty_subplot(ilp_ilp_stat)
  ilp_zon_p <- create_empty_subplot(ilp_ilp_stat)
  ilp_ilp_p <- create_subplot(ilp_ilp_stat, axis_titles = FALSE,
                              axis_text = FALSE, ...)

  p <- grid_arrange_shared_legend(rwr_ilp_p, rwr_zon_p, rwr_rwr_p,
                                  zon_ilp_p, zon_zon_p, zon_rwr_p,
                                  ilp_ilp_p, ilp_zon_p, ilp_rwr_p,
                                  ncol = 3, nrow = 3, position = "right")
  return(p)
}

# Read in and arrange the data ---------------------------------------------

## Kendall tau rank correlation

cors <- readr::read_csv("analyses/comparison/cross_correlation.csv") %>%
  dplyr::mutate(f1_method = match_method(feature1), f1_type = match_type(feature1),
                f2_method = match_method(feature2), f2_type = match_type(feature2),
                key = paste(f1_method, f1_type, f2_method, f2_type,
                            sep = "_")) %>%
  dplyr::select(key, f1_method, f1_type, f2_method, f2_type, tau) %>%
  # Manully fill in the diagonal values (i.e. self-cross, 1.0) since they
  # have not been included in the data
  dplyr::bind_rows(., generate_self_cross("tau", 1.0)) %>%
  dplyr::arrange(key)

## Map comparison statistic

# NOTE: use the complement of MCS: CMCS = 1 - MCS
mcss <- readr::read_csv("analyses/comparison/cross_mcs.csv") %>%
  dplyr::mutate(f1_method = match_method(feature1), f1_type = match_type(feature1),
                f2_method = match_method(feature2), f2_type = match_type(feature2),
                key = paste(f1_method, f1_type, f2_method, f2_type,
                            sep = "_"), cmcs = 1 - mcs) %>%
  dplyr::select(key, f1_method, f1_type, f2_method, f2_type, cmcs) %>%
  dplyr::bind_rows(., generate_self_cross("cmcs", 1.0)) %>%
  dplyr::arrange(key)

## Jaccard coefficients for different thresholds

jac <- readr::read_csv("analyses/comparison/cross_jaccard.csv") %>%
  dplyr::mutate(f1_method = match_method(feature1), f1_type = match_type(feature1),
                f2_method = match_method(feature2), f2_type = match_type(feature2),
                key = paste(f1_method, f1_type, f2_method, f2_type,
                            sep = "_")) %>%
  dplyr::mutate(threshold = gsub("\\(0.0, 0.1, 0.0, 0.1\\)", 0.10, threshold)) %>%
  dplyr::mutate(threshold = gsub("\\(0.9, 1.0, 0.9, 1.0\\)", 0.90, threshold)) %>%
  #dplyr::filter(threshold == 0.10 | threshold == 0.90) %>%
  dplyr::select(key, f1_method, f1_type, f2_method, f2_type, threshold, coef) %>%
  dplyr::mutate(threshold = paste0("jac_", gsub("\\.", "", threshold))) %>%
  tidyr::spread(threshold, coef) %>%
  dplyr::bind_rows(., generate_self_cross(c("jac_01", "jac_09"), 1.0)) %>%
  dplyr::arrange(key)

## Join all stats and do additional filtering

# Join correlation and map comparison statistics
all_stats <- dplyr::left_join(cors, mcss, by = c("key" = "key")) %>%
  dplyr::select(key, f1_method = f1_method.x, f2_method = f2_method.x,
                f1_type = f1_type.x, f2_type = f2_type.x, tau, cmcs)
# Join in also the jaccard coefficients
all_stats <- dplyr::left_join(all_stats, jac, by = c("key" = "key")) %>%
  dplyr::select(f1_method = f1_method.x, f1_type = f1_type.x,
                f2_method = f2_method.x, f2_type = f2_type.x, tau, cmcs, jac_01, jac_09)
# Remove ALL and rename ALL_WGT to ALL. From this point on, "ALL" means all
# features with weights. In same go, make f1_type and f2_type ordered
# factors
all_stats <- all_stats %>%
  # Remove original "ALL" types
  dplyr::filter(f1_type != "ALL") %>%
  dplyr::filter(f2_type != "ALL") %>%
  # Rename original "ALL_WGT" to "ALL"
  dplyr::mutate(f1_type = gsub("ALL_WGT", "ALL", f1_type),
                f2_type = gsub("ALL_WGT", "ALL", f2_type)) %>%
  # Convert f1_type and f2_type to factors
  dplyr::mutate(f1_type = factor(f1_type, levels = c("ALL", "ES", "BD"), ordered = TRUE),
                f2_type = factor(f2_type, levels = c("ALL", "ES", "BD"), ordered = TRUE)) %>%
  dplyr::arrange(f1_method, f1_type, f2_method, f2_type)

# Create the plot ---------------------------------------------------------

tau <- extract_stat(all_stats, "tau")
jac_01 <- extract_stat(all_stats, "jac_01")
jac_09 <- extract_stat(all_stats, "jac_09")
cmcs <- extract_stat(all_stats, "cmcs")

p1 <- plot_stat(tau, title = "COR", min_lim = -0.25, max_lim = 1.0, step = 0.25)
p2 <- plot_stat(cmcs, title = "MCS")
p3 <- plot_stat(jac_09, title = "top10")
p4 <- plot_stat(jac_01, title = "low10")

#p_combined <- grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)

# Save plots --------------------------------------------------------------

img_width <- 7
img_width <- 6.6

ggsave("reports/figures/09_figure_04_A.png", p1, width = img_width, height = img_width)
ggsave("reports/figures/10_figure_04_B.png", p2, width = img_width, height = img_width)
ggsave("reports/figures/11_figure_04_C.png", p3, width = img_width, height = img_width)
ggsave("reports/figures/12_figure_04_D.png", p4, width = img_width, height = img_width)


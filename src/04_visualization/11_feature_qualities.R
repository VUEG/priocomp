library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(zonator)


# Helper functions --------------------------------------------------------

# panel.smooth function is built in.
# panel.cor puts correlation in upper panels, size proportional to correlation
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  if (missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


plot_pairs <- function(dat, taxon_name = NULL) {

  if (is.null(taxon_name)) {
    title <- "All taxa"
    dat <- dat %>%
      dplyr::mutate(log_count = log(count),
                    log_mean_ol = log(mean_ol)) %>%
      dplyr::select(log_count, log_mean_ol, morans_i, pr_rem)
  } else {
    title <- taxon_name
  }
  pairs(dat, lower.panel = panel.smooth, upper.panel = panel.cor,
        pch = 20, main = title)

}

# Load variants and configure groups --------------------------------------

zproject <- zonator::load_zproject('analyses/zonation/priocomp/')
# 04_abf_all_wgt

# Get performance for top 5% and Drop first (pr_lost) and last (cost)
# columns
v04_top5 <- zonator::get_variant(zproject, 4) %>%
  zonator::results() %>%
  zonator::performance(pr.lost = 0.95) %>%
  dplyr::select(-pr_lost, -pop_density_v5) %>%
  tidyr::gather(feature, pr_rem)

# Load auxiliary data -----------------------------------------------------

feature_ranges <- readr::read_csv("data/feature_ranges.csv") %>%
  dplyr::select(feature_long = feature, dplyr::everything()) %>%
  dplyr::mutate(feature = gsub("\\.tif$", "", basename(feature_long))) %>%
  dplyr::select(feature, dplyr::everything(), -feature_long)

feature_autocor <- readr::read_csv("data/morans_I_values_772_features_2017-03-08_04-08-05.csv") %>%
  dplyr::mutate(feature = gsub("\\.tif$", "", basename(feature)))

feature_data <- feature_ranges %>%
  dplyr::left_join(feature_autocor) %>%
  dplyr::left_join(v04_top5)


# Plot cross-correlation matrix -------------------------------------------

plot_pairs(feature_data)

library(yaml)
library(zonator)

source("~/Dropbox/project.J-PriPA/japan-zsetup/R/00_lib/utils.R")

# Setup -------------------------------------------------------------------

# Load data manifest
data_manifest <- yaml::yaml.load_file("data/data_manifest.yml")[[1]]
pro_data <- data_manifest[["https://beehub.nl/environmental-geography-group"]][[2]][[2]]
pro_feat <- sapply(1:length(pro_data), function(x) pro_data[[x]][[1]]$resources)

udr_data <- data_manifest[["https://beehub.nl/environmental-geography-group"]][[5]][[2]]
amp_feat <- udr_data[[1]]$european_tetrapods[[1]]$amphibians$resources
bir_feat <- udr_data[[1]]$european_tetrapods[[2]]$birds$resources
mam_feat <- udr_data[[1]]$european_tetrapods[[3]]$mammals$resources
rep_feat <- udr_data[[1]]$european_tetrapods[[4]]$reptiles$resources

# Helper functions --------------------------------------------------------

get_group_stats <- function(x, variant_id) {

  # Add one to ES (wood production)
  new_groups <- c(rep(1, length(pro_feat) + 1),
                  rep(2, length(amp_feat)), rep(3, length(bir_feat)),
                  rep(4, length(mam_feat)), rep(5, length(rep_feat)),
                  6)
  new_group_names <- c("1" = "ES", "2" = "Amphibians", "3" = "Birds",
                       "4" = "Mammals", "5" = "Reptiles", "6" = "Cost")

  x_variant <- zonator::get_variant(x, variant_id)
  zonator::groups(x_variant) <- new_groups
  zonator::groupnames(x_variant) <- new_group_names

  x_variant_cstat <- x_variant %>%
    calc_row_stats(groups = TRUE) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(variant = factor(variant, levels = unique(variant),
                                   ordered = TRUE))
  return(x_variant_cstat)
}

# Load variants and configure groups --------------------------------------

zproject <- zonator::load_zproject('analyses/zonation/priocomp/')

highlights <- c(0.1, 0.25)

labels <- c("Ecosystem services", "Amphibian", "Birds", "Mammals",
            "Reptiles")

# 04_abf_all_wgt
abf_all_wgt_cstat <- get_group_stats(zproject, 4)
# 06_abf_all_wgt_cst
abf_all_wgt_cst_cstat <- get_group_stats(zproject, 6)

p1 <- plot_curves(abf_all_wgt_cstat, title = "", non_param = TRUE,
                  invert_x = TRUE, nrow = 3, ncol = 2,
                  highlights = highlights, labels = labels)

p2 <- plot_curves(abf_all_wgt_cst_cstat, title = "", non_param = TRUE,
                  invert_x = TRUE, nrow = 3, ncol = 2,
                  highlights = highlights, labels = labels)

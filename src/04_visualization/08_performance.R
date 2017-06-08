library(dplyr)
library(DT)
library(tidyr)
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

get_groups <- function(group="all") {

  if (group == "all") {
    # Add one to ES (wood production)
    grps <- c(rep(1, length(pro_feat) + 1),
                    rep(2, length(amp_feat)), rep(3, length(bir_feat)),
                    rep(4, length(mam_feat)), rep(5, length(rep_feat)),
                    6)
  } else if (group == "es") {
    grps <- rep(1, length(pro_feat) + 1)
  } else if (group == "bd") {
    grps <- c(rep(2, length(amp_feat)), rep(3, length(bir_feat)),
              rep(4, length(mam_feat)), rep(5, length(rep_feat)))
  } else {
    stop("Unknown group: ", group)
  }
  return(grps)
}

get_group_names <- function() {
  return(c("1" = "ES", "2" = "Amphibians", "3" = "Birds",
           "4" = "Mammals", "5" = "Reptiles", "6" = "Cost"))
}

get_group_stats <- function(x, variant_id) {

  new_groups <- get_groups("all")
  new_group_names <- get_group_names()

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

get_bd_feature_perf <- function(variant) {

  bd_groups <- get_groups("bd")
  group_names <- get_group_names()
  bd_sum <- zonator::curves(variant) %>%
    dplyr::select(alytes_cisternasii:zamenis_situla) %>%
    dplyr::summarise_each(funs(sum)) %>%
    tidyr::gather(feature, pr_rem_sum) %>%
    dplyr::mutate(taxon = factor(group_names[bd_groups])) %>%
    dplyr::arrange(desc(pr_rem_sum))

  top10 <- bd_sum %>%
    dplyr::top_n(10, wt = pr_rem_sum) %>%
    dplyr::mutate(rank = 1:n()) %>%
    dplyr::select(rank, feature, taxon, pr_rem_sum)

  low10 <- bd_sum %>%
    dplyr::top_n(-10, wt = pr_rem_sum) %>%
    dplyr::mutate(rank = nrow(bd_sum) - (n():1 - 1)) %>%
    dplyr::select(rank, feature, taxon, pr_rem_sum)

  return(list("all_features" = bd_sum,
              "top10" = top10,
              "low10" = low10))

}

get_es_feature_perf <- function(variant) {
  es_sum <- zonator::curves(variant) %>%
    dplyr::select(woodprod_average:species_richness_vascular_plants) %>%
    dplyr::summarise_each(funs(sum)) %>%
    tidyr::gather(feature, pr_rem_sum) %>%
    dplyr::arrange(desc(pr_rem_sum))
  return(es_sum)
}

# Load variants and configure groups --------------------------------------

zproject <- zonator::load_zproject('analyses/zonation/priocomp/')

highlights <- c(0.1)

labels <- c("Ecosystem services", "Amphibian", "Birds", "Mammals",
            "Reptiles")


# Over all group performance ----------------------------------------------

# 04_abf_all_wgt
zon_all_wgt_cstat <- get_group_stats(zproject, 4)
# 06_abf_all_wgt_cst
zon_all_wgt_cst_cstat <- get_group_stats(zproject, 6)

# 23_load_rwr_all
rwr_all_wgt_cstat <- get_group_stats(zproject, 23)
# 25_load_rwr_all_cst
rwr_all_wgt_cstat <- get_group_stats(zproject, 25)

# 24_load_ilp_all
ilp_all_wgt_cstat <- get_group_stats(zproject, 24)
# 26_load_ilp_all_cst
ilp_all_wgt_cstat <- get_group_stats(zproject, 26)


# No coverage -------------------------------------------------------------

# No costs, how many features are not covered by 10%?

get_zeros <- function(x, fraction) {
  no_cov <- x %>%
    dplyr::group_by(variant) %>%
    dplyr::select(pr_lost, zeros, variant) %>%
    dplyr::filter(pr_lost >= fraction) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(fraction = fraction) %>%
    dplyr::ungroup()
  return(no_cov)
}

rwr_no_cov_10 <- get_zeros(rwr_all_wgt_cstat, 0.9)
zon_no_cov_10 <- get_zeros(zon_all_wgt_cstat, 0.9)
ilp_no_cov_10 <- get_zeros(ilp_all_wgt_cstat, 0.9)

no_cov_10 <- dplyr::bind_rows(rwr_no_cov_10, zon_no_cov_10, ilp_no_cov_10)

no_cov_10$method <- ifelse(grepl("^04", no_cov_10$variant), "ZON",
                           ifelse(grepl("^25", no_cov_10$variant), "RWR", "ILP"))
group_names <- get_group_names()
no_cov_10$group <- sapply(no_cov_10$variant, function(x) {
                                                    items <- unlist(strsplit(x, "_"))
                                                    group_no <- items[length(items)]
                                                    return(group_names[group_no])
                                                  })
no_cov_10 <- no_cov_10 %>%
  dplyr::select(method, group, zeros) %>%
  tidyr::spread(group, zeros) %>%
  dplyr::select(method, ES, Amphibians, Birds, Mammals, Reptiles) %>%
  dplyr::slice(match(c("RWR", "ZON", "ILP"), method))

# Feature-specific performance --------------------------------------------

# RWR
rwr_all_nc <- zonator::get_variant(zproject, 23)
## ES
rwr_all_nc_es_sum <- get_es_feature_perf(rwr_all_nc)
## BD
rwr_all_nc_bd <- get_bd_feature_perf(rwr_all_nc)

# ZON
zon_all_nc <- zonator::get_variant(zproject, 4)
## ES
zon_all_nc_es_sum <- get_es_feature_perf(zon_all_nc)
## BD
zon_all_nc_bd <- get_bd_feature_perf(zon_all_nc)

# ILP
ilp_all_nc <- zonator::get_variant(zproject, 24)
## ES
ilp_all_nc_es_sum <- get_es_feature_perf(ilp_all_nc)
## BD
ilp_all_nc_bd <- get_bd_feature_perf(ilp_all_nc)

# Top and bottom performers -----------------------------------------------

top_bd_features <- dplyr::bind_cols(rwr_all_nc_bd[["top10"]],
                                    zon_all_nc_bd[["top10"]],
                                    ilp_all_nc_bd[["top10"]])

low_bd_features <- dplyr::bind_cols(rwr_all_nc_bd[["low10"]],
                                    zon_all_nc_bd[["low10"]],
                                    ilp_all_nc_bd[["low10"]])

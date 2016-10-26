# NOTE: you will need the latest version for this to work
# zonator >= 0.5.3
if (!require("zonator")) {
  devtools::install_github("cbig/zonator", dependencies = TRUE)
}
library(raster)
library(zonator)

source("src/00_lib/utils.R")

# Number of features in groups ------------------------------------------------
# NOTE: the numbers are hard coded and are not updated if the data change.

NAMPHIBIANS <- 83
NBIRDS <- 404
NMAMMALS <- 164
NREPTILES <- 112

NESFEATURES <- 9
NBDFEATURES <- NAMPHIBIANS + NBIRDS + NMAMMALS + NREPTILES

# Generate variants for all taxa ----------------------------------------------

variants <- c("01_caz", "02_abf", "03_caz_wgt", "04_abf_wgt",
              "05_caz_es", "06_abf_es", "07_caz_bd", "08_abf_bd")

zsetup_root <- "analyses/zonation"

ppa_raster_file <- "../../../data/processed/eurostat/nuts_level2/NUTS_RG_01M_2013_level2.tif"
ppa_config_file <- "ppa_config.txt"

project_name <- "priocomp"

zonator::create_zproject(name = project_name, dir = zsetup_root, variants = variants,
                         dat_template_file = "analyses/zonation/template.dat",
                         spp_template_dir = "data/processed/features",
                         override_path = "../../../data/processed/features",
                         recursive = TRUE, overwrite = TRUE, debug = TRUE)
priocomp_zproject <- load_zproject(zsetup_root)

# Set run configuration parameters --------------------------------------------

## 01_caz ---------------------------------------------------------------------

variant1 <- get_variant(priocomp_zproject, 1)

# Set post-processing (LSM). First, let's create the file itself (zonator
# can't handle this yet). The file needs to be created only once per raxon
# since all the variants can use the same file.
ppa_file_name <- file.path(zsetup_root, project_name, ppa_config_file)
ppa_cmd_string <- paste(c("LSM", ppa_raster_file, 0, -1, 0), collapse = " ")
write(ppa_cmd_string, ppa_file_name)

# Need to define ppa_config.txt relative to the bat-file (same dir)-
variant1 <- set_dat_param(variant1, "post-processing list file",
                          ppa_config_file)

# Define groups. Features 1-12 are ecosystem services, Features 13-771 are
# species features.
groups(variant1) <- c(rep(1, NESFEATURES),
                      rep(2, NBDFEATURES))
groupnames(variant1) <- c("1" = "ecosystem_services", "2" = "species")
# Set groups use and groups file
variant1 <- set_dat_param(variant1, "use groups", 1)
# Note that groups file location is always relative to the bat file
groups_file <- file.path(variant1@name, paste0(variant1@name, "_groups.txt"))
variant1 <- set_dat_param(variant1, "groups file", groups_file)

# Save variant
save_zvariant(variant1, dir = file.path(zsetup_root, project_name),
              overwrite = TRUE, debug_msg = FALSE)

## 02_abf ---------------------------------------------------------------------

variant2 <- get_variant(priocomp_zproject, 2)
variant2 <- set_dat_param(variant2, "removal rule", 2)
variant2 <- set_dat_param(variant2, "post-processing list file",
                          ppa_config_file)

# Define groups. Features 1-12 are ecosystem services, Features 13-771 are
# species features.
groups(variant2) <- c(rep(1, NESFEATURES),
                      rep(2, NBDFEATURES))
groupnames(variant2) <- c("1" = "ecosystem_services", "2" = "species")
# Set groups use and groups file
variant2 <- set_dat_param(variant2, "use groups", 1)
# Note that groups file location is always relative to the bat file
groups_file <- file.path(variant2@name, paste0(variant2@name, "_groups.txt"))
variant2 <- set_dat_param(variant2, "groups file", groups_file)

save_zvariant(variant2, dir = file.path(zsetup_root, project_name),
              overwrite = TRUE, debug_msg = FALSE)

## 03_caz_wgt ----------------------------------------------------------------
variant3 <- get_variant(priocomp_zproject, 3)
variant3 <- set_dat_param(variant3, "removal rule", 1)
variant3 <- set_dat_param(variant3, "post-processing list file",
                          ppa_config_file)

# Define groups. Features 1-12 are ecosystem services, Features 13-771 are
# species features.
groups(variant3) <- c(rep(1, NESFEATURES),
                      rep(2, NBDFEATURES))
groupnames(variant3) <- c("1" = "ecosystem_services", "2" = "species")
# Set groups use and groups file
variant3 <- set_dat_param(variant3, "use groups", 1)
# Note that groups file location is always relative to the bat file
groups_file <- file.path(variant3@name, paste0(variant3@name, "_groups.txt"))
variant3 <- set_dat_param(variant3, "groups file", groups_file)

# Give weight of 759 / 11 to each ES. Altogether, there are 759 species
# features. Give weight 1 to each.
sppweights(variant3) <- c(rep(nfeatures(variant3) / NESFEATURES, NESFEATURES),
                          rep(1, nfeatures(variant3) - NESFEATURES))

save_zvariant(variant3, dir = file.path(zsetup_root, project_name),
              overwrite = TRUE, debug_msg = FALSE)

## 04_abf_wgt ----------------------------------------------------------------
variant4 <- get_variant(priocomp_zproject, 4)
variant4 <- set_dat_param(variant4, "removal rule", 2)
variant4 <- set_dat_param(variant4, "post-processing list file",
                          ppa_config_file)

# Define groups. Features 1-12 are ecosystem services, Features 13-771 are
# species features.
groups(variant4) <- c(rep(1, 12), rep(2, (nfeatures(variant4) - 12)))
groupnames(variant4) <- c("1" = "ecosystem_services", "2" = "species")
# Set groups use and groups file
variant4 <- set_dat_param(variant4, "use groups", 1)
# Note that groups file location is always relative to the bat file
groups_file <- file.path(variant4@name, paste0(variant4@name, "_groups.txt"))
variant4 <- set_dat_param(variant4, "groups file", groups_file)

# Give weight of 759 / 12 to each. Altogether, there are 759 spcies
# features. Give weight 1 to each.
sppweights(variant4) <- c(rep(nfeatures(variant4) / NESFEATURES, NESFEATURES),
                          rep(1, nfeatures(variant4) - NESFEATURES))

save_zvariant(variant4, dir = file.path(zsetup_root, project_name),
              overwrite = TRUE, debug_msg = FALSE)


# Just ecoystem services ----------------------------------------------------

## 05_caz_es ----------------------------------------------------------------

variant5 <- get_variant(priocomp_zproject, 5)

ppa_file_name <- file.path(zsetup_root, project_name, ppa_config_file)
ppa_cmd_string <- paste(c("LSM", ppa_raster_file, 0, -1, 0), collapse = " ")
write(ppa_cmd_string, ppa_file_name)
variant5 <- set_dat_param(variant5, "post-processing list file",
                          ppa_config_file)

# Select ONLY ES features
sppdata(variant5) <- sppdata(variant5)[1:NESFEATURES,]

# Save variant
save_zvariant(variant5, dir = file.path(zsetup_root, project_name),
              overwrite = TRUE, debug_msg = FALSE)

## 06_abf ----------------------------------------------------------------

variant6 <- get_variant(priocomp_zproject, 6)
variant6 <- set_dat_param(variant6, "removal rule", 2)
variant6 <- set_dat_param(variant6, "post-processing list file",
                          ppa_config_file)

ppa_file_name <- file.path(zsetup_root, project_name, ppa_config_file)
ppa_cmd_string <- paste(c("LSM", ppa_raster_file, 0, -1, 0), collapse = " ")
write(ppa_cmd_string, ppa_file_name)
variant6 <- set_dat_param(variant6, "post-processing list file",
                          ppa_config_file)

# Select ONLY ES features
sppdata(variant6) <- sppdata(variant6)[1:NESFEATURES,]

# Save variant
save_zvariant(variant6, dir = file.path(zsetup_root, project_name),
              overwrite = TRUE, debug_msg = FALSE)

# Just biodiversity features ------------------------------------------------

## 07_caz_bd ----------------------------------------------------------------

variant7 <- get_variant(priocomp_zproject, 7)

ppa_file_name <- file.path(zsetup_root, project_name, ppa_config_file)
ppa_cmd_string <- paste(c("LSM", ppa_raster_file, 0, -1, 0), collapse = " ")
write(ppa_cmd_string, ppa_file_name)
variant7 <- set_dat_param(variant7, "post-processing list file",
                          ppa_config_file)

# Select ONLY BD features
sppdata(variant7) <- sppdata(variant7)[(NESFEATURES + 1):nfeatures(variant7),]

# Define groups based on taxon. NOTE: for now, groups hard coded
groups(variant7) <- c(rep(1, NAMPHIBIANS),
                      rep(2, NBIRDS),
                      rep(3, NMAMMALS),
                      rep(4, NREPTILES))
groupnames(variant7) <- c("1" = "amphibians", "2" = "birds", "3" = "mammals",
                          "4" = "reptiles")
# Set groups use and groups file
variant7 <- set_dat_param(variant7, "use groups", 1)
# Note that groups file location is always relative to the bat file
groups_file <- file.path(variant7@name, paste0(variant7@name, "_groups.txt"))
variant7 <- set_dat_param(variant7, "groups file", groups_file)

# Save variant
save_zvariant(variant7, dir = file.path(zsetup_root, project_name),
              overwrite = TRUE, debug_msg = FALSE)

## 08_abf_bd ----------------------------------------------------------------

variant8 <- get_variant(priocomp_zproject, 8)
variant8 <- set_dat_param(variant8, "removal rule", 2)
variant8 <- set_dat_param(variant8, "post-processing list file",
                          ppa_config_file)

ppa_file_name <- file.path(zsetup_root, project_name, ppa_config_file)
ppa_cmd_string <- paste(c("LSM", ppa_raster_file, 0, -1, 0), collapse = " ")
write(ppa_cmd_string, ppa_file_name)
variant8 <- set_dat_param(variant8, "post-processing list file",
                          ppa_config_file)

# Select ONLY BD features
sppdata(variant8) <- sppdata(variant8)[(NESFEATURES + 1):nfeatures(variant8),]

# Define groups based on taxon. NOTE: for now, groups hard coded
groups(variant8) <- c(rep(1, NAMPHIBIANS),
                      rep(2, NBIRDS),
                      rep(3, NMAMMALS),
                      rep(4, NREPTILES))
groupnames(variant8) <- c("1" = "amphibians", "2" = "birds", "3" = "mammals",
                          "4" = "reptiles")
# Set groups use and groups file
variant8 <- set_dat_param(variant8, "use groups", 1)
# Note that groups file location is always relative to the bat file
groups_file <- file.path(variant8@name, paste0(variant8@name, "_groups.txt"))
variant8 <- set_dat_param(variant8, "groups file", groups_file)

# Save variant
save_zvariant(variant8, dir = file.path(zsetup_root, project_name),
              overwrite = TRUE, debug_msg = FALSE)

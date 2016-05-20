# NOTE: you will need the latest version for this to work
# zonator > 0.5.0
if (!require("zonator")) {
  devtools::install_github("cbig/zonator", dependencies = TRUE)
}
library(raster)
library(zonator)

source("src/utils.R")

# Generate variants for all taxa ------------------------------------------

# Define names for variants. "[ID]" is a placeholder for running id, "[TX]" is
# for taxon codes.
variants <- c("01_caz", "02_abf", "03_caz_wgt")

zsetup_root <- "analyses/zonation"

ppa_raster_file <- "../../../data/processed/nuts/NUTS_RG_01M_2013/level2/NUTS_RG_01M_2013_level2_subset.tif"
ppa_config_file <- "ppa_config.txt"

project_name <- "priocomp"

create_zproject(name = project_name, dir = zsetup_root, variants = variants,
                dat_template_file = "analyses/zonation/template.dat",
                spp_template_dir = "data/processed/features",
                override_path = "../../../data/processed/features",
                recursive = TRUE, overwrite = TRUE, debug = TRUE)
priocomp_zproject <- load_zproject(zsetup_root)

# Set run configuration parameters ----------------------------------------

## 01_caz

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

# Save variant
save_zvariant(variant1, dir = file.path(zsetup_root, project_name),
              overwrite = TRUE, debug_msg = FALSE)

## 02_abf

variant2 <- get_variant(priocomp_zproject, 2)
variant2 <- set_dat_param(variant2, "removal rule", 2)
variant2 <- set_dat_param(variant2, "post-processing list file",
                          ppa_config_file)
save_zvariant(variant2, dir = file.path(zsetup_root, project_name),
              overwrite = TRUE, debug_msg = FALSE)

library(R.cache)
library(zonator)

# Common data -------------------------------------------------------------

# Red-list status used as a basis for grouping
rl_groups <- list(LC = 1, NT = 2, VU = 3, EN = 4, CR = 5, EW = 6, EX = 7, DD = 8)

# Common functions ---------------------------------------------------------

#' Create a set of variant on file system.
#'
#' @param zsetup_root Character path to the parent dir of variants.
#' @param variant_templates Character vector containing the templates for variants
#'   to be created.
#' @param spp_data List holding auxiliary data for taxa.
#' @param data_dir Character path to data directory.
#' @param prefix_spp_paths Prefix path inserted into spp files.
#' @param dat_template_file Character path to template dat file used.
#'
#' @return Zproject object
#'
initiate_zproject <- function(zsetup_root = "zsetup", variant_templates,
                              spp_data, data_dir = "../Data.150928",
                              prefix_spp_paths = "../..",
                              dat_template_file = "templates/template.dat") {
  variant_id <- 1

  for (taxon in names(spp_data)) {
    variant_names <- variant_templates
    sub_project_dir <- file.path(zsetup_root, taxon)
    message("Creating sub-project in ", sub_project_dir)

    # Generate ids dynamically
    ids <- sprintf("%02d_",
                   variant_id:(variant_id + length(variant_names) - 1))

    # Generate actual variant names
    variant_names <- paste0(ids,
                            gsub(pattern = "\\[TX\\]",
                                 replacement = spp_data[[taxon]]$code,
                                 variant_names))

    # NOTE: set the data path manually
    create_zproject(name = sub_project_dir, dir = ".",
                    variants = variant_names,
                    dat_template_file = dat_template_file,
                    spp_template_dir = file.path(data_dir, spp_data[[taxon]]$sheet),
                    override_path = file.path(prefix_spp_paths, data_dir,
                                              spp_data[[taxon]]$sheet),
                    overwrite = TRUE, debug = FALSE)

    # Update variant ID for the next taxon
    variant_id <- variant_id + length(variant_names)
    message(" Created variants: \n  ", paste(variant_names, collapse = "\n  "))
  }

  # Finally, create variants with all species together
  sub_project_dir <- file.path(zsetup_root, "taxa_all")
  message("Creating sub-project in ", sub_project_dir)
  # No need for the taxa names anymore
  variant_names <- gsub("\\[TX\\]_", "", variant_templates)
  ids <- sprintf("%02d_", variant_id:(variant_id + length(variant_names) - 1))
  variant_names <- paste0(ids, variant_names)
  # Use "recursive = TRUE" with specific spp name template to get all the
  # rasters from different input folders
  create_zproject(name = sub_project_dir,
                  dir = ".",
                  variants = variant_names,
                  dat_template_file = dat_template_file,
                  spp_template_dir = data_dir,
                  override_path = file.path(prefix_spp_paths, data_dir),
                  overwrite = TRUE,
                  debug = FALSE,
                  recursive = TRUE,
                  # The pattern will only take file names with upper case first
                  # character
                  spp_file_pattern = "[A-Z][a-z]+_.+\\.(tif|img)$")
  message(" Created variants: \n  ", paste(variant_names, collapse = "\n  "))

  return(load_zproject(zsetup_root))
}

#' Load Zonation results.
#'
#' @param zsetup_root Character path to the parent dir of variants.
#' @param cache
#' @param ... other arguments passed on to \code{load_zproject()}.
#'
#' @return Zproject object
#'
.load_zproject <- function(zsetup_root = "zsetup", cache = TRUE, ...) {


  if (cache) {
    # Create the key
    key <- list(zsetup_root)
    zproject <- loadCache(key)
    if (!is.null(zproject)) {
      message("Loaded cached data")
      return(zproject)
    } else {
      warning("Couldn't find cached data, reloading...")
      zproject <- zonator::load_zproject(zsetup_root, ...)
      saveCache(zproject, key = key)
    }
  } else {
    zproject <- zonator::load_zproject(zsetup_root, ...)
  }
  return(zproject)
}

# Lookup a numeric value for red-list group for a species.
# @param table Dataframe table containing the data for a given taxon.
# @param species String character defining species.
lookup_rl_group <-
  Vectorize(function(table, species) {

    if (!"st.species" %in% names(table) | !"category" %in% names(table)) {
      stop("Table must have columns 'st.species' and 'category'")
    }
    if (!species %in% table$st.species) {
      warning(paste0("Species provided (", species, ") not found in the table."))
      return(NA)
    }
    if (!all(table$category %in% names(rl_groups))) {
      stop("All category values in table must be in: ", paste(names(rl_groups),
                                                              collapse = ", "))
    }
    rl_cat <- table[which(table$st.species == species), ]$category
    rl_cat <- unlist(rl_groups[rl_cat])
    # Remove vector names
    rl_cat <- as.vector(rl_cat)
    return(rl_cat)
  },
  c("species"), USE.NAMES = TRUE)

# Lookup a numeric value for weight for a species.
# @param table Dataframe table containing the data for a given taxon.
# @param species String character defining species.
lookup_weight <-
  Vectorize(function(table, species) {

    if (!"st.species" %in% names(table) | !"WC" %in% names(table)) {
      stop("Table must have columns 'st.species' and 'WC'")
    }
    if (!species %in% table$st.species) {
      warning(paste0("Species provided (", species, ") not found in the table."))
      return(NA)
    }
    if (!species %in% table$st.species) {
      warning(paste0("Species provided (", species, ") not found in the table."))
      return(NA)
    }

    wc <- table[which(table$st.species == species), ]$WC
    if (length(wc) > 1) {
      wc <- wc[1]
      warning("Multiple weights found for a single species, using the first: ",
              wc)
    }
    return(wc)
  },
  c("species"), USE.NAMES = TRUE)

#' Write Zonation groups file based on a Zvariant object.
#'
#' Function accesses \code{groups} slot of a \code{Zvariant} object and writes
#' the content into a file (if groups are assigned). Accessing Zvariant object's
#' slost directly not a generally good idea. Done here until zonator actually
#' supports writing out files properly (does not yet in 0.4.1).
#'
#' @note This functionality should eventually be incorporated into
#'  \code{zonator}.
#'
#' @param x Zvariant object (with groups).
#' @param filename String character file path to be written.
#' @parma overwrite Logical indicating if an existing file should be overwritten.
#'
#' @return Invisible NULL.
#'
#' @author Joona Lehtomaki \email{joona.lehtomaki@@gmail.com}
#'
write_groups <- function(x, filename, overwrite = FALSE) {
  if (!class(x) == "Zvariant") {
    stop("Object provided must be of Zvariant class")
  }
  if (all(dim(x@groups) == c(0, 0))) {
    stop("Zvariant must have groups.")
  }
  groups <- x@groups
  # Strip out "name" (last) column
  groups <- groups[, -ncol(groups)]
  if (!file.exists(filename) & !overwrite) {
    stop("File exists and overwrite is off.")
  }
  write.table(groups, file = filename, sep = "\t", row.names = FALSE,
              col.names = FALSE)
  message("Wrote groups file ", filename)
  return(invisible(NULL))
}

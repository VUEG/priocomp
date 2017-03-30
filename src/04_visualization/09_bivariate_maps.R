library(viridis)
library(zonator)

source("src/00_lib/bivariate_map.R")

# Setup -------------------------------------------------------------------

zproject <- zonator::load_zproject('analyses/zonation/priocomp/')

# 04_abf_all_wgt
v04 <- zonator::get_variant(zproject, 4)
v04_rank_raster <- zonator::rank_raster(v04)

# 08_abf_es
v08 <- zonator::get_variant(zproject, 8)
v08_rank_raster <- zonator::rank_raster(v08)

# 12_abf_bd
v12 <- zonator::get_variant(zproject, 12)
v12_rank_raster <- zonator::rank_raster(v12)


# Plot maps ---------------------------------------------------------------

n_quantiles <- 10

# Create a color matrix
col_mtrx <- colmat(nquantiles = n_quantiles,
                   upperleft = "#B2AF09", upperright = "#FF000B",
                   bottomleft = "grey", bottomright = "#00A3FF",
                   xlab = "ES", ylab = "BD")

bivmap <- bivariate.map(v08_rank_raster, v12_rank_raster,
                        colormatrix = col_mtrx, nquantiles = n_quantiles)

# Plot the bivariate map:
p1 <- plot(bivmap, frame.plot = FALSE, axes = FALSE, box = FALSE, add = FALSE,
           legend = FALSE, col = as.vector(col_mtrx))

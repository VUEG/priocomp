library(raster)
library(rasterVis)
library(imagemetrics)
library(viridis)

path1 <- "data/processed/features/provide/carbon_sequestration/carbon_sequestration.tif"
path2 <- "data/processed/features/provide/cultural_landscape_index_agro/cultural_landscape_index_agro.tif"
path3 <- "data/processed/features/provide/cultural_landscape_index_forest/cultural_landscape_index_forest.tif"

#############################################################################
#### STEP 1 : Load or calculate all image layers we need to work with.
#############################################################################

# Load 3 rasters from the target image,
# one for red, one for green and one for the blue channel
R1 <- raster(path1, band = 1)
R2 <- raster(path2, band = 1)
R3 <- raster(path3, band = 1)

Rstack <- stack(R1, R2, R3)

levelplot(Rstack, margin = FALSE, col.regions = viridis)

#############################################################################
#### STEP 2 : Calculate neighbor matrices and histograms on each layer.
####          We will then feed these into the actual analysis functions
#############################################################################

# On the four channels, find either right (1), diagonal (2) or below (3)
# neighbors for histogram calculations
v_R1_1 <- getImagePixels(R1, side = 1)
v_R2_1 <- getImagePixels(R2, side = 1)
v_R3_1 <- getImagePixels(R3, side = 1)

v_R1_2 <- getImagePixels(R1, side = 2)
v_R2_2 <- getImagePixels(R2, side = 2)
v_R3_2 <- getImagePixels(R3, side = 2)

v_R1_3 <- getImagePixels(R1, side = 3)
v_R2_3 <- getImagePixels(R2, side = 3)
v_R3_3 <- getImagePixels(R3, side = 3)

# Calculate histograms from neighbor vectors
nbins <- 15

prob_R1_1 <- calculateHisto(
  reference_vector = v_R1_1$reference_vector,
  neighbour_vector = v_R1_1$neighbour_vector,
  nbins = nbins
)
prob_R2_1 <- calculateHisto(
  reference_vector = v_R2_1$reference_vector,
  neighbour_vector = v_R2_1$neighbour_vector,
  nbins = nbins
)
prob_R3_1 <- calculateHisto(
  reference_vector = v_R3_1$reference_vector,
  neighbour_vector = v_R3_1$neighbour_vector,
  nbins = nbins
)

prob_R1_2 <- calculateHisto(
  reference_vector = v_R1_2$reference_vector,
  neighbour_vector = v_R1_2$neighbour_vector,
  nbins = nbins
)
prob_R2_2 <- calculateHisto(
  reference_vector = v_R2_2$reference_vector,
  neighbour_vector = v_R2_2$neighbour_vector,
  nbins = nbins
)
prob_R3_2 <- calculateHisto(
  reference_vector = v_R3_2$reference_vector,
  neighbour_vector = v_R3_2$neighbour_vector,
  nbins = nbins
)

prob_R1_3 <- calculateHisto(
  reference_vector = v_R1_3$reference_vector,
  neighbour_vector = v_R1_3$neighbour_vector,
  nbins = nbins
)
prob_R2_3 <- calculateHisto(
  reference_vector = v_R2_3$reference_vector,
  neighbour_vector = v_R2_3$neighbour_vector,
  nbins = nbins
)
prob_R3_3 <- calculateHisto(
  reference_vector = v_R3_3$reference_vector,
  neighbour_vector = v_R3_3$neighbour_vector,
  nbins = nbins
)

#############################################################################
#### STEP 3 : Calculate the metrics.
#############################################################################


# Calculate mean information gain
meanInformationGain(prob_R1_2)
meanInformationGain(prob_R2_2)
meanInformationGain(prob_R3_2)

# Calculate anisotropy
meanInformationGain(prob_R1_1) / meanInformationGain(prob_R1_3)
meanInformationGain(prob_R2_1) / meanInformationGain(prob_R2_3)
meanInformationGain(prob_R2_1) / meanInformationGain(prob_R2_3)

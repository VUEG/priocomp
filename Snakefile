import dotenv
import logging
import numpy as np
import os
import rasterio
import requests
import yaml
from snakemake import logger
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from importlib.machinery import SourceFileLoader

rescale = SourceFileLoader("data_processing.rescale", "src/data_processing/rescale.py").load_module()

## GLOBALS ---------------------------------------------------------------------

# dotenv project variables
dotenv_path = '.env'
dotenv.load_dotenv(dotenv_path)

# Set up a remote provider
HTTP = HTTPRemoteProvider()

DATADRYAD_DATASETS = ['woodprod_average']

PROVIDE_DATASETS = ['agrobiodiversity_species_richness',
                    'carbon_sequestration',
                    'cultural_landscape_index_agro',
                    'cultural_landscape_index_forest',
                    'erosion_prevention',
                    'flood_regulation_supply',
                    'floodregulation',
                    'megafauna',
                    'nature_tourism',
                    'pollination_flows',
                    'pollination_visitprob',
                    'species_richness_farmland_birds_original1',
                    'species_richness_vascular_plants']

# PROJECT RULES ----------------------------------------------------------------

rule all:
    input:
        expand(["data/external/{dataset}/datapackage.json", "data/external/{dataset}/{dataset}.tif"], dataset=PROVIDE_DATASETS)

## Get data --------------------------------------------------------------------

rule get_datadryad_data:
    input:
        HTTP.remote(expand(["beehub.nl/environmental-geography-group/datadryad/forest_production_europe/datapackage.json",
                            "beehub.nl/environmental-geography-group/datadryad/forest_production_europe/README.txt",
                            "beehub.nl/environmental-geography-group/datadryad/forest_production_europe/{dataset}.tif"], dataset=DATADRYAD_DATASETS),
                    username=os.environ.get("BEEHUB_USERNAME"), password=os.environ.get("BEEHUB_PASSWORD"), keep_local=False)
    output:
        expand(["data/external/datadryad/forest_production_europe/datapackage.json",
                "data/external/datadryad/forest_production_europe/README.txt",
                "data/external/datadryad/forest_production_europe/{dataset}.tif"], dataset=DATADRYAD_DATASETS)
    log:
        "logs/data_datadryad.log"
    run:
        # Configure logger
        fileHandler = logging.FileHandler(log[0])
        fileHandler.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
        logger.logger.addHandler(fileHandler)

        for i in range(0, len(input)):
            shell("mv {0} {1}".format(input[i], output[i]))
            logger.info("Downloaded {0} to {1}".format(input[i], output[i]))

rule get_provide_data:
    input:
        HTTP.remote(expand(["beehub.nl/environmental-geography-group/provide/{dataset}/datapackage.json",
                            "beehub.nl/environmental-geography-group/provide/{dataset}/{dataset}.tif"], dataset=PROVIDE_DATASETS),
                    username=os.environ.get("BEEHUB_USERNAME"), password=os.environ.get("BEEHUB_PASSWORD"), keep_local=False)
    output:
        expand(["data/external/provide/{dataset}/datapackage.json", "data/external/provide/{dataset}/{dataset}.tif"], dataset=PROVIDE_DATASETS)
    log:
        "logs/data_provide.log"
    run:
        # Configure logger
        fileHandler = logging.FileHandler(log[0])
        fileHandler.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
        logger.logger.addHandler(fileHandler)

        for i in range(0, len(input)):
            shell("mv {0} {1}".format(input[i], output[i]))
            logger.info("Downloaded {0} to {1}".format(input[i], output[i]))

## Data pre-processing ---------------------------------------------------------

rule build_data_coverage:
    input:
        expand("data/interim/rescaled/provide/{dataset}/{dataset}.tif", dataset=PROVIDE_DATASETS) + \
        expand("data/interim/rescaled/datadryad/forest_production_europe/{dataset}.tif", dataset=DATADRYAD_DATASETS)
    output:
        "data/processed/common/data_coverage.tif"
    log:
        "logs/build_data_coverage.log"
    message:
        "Building data coverage raster..."
    run:
        # In the first round, get the shape of the rasters and construct
        # a np.array to hold the binary data maskas. NOTE: the assumption
        # here is the all data have the same dimension
        dims = None
        masks = None
        profile = None
        for i, s_raster in enumerate(input):
            with rasterio.open(s_raster) as in_src:
                if i == 0:
                    dims = in_src.shape
                    # Also store template profile
                    profile = in_src.profile
                    masks = np.zeros((len(input), dims[0], dims[1]), dtype=np.uint16)
                else:
                    if in_src.shape != dims:
                        logger.warning(" WARNING: Dimensions for {} don't match, skipping".format(s_raster))
                        continue
                # Place a binary version of the current raster mask into the
                # container
                logger.info(" [{0}/{1}] Reading mask from {2}".format(i+1, len(input), s_raster))
                mask = in_src.read_masks(1)
                # We're using the GDAL type mask where any values above from
                # 0 are actual data values. Convert (in-place) to binary.
                np.place(mask, mask > 0, 1)
                masks[i] = mask

        logger.info(" Summing mask information and saving to {}".format(output))
        data_coverage = masks.sum(axis=0)
        # Write the product.
        profile.update(dtype=rasterio.uint16, compress="DEFLATE")

        with rasterio.open(output[0], 'w', **profile) as dst:
            dst.write(data_coverage.astype(rasterio.uint16), 1)

rule harmonize_data:
    input:
        expand("data/external/provide/{dataset}/{dataset}.tif", dataset=PROVIDE_DATASETS) + \
        expand("data/external/datadryad/forest_production_europe/{dataset}.tif", dataset=DATADRYAD_DATASETS)
    output:
        expand("data/interim/harmonized/provide/{dataset}/{dataset}.tif", dataset=PROVIDE_DATASETS) + \
        expand("data/interim/harmonized/datadryad/forest_production_europe/{dataset}.tif", dataset=DATADRYAD_DATASETS)
    params:
        # Snap raster
        like_raster = "data/external/provide/carbon_sequestration/carbon_sequestration.tif",
        # Target CRS
        dst_src = 3035
    message:
        "Harmonizing datasets..."
    run:
        for i, s_raster in enumerate(input):
            # No need to process the snap raster, just copy it
            if s_raster == params['like_raster']:
                logger.info(" [{0}/{1}] Copying dataset {2}".format(i+1, len(input), s_raster))
                shell("cp {0} {1}".format(input[i], output[i]))
                continue
            # No need to process the snap raster
            logger.info(" [{0}/{1}] Warping dataset {2}".format(i+1, len(input), s_raster))
            shell("rio warp " + input[i] + " --like " + params['like_raster'] + \
                  " " + output[i] + " --dst-crs " + str(params['dst_src']) + \
                  " --co 'COMPRESS=DEFLATE' --threads {threads}")

rule rescale_data:
    input:
        expand("data/interim/harmonized/provide/{dataset}/{dataset}.tif", dataset=PROVIDE_DATASETS) + \
        expand("data/interim/harmonized/datadryad/forest_production_europe/{dataset}.tif", dataset=DATADRYAD_DATASETS)
    output:
        expand("data/interim/rescaled/provide/{dataset}/{dataset}.tif", dataset=PROVIDE_DATASETS) + \
        expand("data/interim/rescaled/datadryad/forest_production_europe/{dataset}.tif", dataset=DATADRYAD_DATASETS)
    message:
        "Rescaling data..."
    run:
        for i, s_raster in enumerate(input):
            # No need to process the snap raster
            logger.info(" [{0}/{1}] Rescaling dataset {2}".format(i+1, len(input), s_raster))
            rescale.rescale_raster(input[i], output[i], method="normalize", verbose=True)

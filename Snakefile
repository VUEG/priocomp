import dotenv
import logging
import os
import requests
import yaml
from snakemake import logger
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider


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

rule harmonize_data:
    input:
        expand("data/external/provide/{dataset}/{dataset}.tif", dataset=PROVIDE_DATASETS) + \
        expand("data/external/datadryad/forest_production_europe/{dataset}.tif", dataset=DATADRYAD_DATASETS)
    output:
        expand("data/interim/provide/{dataset}/{dataset}.tif", dataset=PROVIDE_DATASETS) + \
        expand("data/interim/datadryad/forest_production_europe/{dataset}.tif", dataset=DATADRYAD_DATASETS)
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

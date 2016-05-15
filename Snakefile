import dotenv
import fiona
import gdal
import logging
import numpy as np
import os
import rasterio
import requests
import sys
import yaml
from snakemake import logger
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from importlib.machinery import SourceFileLoader

utils = SourceFileLoader("src.utils", "src/utils.py").load_module()
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

SHP_COMPONENTS = ["dbf", "prj", "shx", "shp"]

# Analysis extent
PROJECT_EXTENT = {"bottom": 1000000.0, "left": 2000000.0, "right": 6526000.0,
                  "top": 5410000.0}
# EPSG for project
PROJECT_CRS = 3035

# Offset the bounds given in extent_yml. Values are in order
# (left, bottom, right, top) and interpreted in the CRS units. values
# are added to bounds given by extent.yml
OFFSET = (100000, 100000, 0, 0)
# Which eurostat countries are included in the processed output? The
# following countries have been removed:
# "CY", "CH", "IS", "HR", "NO", "ME", "MT", "MK", "TR"
PROJECT_COUNTRIES = ["AT", "BE", "BG", "CZ", "DE", "DK", "ES", "EL", "EE",
                     "FR", "FI", "IT", "HU", "IE", "NL", "LU", "LI", "LT",
                     "LV", "PL", "SE", "RO", "PT", "SK", "SI", "UK"]

# PROJECT RULES ----------------------------------------------------------------

rule all:
    input:
        expand("data/processed/provide/{dataset}/{dataset}.tif", dataset=PROVIDE_DATASETS) + \
        expand("data/processed/datadryad/forest_production_europe/{dataset}.tif", dataset=DATADRYAD_DATASETS)

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

rule get_nuts_data:
    input:
        HTTP.remote(expand(["beehub.nl/environmental-geography-group/nuts/NUTS_RG_01M_2013/level0/NUTS_RG_01M_2013_level0.{ext}",
                            "beehub.nl/environmental-geography-group/nuts/NUTS_RG_01M_2013/level2/NUTS_RG_01M_2013_level2.{ext}"],
                            ext=SHP_COMPONENTS),
                    username=os.environ.get("BEEHUB_USERNAME"), password=os.environ.get("BEEHUB_PASSWORD"), keep_local=False)
    output:
        expand(["data/external/nuts/NUTS_RG_01M_2013/level0/NUTS_RG_01M_2013_level0.{ext}",
                "data/external/nuts/NUTS_RG_01M_2013/level2/NUTS_RG_01M_2013_level2.{ext}"],
                ext=SHP_COMPONENTS)
    log:
        "logs/data_nuts.log"
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
        expand("data/processed/provide/{dataset}/{dataset}.tif", dataset=PROVIDE_DATASETS) + \
        expand("data/processed/datadryad/forest_production_europe/{dataset}.tif", dataset=DATADRYAD_DATASETS)
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
        external=expand("data/external/provide/{dataset}/{dataset}.tif", dataset=PROVIDE_DATASETS) + \
                 expand("data/external/datadryad/forest_production_europe/{dataset}.tif", dataset=DATADRYAD_DATASETS),
        like_raster="data/external/datadryad/forest_production_europe/woodprod_average.tif",
        clip_shp=expand("data/processed/nuts/NUTS_RG_01M_2013/level0/NUTS_RG_01M_2013_level0_subset.{ext}",
                         ext=SHP_COMPONENTS)
    output:
        warped=temp(expand("data/interim/warped/provide/{dataset}/{dataset}.tif", dataset=PROVIDE_DATASETS) + \
                     expand("data/interim/warped/datadryad/forest_production_europe/{dataset}.tif", dataset=DATADRYAD_DATASETS)),
        harmonized=expand("data/interim/harmonized/provide/{dataset}/{dataset}.tif", dataset=PROVIDE_DATASETS) + \
                   expand("data/interim/harmonized/datadryad/forest_production_europe/{dataset}.tif", dataset=DATADRYAD_DATASETS)
    message:
        "Harmonizing datasets..."
    run:
        ndatasets = len(input.external)
        for i, s_raster in enumerate(input.external):
            ## WARP
            # Target raster
            warped_raster = s_raster.replace("external", "interim/warped")
            # No need to process the snap raster, just copy it
            if s_raster == input.like_raster:
                logger.info(" [{0}/{1} step 1] Copying dataset {2}".format(i+1, ndatasets, s_raster))
                logger.debug(" [{0}/{1} step 1] Target dataset {2}".format(i+1, ndatasets, warped_raster))
                shell("cp {s_raster} {warped_raster}")
            else:
                logger.info(" [{0}/{1} step 1] Warping dataset {2}".format(i+1, ndatasets, s_raster))
                logger.debug(" [{0}/{1} step 1] Target dataset {2}".format(i+1, ndatasets, warped_raster))
                shell("rio warp " + s_raster + " --like " + input.like_raster + \
                      " " + warped_raster + " --dst-crs " + str(PROJECT_CRS) + \
                      " --co 'COMPRESS=DEFLATE' --threads {threads}")

            ## CLIP
            harmonized_raster = warped_raster.replace("interim/warped", "interim/harmonized")
            clip_shp = utils.pick_from_list(input.clip_shp, ".shp")
            logger.info(" [{0}/{1} step 2] Clipping dataset {2}".format(i+1, ndatasets, warped_raster))
            logger.debug(" [{0}/{1} step 2] Target dataset {2}".format(i+1, ndatasets, harmonized_raster))
            shell("gdalwarp -cutline {clip_shp} {warped_raster} {harmonized_raster} -co COMPRESS=DEFLATE")

rule preprocess_nuts_level0_data:
    input:
        shp=expand("data/external/nuts/NUTS_RG_01M_2013/level0/NUTS_RG_01M_2013_level0.{ext}",
                   ext=SHP_COMPONENTS)
    output:
        reprojected=temp(expand("data/interim/nuts/NUTS_RG_01M_2013/level0/NUTS_RG_01M_2013_level0.{ext}",
                                ext=SHP_COMPONENTS)),
        processed=expand("data/processed/nuts/NUTS_RG_01M_2013/level0/NUTS_RG_01M_2013_level0_subset.{ext}",
                         ext=SHP_COMPONENTS)
    message:
        "Pre-processing NUTS level 0 data..."
    run:
        # Read in the bounds as used in harmonize_data rule
        bleft = PROJECT_EXTENT["left"] + OFFSET[0]
        bbottom = PROJECT_EXTENT["bottom"] + OFFSET[1]
        bright = PROJECT_EXTENT["right"] + OFFSET[2]
        btop = PROJECT_EXTENT["top"] + OFFSET[3]
        bounds = "{0} {1} {2} {3}".format(bleft, bbottom, bright, btop)
        # Reproject to EPSG:3035 from EPSG:4258
        input_shp = utils.pick_from_list(input.shp, ".shp")
        reprojected_shp = utils.pick_from_list(output.reprojected, ".shp")
        shell('ogr2ogr {reprojected_shp} -t_srs "EPSG:{PROJECT_CRS}" {input_shp}')
        logger.debug("Reprojected NUTS data from EPSG:4258 to EPSG:3035")
        # Clip shapefile using ogr2ogr, syntax:
        # ogr2ogr output.shp input.shp -clipsrc <left> <bottom> <right> <top>
        processed_shp = utils.pick_from_list(output.processed, ".shp")
        # Do 2 things at the same time:
        #  1. Select a subset of counties (defined by params.countries)
        #  2. Clip output to an extent (given by bounds)
        # Build the -where clause for ogr2ogr
        where_clause = "NUTS_ID IN ({})".format(", ".join(["'" + item + "'" for item in PROJECT_COUNTRIES]))
        shell('ogr2ogr -where "{where_clause}" {processed_shp} {reprojected_shp} -clipsrc {bounds}')
        logger.debug("Clipped NUTS data to analysis bounds: {}".format(bounds))
        logger.debug("Selected only a subset of eurostat countries")

rule preprocess_nuts_level2_data:
    input:
        shp=expand("data/external/nuts/NUTS_RG_01M_2013/level2/NUTS_RG_01M_2013_level2.{ext}",
                   ext=SHP_COMPONENTS)
    output:
        reprojected=temp(expand("data/interim/nuts/NUTS_RG_01M_2013/level2/NUTS_RG_01M_2013_level2.{ext}",
                                ext=SHP_COMPONENTS)),
        enhanced=temp(expand("data/interim/nuts/NUTS_RG_01M_2013/level2/NUTS_RG_01M_2013_level2_enhanced.{ext}",
                                ext=SHP_COMPONENTS)),
        processed=expand("data/processed/nuts/NUTS_RG_01M_2013/level2/NUTS_RG_01M_2013_level2_subset.{ext}",
                         ext=SHP_COMPONENTS)
    message:
        "Pre-processing NUTS level 2 data..."
    run:
        # Read in the bounds as used in harmonize_data rule
        bleft = PROJECT_EXTENT["left"] + OFFSET[0]
        bbottom = PROJECT_EXTENT["bottom"] + OFFSET[1]
        bright = PROJECT_EXTENT["right"] + OFFSET[2]
        btop = PROJECT_EXTENT["top"] + OFFSET[3]
        bounds = "{0} {1} {2} {3}".format(bleft, bbottom, bright, btop)
        # Reproject to EPSG:3035 from EPSG:4258
        input_shp = utils.pick_from_list(input.shp, ".shp")
        reprojected_shp = utils.pick_from_list(output.reprojected, ".shp")
        shell('ogr2ogr {reprojected_shp} -t_srs "EPSG:{PROJECT_CRS}" {input_shp}')
        logger.debug("Reprojected NUTS level 2 data from EPSG:4258 to EPSG:3035")

        # The Pre-processing steps need to be done:
        #  1. Tease apart country code from field NUTS_ID
        #  2. Create a running ID field that can be used as value in the
        #     rasterized version
        enhanced_shp = utils.pick_from_list(output.enhanced, ".shp")
        with fiona.drivers():
            with fiona.open(reprojected_shp) as source:
                meta = source.meta
                meta['schema']['geometry'] = 'Polygon'
                # Insert new fields
                meta['schema']['properties']['ID'] = 'int'
                meta['schema']['properties']['country'] = 'str'

                ID = 1
                with fiona.open(enhanced_shp, 'w', **meta) as sink:
                    # Loop over features
                    for f in source:
                        # Check the country code part of NUTS_ID (2 first
                        # charatcters). NOTE: we're effectively doing filtering
                        # here.
                        country_code = f['properties']['NUTS_ID'][0:2]
                        if country_code in PROJECT_COUNTRIES:
                            f['properties']['ID'] = ID
                            ID += 1
                            f['properties']['country'] = country_code
                            # Write the record out.
                            sink.write(f)

        # Clip shapefile using ogr2ogr, syntax:
        # ogr2ogr output.shp input.shp -clipsrc <left> <bottom> <right> <top>
        processed_shp = utils.pick_from_list(output.processed, ".shp")
        # Clip output to an extent (given by bounds)
        shell('ogr2ogr {processed_shp} {enhanced_shp} -clipsrc {bounds}')
        logger.debug("Clipped NUTS level 2 data to analysis bounds: {}".format(bounds))
        logger.debug("Selected only a subset of eurostat countries")

rule rescale_data:
    input:
        expand("data/interim/harmonized/provide/{dataset}/{dataset}.tif", dataset=PROVIDE_DATASETS) + \
        expand("data/interim/harmonized/datadryad/forest_production_europe/{dataset}.tif", dataset=DATADRYAD_DATASETS)
    output:
        expand("data/processed/provide/{dataset}/{dataset}.tif", dataset=PROVIDE_DATASETS) + \
        expand("data/processed/datadryad/forest_production_europe/{dataset}.tif", dataset=DATADRYAD_DATASETS)
    message:
        "Rescaling data..."
    run:
        for i, s_raster in enumerate(input):
            # No need to process the snap raster
            logger.info(" [{0}/{1}] Rescaling dataset {2}".format(i+1, len(input), s_raster))
            # NOTE: looping over input and output only works if they have
            # exactly the same definition. Otherwise order may vary.
            rescale.rescale_raster(input[i], output[i], method="normalize", verbose=True)

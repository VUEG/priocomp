import fiona
import geopandas as gpd
import gdal
import numpy as np
import numpy.ma as ma
import os
import pandas as pd
import pdb
import pypandoc
import rasterio
import rasterio.tools.mask
import rasterstats
import sys
from importlib.machinery import SourceFileLoader

## LOAD MODULES --------------------------------------------------------------
utils = SourceFileLoader("lib.utils", "src/00_lib/utils.py").load_module()
spatutils = SourceFileLoader("lib.spatutils", "src/00_lib/spatutils.py").load_module()
gurobi = SourceFileLoader("analysis.gurobi", "src/02_analysis/gurobi.py").load_module()
rwr = SourceFileLoader("analysis.rwr", "src/02_analysis/rwr.py").load_module()
similarity = SourceFileLoader("results.similarity", "src/03_post_processing/similarity.py").load_module()
coverage = SourceFileLoader("results.coverage", "src/03_post_processing/data_coverage.py").load_module()

## GLOBALS --------------------------------------------------------------------

# Analysis extent
PROJECT_EXTENT = {"bottom": 1000000.0, "left": 2000000.0, "right": 6526000.0,
                  "top": 5410000.0}
# EPSG for project
PROJECT_CRS = 3035

# Project resolution in units of PROJECT_CRS
PROJECT_RES = 100

# Offset the bounds given in extent_yml. Values are in order
# (left, bottom, right, top) and interpreted in the CRS units. values
# are added to bounds given by PROJECT_EXTENT
OFFSET = (100000, 100000, 0, 0)
# Which eurostat countries are included in the processed output? The
# following countries have been removed:
# "CY", "CH", "IS", "HR", "NO", "ME", "MT", "MK", "TR", "LI"
PROJECT_COUNTRIES = ["AT", "BE", "BG", "CZ", "DE", "DK", "ES", "EL", "EE",
                     "FR", "FI", "IT", "HU", "IE", "NL", "LU", "LT",
                     "LV", "PL", "SE", "RO", "PT", "SK", "SI", "UK"]

# Load the content of the data manifest file into a DataManager object.
dm = utils.DataManager("data/data_manifest.yml")

external_data = "data/external"
# Final feature data
feature_data = "data/processed/features"
beehub_url = "https://beehub.nl/environmental-geography-group"

# Define source and desination datasets. NOTE: data/Snakefile must be run
# before this Snakefile will work
DATADRYAD_SRC_DATASETS = [url.replace(beehub_url, external_data) for url in dm.get_resources(provider="datadryad", full_path=True)]

EEA_SRC_DATASETS = [url.replace(beehub_url, external_data) for url in dm.get_resources(provider="eea", full_path=True)]

# Get specific NUTS collections from Eurostat
NUTS_LEVEL0_DATA = [url.replace(beehub_url, external_data) for url in dm.get_resources(collection="nuts_level0", full_path=True)]
NUTS_LEVEL2_DATA = [url.replace(beehub_url, external_data) for url in dm.get_resources(collection="nuts_level2", full_path=True)]

PROVIDE_SRC_DATASETS = [url.replace(beehub_url, external_data) for url in dm.get_resources(provider="provide", full_path=True)]

# UDR collection "european_tetrapods" is already correctly formatted, place
# it directly to "processed/features"
UDR_SRC_DATASETS = [url.replace(beehub_url, external_data) for url in dm.get_resources(provider="udr", full_path=True)]

# Construct lists for processed BD and ES datasets
BD_DST_DATASETS = [url.replace(beehub_url, feature_data) for url in dm.get_resources(category="biodiversity", full_path=True)]
ES_DST_DATASETS = [url.replace(beehub_url, feature_data) for url in dm.get_resources(category="ecosystemservices", full_path=True)]

ALL_SRC_DATASETS = DATADRYAD_SRC_DATASETS + PROVIDE_SRC_DATASETS + UDR_SRC_DATASETS

# Let's build weight vectors needed for the weighted analysis variants

# Set the weight of each feature in category "biodiversity" to 1.0 and the
# weight of each feature in "ecosystemservices" to sum(bd_wights) / n_es.
# NOTE: the above weighting scheme is needed to avoid small weights values that
# e.g. 1.0 / n_bd would produce. For some reason the ILP implementation
# doesn't like small weight values
# (see: https://github.com/VUEG/priocomp/issues/8)
# NOTE: the order matters here greatly: ecosystem services
# need to come first.
N_ES = dm.count(category="ecosystemservices")
N_BD = dm.count(category="biodiversity")
WEIGHTS = [N_BD / N_ES] * N_ES + [1.0] * N_BD

# Define a group of rasters that need to be rescaled (normalized) for the
# following reasons:
#
# carbon_sequestration.tif = Values can have negative values (area is a carbon
#                            source instead of sink)
NORMALIZED_DATASETS = {"carbon_sequestration.tif":
                       "carbon_sequestration_rescaled.tif"}


# Define a group of rasters that need to be smoothed (9x9 median filter + log transformation)
# for the following reasons:
#
# pop_density_v5.tif = smooth the differences in population density so that the raster is
#                      more reasonable proxy for land cost
#
SMOOTHED_DATASETS = {"pop_density_v5.tif":
                     "pop_density_v5_smoothed.tif"}

# PROJECT RULES ----------------------------------------------------------------

rule all:
    input:
        "analyses/comparison/cross_correlation.csv",
        "analyses/comparison/cross_jaccard.csv",
        "analyses/comparison/cross_mcs.csv",
        "analyses/comparison/nuts2_rank_variation.shp"

# ## Data pre-processing ---------------------------------------------------------
#
# rule build_data_coverage:
#     input:
#         expand("data/processed/features/provide/{dataset}/{dataset}.tif", dataset=PROVIDE_DATASETS) + \
#         expand("data/processed/features/datadryad/forest_production_europe/{dataset}.tif", dataset=DATADRYAD_DATASETS)
#     output:
#         "data/processed/common/data_coverage.tif"
#     log:
#         "logs/build_data_coverage.log"
#     message:
#         "Building data coverage raster..."
#     run:
#         # In the first round, get the shape of the rasters and construct
#         # a np.array to hold the binary data maskas. NOTE: the assumption
#         # here is the all data have the same dimension
#         dims = None
#         masks = None
#         profile = None
#         for i, s_raster in enumerate(input):
#             with rasterio.open(s_raster) as in_src:
#                 if i == 0:
#                     dims = in_src.shape
#                     # Also store template profile
#                     profile = in_src.profile
#                     masks = np.zeros((len(input), dims[0], dims[1]), dtype=np.uint16)
#                 else:
#                     if in_src.shape != dims:
#                         llogger.warning(" WARNING: Dimensions for {} don't match, skipping".format(s_raster))
#                         continue
#                 # Place a binary version of the current raster mask into the
#                 # container
#                 llogger.info(" [{0}/{1}] Reading mask from {2}".format(i+1, len(input), s_raster))
#                 mask = in_src.read_masks(1)
#                 # We're using the GDAL type mask where any values above from
#                 # 0 are actual data values. Convert (in-place) to binary.
#                 np.place(mask, mask > 0, 1)
#                 masks[i] = mask
#
#         llogger.info(" Summing mask information and saving to {}".format(output))
#         data_coverage = masks.sum(axis=0)
#         # Write the product.
#         profile.update(dtype=rasterio.uint16, compress="DEFLATE")
#
#         with rasterio.open(output[0], 'w', **profile) as dst:
#             dst.write(data_coverage.astype(rasterio.uint16), 1)
#
rule preprocess_nuts_level0_data:
    input:
        shp=NUTS_LEVEL0_DATA
    output:
        reprojected=temp([path.replace("external", "interim/reprojected") for path in NUTS_LEVEL0_DATA]),
        enhanced=temp([path.replace("external", "interim/enhanced") for path in NUTS_LEVEL0_DATA]),
        processed=[path.replace("external", "processed") for path in NUTS_LEVEL0_DATA]
    log:
        "logs/preprocess_nuts_level0_data.log"
    message:
        "Pre-processing NUTS level 0 data..."
    run:
        llogger = utils.get_local_logger("pprocess_nuts0", log[0])
        # Read in the bounds as used in harmonize_data rule
        bleft = PROJECT_EXTENT["left"] + OFFSET[0]
        bbottom = PROJECT_EXTENT["bottom"] + OFFSET[1]
        bright = PROJECT_EXTENT["right"] + OFFSET[2]
        btop = PROJECT_EXTENT["top"] + OFFSET[3]
        bounds = "{0} {1} {2} {3}".format(bleft, bbottom, bright, btop)
        # Reproject to EPSG:3035 from EPSG:4258
        input_shp = utils.pick_from_list(input.shp, ".shp")
        reprojected_shp = utils.pick_from_list(output.reprojected, ".shp")
        cmd_str = 'ogr2ogr {0} -t_srs "EPSG:{1}" {2}'.format(reprojected_shp, PROJECT_CRS, input_shp)
        shell(cmd_str)
        llogger.info("Reprojected NUTS level 0 data from EPSG:4258 to EPSG:3035")
        llogger.debug(cmd_str)

        # NUTS 0 data has "NUTS_ID" field, but it's character. Convert to
        # integer for raserization
        enhanced_shp = utils.pick_from_list(output.enhanced, ".shp")
        with fiona.drivers():
            with fiona.open(reprojected_shp) as source:
                meta = source.meta
                meta['schema']['geometry'] = 'Polygon'
                # Insert new fields
                meta['schema']['properties']['ID'] = 'int'
                meta['schema']['properties']['mask'] = 'int'

                ID = 1
                with fiona.open(enhanced_shp, 'w', **meta) as sink:
                    # Loop over features
                    for f in source:
                        f['properties']['ID'] = ID
                        ID += 1
                        # Create a mask ID (same for each feature) that can
                        # later be used in creating a mask.
                        f['properties']['mask'] = 1
                        # Write the record out.
                        sink.write(f)

        # Clip shapefile using ogr2ogr, syntax:
        # ogr2ogr output.shp input.shp -clipsrc <left> <bottom> <right> <top>
        processed_shp = utils.pick_from_list(output.processed, ".shp")
        # Do 2 things at the same time:
        #  1. Select a subset of counties (defined by params.countries)
        #  2. Clip output to an extent (given by bounds)
        # Build the -where clause for ogr2ogr
        where_clause = "NUTS_ID IN ({})".format(", ".join(["'" + item + "'" for item in PROJECT_COUNTRIES]))
        shell('ogr2ogr -where "{where_clause}" {processed_shp} {enhanced_shp} -clipsrc {bounds}')
        llogger.debug("Clipped NUTS data to analysis bounds: {}".format(bounds))
        llogger.debug("Selected only a subset of eurostat countries:")
        llogger.debug(" " + ", ".join(PROJECT_COUNTRIES))
        llogger.debug("Resulting file: {}".format(processed_shp))

rule preprocess_nuts_level2_data:
    input:
        shp=NUTS_LEVEL2_DATA
    output:
        reprojected=temp([path.replace("external", "interim/reprojected") for path in NUTS_LEVEL2_DATA]),
        clipped=temp([path.replace("external", "interim/clipped") for path in NUTS_LEVEL2_DATA]),
        processed=[path.replace("external", "processed") for path in NUTS_LEVEL2_DATA]
    log:
        "logs/preprocess_nuts_level2_data.log"
    message:
        "Pre-processing NUTS level 2 data..."
    run:
        llogger = utils.get_local_logger("pprocess_nuts2", log[0])
        # Read in the bounds as used in harmonize_data rule
        bleft = PROJECT_EXTENT["left"] + OFFSET[0]
        bbottom = PROJECT_EXTENT["bottom"] + OFFSET[1]
        bright = PROJECT_EXTENT["right"] + OFFSET[2]
        btop = PROJECT_EXTENT["top"] + OFFSET[3]
        bounds = "{0} {1} {2} {3}".format(bleft, bbottom, bright, btop)
        # Reproject to EPSG:3035 from EPSG:4258 and clip
        input_shp = utils.pick_from_list(input.shp, ".shp")
        reprojected_shp = utils.pick_from_list(output.reprojected, ".shp")
        cmd_str = 'ogr2ogr {0} -t_srs "EPSG:{1}" {2}'.format(reprojected_shp, PROJECT_CRS, input_shp)
        shell(cmd_str)
        llogger.debug("Reprojected NUTS level 2 data from EPSG:4258 to EPSG:3035")
        llogger.debug(cmd_str)

        # Clip shapefile using ogr2ogr, syntax:
        # ogr2ogr output.shp input.shp -clipsrc <left> <bottom> <right> <top>
        clipped_shp = utils.pick_from_list(output.clipped, ".shp")
        # Clip output to an extent (given by bounds)
        shell('ogr2ogr {clipped_shp} {reprojected_shp} -clipsrc {bounds}')

        # The Pre-processing steps need to be done:
        #  1. Tease apart country code from field NUTS_ID
        #  2. Create a running ID field that can be used as value in the
        #     rasterized version
        processed_shp = utils.pick_from_list(output.processed, ".shp")
        with fiona.drivers():
            with fiona.open(clipped_shp) as source:
                meta = source.meta
                meta['schema']['geometry'] = 'Polygon'
                # Insert new fields
                meta['schema']['properties']['ID'] = 'int'
                meta['schema']['properties']['country'] = 'str'

                ID = 1
                with fiona.open(processed_shp, 'w', **meta) as sink:
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

        llogger.debug("Clipped NUTS level 2 data to analysis bounds: {}".format(bounds))
        llogger.debug("Selected only a subset of eurostat countries:")
        llogger.debug(" " + ", ".join(PROJECT_COUNTRIES))
        llogger.debug("Resulting file: {}".format(processed_shp))

rule rasterize_nuts_level0_data:
    # Effectively, create a 1) land mask and 2) common data mask
    input:
        rules.preprocess_nuts_level0_data.output.processed
    output:
        land_mask=utils.pick_from_list(rules.preprocess_nuts_level0_data.output.processed, ".shp").replace(".shp", ".tif"),
        data_mask=utils.pick_from_list(rules.preprocess_nuts_level0_data.output.processed, ".shp").replace(".shp", "_data_mask.tif"),
    log:
        "logs/rasterize_nuts_level0_data.log"
    message:
        "Rasterizing NUTS level 0 data..."
    run:
        llogger = utils.get_local_logger("rasterize_nuts0", log[0])
        input_shp = utils.pick_from_list(input, ".shp")
        layer_shp = os.path.basename(input_shp).replace(".shp", "")

        # Construct extent
        bounds = "{0} {1} {2} {3}".format(PROJECT_EXTENT["left"],
                                          PROJECT_EXTENT["bottom"],
                                          PROJECT_EXTENT["right"],
                                          PROJECT_EXTENT["top"])
        # 1) Rasterize land mask
        cmd_str = "gdal_rasterize -l {} ".format(layer_shp) + \
                  "-a ID -tr 1000 1000 -te {} ".format(bounds) + \
                  "-ot Int16 -a_nodata -32768 -co COMPRESS=DEFLATE " + \
                  "{0} {1}".format(input_shp, output.land_mask)
        llogger.debug(cmd_str)
        for line in utils.process_stdout(shell(cmd_str, read=True)):
            llogger.debug(line)
        # 2) Rasterize common data mask
        cmd_str = "gdal_rasterize -l {} ".format(layer_shp) + \
                  "-a mask -tr 1000 1000 -te {} ".format(bounds) + \
                  "-ot Int8 -a_nodata -128 -co COMPRESS=DEFLATE " + \
                  "{0} {1}".format(input_shp, output.data_mask)
        llogger.debug(cmd_str)
        for line in utils.process_stdout(shell(cmd_str, read=True)):
            llogger.debug(line)

rule rasterize_nuts_level2_data:
    input:
        rules.preprocess_nuts_level2_data.output.processed
    output:
        utils.pick_from_list(rules.preprocess_nuts_level2_data.output.processed, ".shp").replace(".shp", ".tif")
    log:
        "logs/rasterize_nuts_level2_data.log"
    message:
        "Rasterizing NUTS level 2 data..."
    run:
        llogger = utils.get_local_logger("rasterize_nuts2", log[0])
        input_shp = utils.pick_from_list(input, ".shp")
        layer_shp = os.path.basename(input_shp).replace(".shp", "")

        # Construct extent
        bounds = "{0} {1} {2} {3}".format(PROJECT_EXTENT["left"],
                                          PROJECT_EXTENT["bottom"],
                                          PROJECT_EXTENT["right"],
                                          PROJECT_EXTENT["top"])
        # Rasterize
        cmd_str = "gdal_rasterize -l {} ".format(layer_shp) + \
                  "-a ID -tr 1000 1000 -te {} ".format(bounds) + \
                  "-ot Int16 -a_nodata -32768 -co COMPRESS=DEFLATE " + \
                  "{0} {1}".format(input_shp, output[0])
        llogger.debug(cmd_str)
        for line in utils.process_stdout(shell(cmd_str, read=True)):
            llogger.debug(line)

rule clip_udr_data:
    input:
        external=UDR_SRC_DATASETS,
        clip_shp=utils.pick_from_list(rules.preprocess_nuts_level0_data.output.processed, ".shp")
    output:
        clipped=[path.replace("external", "processed/features") for path in UDR_SRC_DATASETS]
    log:
        "logs/clip_udr_data.log"
    message:
        "Clipping UDR data..."
    run:
        llogger = utils.get_local_logger("clip_udr_data", log[0])
        nsteps = len(input.external)
        for i, s_raster in enumerate(input.external):
            # Target raster
            clipped_raster = s_raster.replace("external", "processed/features")
            prefix = utils.get_iteration_prefix(i+1, nsteps)

            llogger.info("{0} Clipping dataset {1}".format(prefix, s_raster))
            llogger.debug("{0} Target dataset {1}".format(prefix, clipped_raster))
            # Clip data. NOTE: UDR species rasters do not have a SRS defined,
            # but they are in EPSG:3035
            cmd_str = 'gdalwarp -s_srs EPSG:3035 -t_srs EPSG:3035 -cutline {0} {1} {2} -co COMPRESS=DEFLATE'.format(input.clip_shp, s_raster, clipped_raster)
            for line in utils.process_stdout(shell(cmd_str, read=True), prefix=prefix):
                llogger.debug(line)

rule harmonize_data:
    input:
        external=DATADRYAD_SRC_DATASETS+PROVIDE_SRC_DATASETS+EEA_SRC_DATASETS,
        like_raster=[path for path in DATADRYAD_SRC_DATASETS if "woodprod_average" in path][0],
        clip_shp=utils.pick_from_list(rules.preprocess_nuts_level0_data.output.processed, ".shp")
    output:
        # NOTE: UDR_SRC_DATASETS do not need to processed
        warped=temp([path.replace("external", "interim/warped") for path in DATADRYAD_SRC_DATASETS+PROVIDE_SRC_DATASETS+EEA_SRC_DATASETS]),
        harmonized=[path.replace("external", "processed/features") for path in DATADRYAD_SRC_DATASETS+PROVIDE_SRC_DATASETS+EEA_SRC_DATASETS]
    log:
        "logs/harmonize_data.log"
    message:
        "Harmonizing datasets..."
    run:
        llogger = utils.get_local_logger("harmonize_data", log[0])
        nsteps = len(input.external)
        for i, s_raster in enumerate(input.external):
            ## WARP
            # Target raster
            warped_raster = s_raster.replace("external", "interim/warped")
            # No need to process the snap raster, just copy it
            prefix = utils.get_iteration_prefix(i+1, nsteps)
            if s_raster == input.like_raster:
                llogger.info("{0} Copying dataset {1}".format(prefix, s_raster))
                llogger.debug("{0} Target dataset {1}".format(prefix, warped_raster))
                ret = shell("cp {s_raster} {warped_raster}", read=True)
            else:
                llogger.info("{0} Warping dataset {1}".format(prefix, s_raster))
                llogger.debug("{0} Target dataset {1}".format(prefix, warped_raster))
                ret = shell("rio warp " + s_raster + " --like " + input.like_raster + \
                            " " + warped_raster + " --dst-crs " + str(PROJECT_CRS) + \
                            " --res " + str(PROJECT_RES) + \
                            " --co 'COMPRESS=DEFLATE' --threads {threads}")
            for line in utils.process_stdout(ret, prefix=prefix):
                llogger.debug(line)
            ## CLIP
            harmonized_raster = warped_raster.replace("data/interim/warped", "data/processed/features")
            llogger.info("{0} Clipping dataset {1}".format(prefix, warped_raster))
            llogger.debug("{0} Target dataset {1}".format(prefix, harmonized_raster))
            cmd_str = "gdalwarp -cutline {0} {1} {2} -co COMPRESS=DEFLATE".format(input.clip_shp, warped_raster, harmonized_raster)
            for line in utils.process_stdout(shell(cmd_str, read=True), prefix=prefix):
                llogger.debug(line)

            # Rescale (normalize) dataset if needed
            org_raster = os.path.basename(harmonized_raster)
            if org_raster in NORMALIZED_DATASETS.keys():
                rescaled_raster = harmonized_raster.replace(org_raster,
                                                            NORMALIZED_DATASETS[org_raster])
                llogger.info("{0} Rescaling dataset {1}".format(prefix, harmonized_raster))
                llogger.debug("{0} Target dataset {1}".format(prefix, rescaled_raster))
                spatutils.rescale_raster(harmonized_raster, rescaled_raster,
                                         method="normalize",
                                         only_positive=True, verbose=False)
                os.remove(harmonized_raster)
                llogger.debug("{0} Renaming dataset {1} to {2}".format(prefix, rescaled_raster, harmonized_raster))
                os.rename(rescaled_raster, harmonized_raster)
                harmonized_raster = rescaled_raster

            # Finally, smooth dataset if needed
            org_raster = os.path.basename(harmonized_raster)
            if org_raster in SMOOTHED_DATASETS.keys():
                smoothed_raster = harmonized_raster.replace(org_raster,
                                                            SMOOTHED_DATASETS[org_raster])
                llogger.info("{0} Smoothing dataset {1}".format(prefix, harmonized_raster))
                llogger.debug("{0} Target dataset {1}".format(prefix, smoothed_raster))
                spatutils.smooth_raster(harmonized_raster, smoothed_raster,
                                        log_transform=True, method="medfilt",
                                        verbose=False)
                os.remove(harmonized_raster)
                llogger.debug("{0} Renaming dataset {1} to {2}".format(prefix, smoothed_raster, harmonized_raster))
                os.rename(smoothed_raster, harmonized_raster)



rule ol_normalize_data:
    input:
        rules.harmonize_data.output.harmonized+UDR_SRC_DATASETS
    output:
        [path.replace("processed/features", "processed/features_ol_normalized") for path in rules.harmonize_data.output.harmonized+UDR_SRC_DATASETS]
    log:
        "logs/ol_normalize_data.log"
    message:
        "Normalizing data based on occurrence levels..."
    run:
        llogger = utils.get_local_logger("ol_normalize_data", log[0])
        for i, s_raster in enumerate(input):
            prefix = utils.get_iteration_prefix(i+1, len(input))
            llogger.info("{0} (OL) Normalizing dataset {1}".format(prefix, s_raster))
            # NOTE: looping over input and output only works if they have
            # exactly the same definition. Otherwise order may vary.
            spatutils.rescale_raster(input[i], output[i], method="ol_normalize",
                                     fill_w_zeros=True, logger=llogger)

rule rescale_data:
    input:
        rules.ol_normalize_data.output
    output:
        [path.replace("features_ol_normalized", "features_ol_normalized_rescaled") for path in rules.ol_normalize_data.output]
    message:
        "Normalizing data..."
    run:
        for i, s_raster in enumerate(input):
            # No need to process the snap raster
            llogger.info(" [{0}/{1}] Rescaling dataset {2}".format(i+1, len(input), s_raster))
            # NOTE: looping over input and output only works if they have
            # exactly the same definition. Otherwise order may vary.
            spatutils.rescale_raster(input[i], output[i], method="normalize",
                                     verbose=False)

## Set up, run and post-process analyses --------------------------------------

# RWR -------------------------------------------------------------------------

rule prioritize_rwr:
    input:
        all=rules.harmonize_data.output.harmonized[:-1]+UDR_SRC_DATASETS,
        es=rules.harmonize_data.output.harmonized[:-1],
        bd=UDR_SRC_DATASETS,
        cost=rules.harmonize_data.output.harmonized[-1]
    output:
        all="analyses/RWR/rwr_all.tif",
        all_w="analyses/RWR/rwr_all_weights.tif",
        all_w_c="analyses/RWR/rwr_all_weights_costs.tif",
        es="analyses/RWR/rwr_es.tif",
        es_c="analyses/RWR/rwr_es_costs.tif",
        bd="analyses/RWR/rwr_bd.tif",
        bd_c="analyses/RWR/rwr_bd_costs.tif"
    log:
        all="logs/calculate_rwr_all.log",
        all_w="logs/calculate_rwr_all_weights.log",
        all_w_c="logs/calculate_rwr_all_weights_costs.log",
        es="logs/calculate_rwr_es.log",
        es_c="logs/calculate_rwr_es_costs.log",
        bd="logs/calculate_rwr_bd.log",
        bd_c="logs/calculate_rwr_bd_costs.log"
    message:
        "Calculating RWR..."
    run:
        # All without weights
        llogger = utils.get_local_logger("calculate_rwr_all", log.all)
        rwr.calculate_rwr(input.all, output.all, logger=llogger)

        # All with weights
        llogger = utils.get_local_logger("calculate_rwr_all_weights", log.all_w)
        rwr.calculate_rwr(input.all, output.all_w, weights=WEIGHTS,
                          logger=llogger)

        # All with weights and costs
        llogger = utils.get_local_logger("calculate_rwr_all_weights_costs",
                                         log.all_w_c)
        rwr.calculate_rwr(input.all, output.all_w_c, cost_raster=input.cost,
                          weights=WEIGHTS, logger=llogger)

        # ES
        llogger = utils.get_local_logger("calculate_rwr_es", log.es)
        rwr.calculate_rwr(input.es, output.es, logger=llogger)

        # ES with costs
        llogger = utils.get_local_logger("calculate_rwr_es_costs", log.es_c)
        rwr.calculate_rwr(input.es, output.es_c, cost_raster=input.cost,
                          logger=llogger)

        # BD
        llogger = utils.get_local_logger("calculate_rwr_bd", log.bd)
        rwr.calculate_rwr(input.bd, output.bd, logger=llogger)

        # BD with costs
        llogger = utils.get_local_logger("calculate_rwr_bd_costs", log.bd_c)
        rwr.calculate_rwr(input.bd, output.bd_c, cost_raster=input.cost,
                          logger=llogger)

rule postprocess_rwr:
    input:
        all=rules.prioritize_rwr.output.all,
        all_w=rules.prioritize_rwr.output.all_w,
        all_w_c=rules.prioritize_rwr.output.all_w_c,
        es=rules.prioritize_rwr.output.es,
        es_c=rules.prioritize_rwr.output.es_c,
        bd=rules.prioritize_rwr.output.bd,
        bd_c=rules.prioritize_rwr.output.bd_c,
        plu=utils.pick_from_list(rules.preprocess_nuts_level2_data.output.processed,
                                 ".shp")
    output:
        all="analyses/RWR/rwr_all_stats.geojson",
        all_w="analyses/RWR/rwr_all_weights_stats.geojson",
        all_w_c="analyses/RWR/rwr_all_weights_costs_stats.geojson",
        es="analyses/RWR/rwr_es_stats.geojson",
        es_c="analyses/RWR/rwr_es_costs_stats.geojson",
        bd="analyses/RWR/rwr_bd_stats.geojson",
        bd_c="analyses/RWR/rwr_bd_costs_stats.geojson"
    log:
        all="logs/postprocess_rwr_all.log",
        all_w="logs/postprocess_rwr_all_weights.log",
        all_w_c="logs/postprocess_rwr_all_weights_costs.log",
        es="logs/postprocess_rwr_es.log",
        es_c="logs/postprocess_rwr_es_costs.log",
        bd="logs/postprocess_rwr_bd.log",
        bd_c="logs/postprocess_rwr_bd_costs.log"
    message:
        "Post-processing RWR results..."
    run:
        llogger = utils.get_local_logger("calculate_rwr_all", log.all)
        llogger.info(" [1/7] Post-processing {}".format(input.all))
        shell("fio cat {input.plu} | rio zonalstats -r {input.all} > {output.all}")

        llogger = utils.get_local_logger("calculate_rwr_all_weights", log.all_w)
        llogger.info(" [2/7] Post-processing {}".format(input.all_w))
        shell("fio cat {input.plu} | rio zonalstats -r {input.all_w} > {output.all_w}")

        llogger = utils.get_local_logger("calculate_rwr_all_weights_costs",
                                         log.all_w)
        llogger.info(" [3/7] Post-processing {}".format(input.all_w_c))
        shell("fio cat {input.plu} | rio zonalstats -r {input.all_w_c} > {output.all_w_c}")

        llogger = utils.get_local_logger("calculate_rwr_es", log.es)
        llogger.info(" [4/7] Post-processing {}".format(input.es))
        shell("fio cat {input.plu} | rio zonalstats -r {input.es} > {output.es}")

        llogger = utils.get_local_logger("calculate_rwr_es_costs", log.es_c)
        llogger.info(" [5/7] Post-processing {}".format(input.es_c))
        shell("fio cat {input.plu} | rio zonalstats -r {input.es_c} > {output.es_c}")

        llogger = utils.get_local_logger("calculate_rwr_bd", log.bd)
        llogger.info(" [6/7] Post-processing {}".format(input.bd))
        shell("fio cat {input.plu} | rio zonalstats -r {input.bd} > {output.bd}")

        llogger = utils.get_local_logger("calculate_rwr_bd_costs", log.bd_c)
        llogger.info(" [7/7] Post-processing {}".format(input.bd_c))
        shell("fio cat {input.plu} | rio zonalstats -r {input.bd_c} > {output.bd_c}")

rule expand_rwr_coverage:
    input:
        zon_all_w="analyses/zonation/priocomp/04_abf_all_wgt/04_abf_all_wgt_out/04_abf_all_wgt.rank.compressed.tif",
        zon_all_w_c="analyses/zonation/priocomp/06_abf_all_wgt_cst/06_abf_all_wgt_cst_out/06_abf_all_wgt_cst.rank.compressed.tif",
        rwr_all_w=rules.prioritize_rwr.output.all_w,
        rwr_all_w_c=rules.prioritize_rwr.output.all_w_c
    output:
        all_w="analyses/RWR/rwr_all_weights_expanded.tif",
        all_w_c="analyses/RWR/rwr_all_weights_costs_expanded.tif"
    log:
        all_w="logs/expand_rwr_all_w.log",
        all_w_c="logs/expand_rwr_all_w_c.log"
    message:
        "Post-processing (expanding) RWR results..."
    run:
        llogger = utils.get_local_logger("expand_rwr_w", log.all_w)
        coverage.expand_value_coverage(input.rwr_all_w, input.zon_all_w,
                                       output.all_w, logger=llogger)

        llogger = utils.get_local_logger("expand_rwr_w_c", log.all_w_c)
        coverage.expand_value_coverage(input.rwr_all_w_c, input.zon_all_w_c,
                                       output.all_w_c, logger=llogger)

# # Zonation ---------------------------------------------------------------------
#
# rule generate_zonation_project:
#     input:
#         expand("data/processed/features/provide/{dataset}/{dataset}.tif", dataset=PROVIDE_DATASETS) + \
#         expand("data/processed/features/datadryad/forest_production_europe/{dataset}.tif", dataset=DATADRYAD_DATASETS),
#         "data/processed/nuts/NUTS_RG_01M_2013/level2/NUTS_RG_01M_2013_level2_subset.tif"
#     output:
#         "analyses/zonation/priocomp"
#     message:
#         "Generating Zonation project..."
#     script:
#         # NOTE: Currently there's a lot of things hardcoded here...
#         "src/zonation/01_create_zonation_project.R"

rule create_zon_coverage:
    input:
        all_w="analyses/zonation/priocomp/04_abf_wgt/04_abf_wgt_out/04_abf_wgt.rank.compressed.tif",
        es="analyses/zonation/priocomp/06_abf_es/06_abf_es_out/06_abf_es.rank.compressed.tif",
        bd="analyses/zonation/priocomp/08_abf_bd/08_abf_bd_out/08_abf_bd.rank.compressed.tif",
    output:
        all_w="analyses/zonation/priocomp/04_abf_wgt/04_abf_wgt_out/04_abf_wgt_data_coverage.tif",
        es="analyses/zonation/priocomp/06_abf_es/06_abf_es_out/06_abf_es_data_coverage.tif",
        bd="analyses/zonation/priocomp/08_abf_bd/08_abf_bd_out/08_abf_bd_data_coverage.tif"
    log:
        all_w="logs/postprocess_zon_04_abf_wgt.log",
        es="logs/postprocess_zon_06_es.log",
        bd="logs/postprocess_zon_08_bd.log"
    message:
        "Post-processing ZON results..."
    run:
        llogger = utils.get_local_logger("coverage_all_w", log.all_w)
        coverage.create_value_coverage(input.all_w, output.all_w,
                                       logger=llogger)
        llogger = utils.get_local_logger("coverage_es", log.es)
        coverage.create_value_coverage(input.es, output.es, logger=llogger)
        llogger = utils.get_local_logger("coverage_bd", log.bd)
        coverage.create_value_coverage(input.es, output.bd, logger=llogger)

rule expand_zon_coverage:
    input:
        all_w="analyses/zonation/priocomp/04_abf_all_wgt/04_abf_all_wgt_out/04_abf_all_wgt.rank.compressed.tif",
        all_w_c="analyses/zonation/priocomp/06_abf_all_wgt_cst/06_abf_all_wgt_cst_out/06_abf_all_wgt_cst.rank.compressed.tif",
        es="analyses/zonation/priocomp/08_abf_es/08_abf_es_out/08_abf_es.rank.compressed.tif",
        es_c="analyses/zonation/priocomp/10_abf_es_cst/10_abf_es_cst_out/10_abf_es_cst.rank.compressed.tif",
        bd="analyses/zonation/priocomp/12_abf_bd/12_abf_bd_out/12_abf_bd.rank.compressed.tif",
        bd_c="analyses/zonation/priocomp/14_abf_bd_cst/14_abf_bd_cst_out/14_abf_bd_cst.rank.compressed.tif"
    output:
        es_expanded="analyses/zonation/priocomp/08_abf_es/08_abf_es_out/08_abf_es.rank_expanded.compressed.tif",
        es_c_expanded="analyses/zonation/priocomp/10_abf_es_cst/10_abf_es_cst_out/10_abf_es_cst.rank_expanded.compressed.tif",
        bd_expanded="analyses/zonation/priocomp/12_abf_bd/12_abf_bd_out/12_abf_bd.rank_expanded.compressed.tif",
        bd_c_expanded="analyses/zonation/priocomp/14_abf_bd_cst/14_abf_bd_cst_out/14_abf_bd_cst.rank_expanded.compressed.tif"
    log:
        es="logs/expand_zon_08_es.log",
        es_c="logs/expand_zon_10_es_cost.log",
        bd="logs/expand_zon_12_bd.log",
        bd_c="logs/expand_zon_14_bd_cost.log",
    message:
        "Post-processing (expanding) ZON results..."
    run:
        llogger = utils.get_local_logger("expand_es", log.es)
        coverage.expand_value_coverage(input.es, input.all_w, output.es_expanded,
                                       logger=llogger)
        llogger = utils.get_local_logger("expand_es_cost", log.es_c)
        coverage.expand_value_coverage(input.es_c, input.all_w_c, output.es_c_expanded,
                                       logger=llogger)

        llogger = utils.get_local_logger("expand_bd", log.bd)
        coverage.expand_value_coverage(input.bd, input.all_w, output.bd_expanded,
                                       logger=llogger)
        llogger = utils.get_local_logger("expand_bd_cost", log.bd_c)
        coverage.expand_value_coverage(input.bd_c, input.all_w_c, output.bd_c_expanded,
                                       logger=llogger)

rule match_zon_esbd_coverages:
    input:
        es="analyses/zonation/priocomp/08_abf_es/08_abf_es_out/08_abf_es.rank.compressed.tif",
        es_c="analyses/zonation/priocomp/10_abf_es_cst/10_abf_es_cst_out/10_abf_es_cst.rank.compressed.tif",
        bd="analyses/zonation/priocomp/12_abf_bd/12_abf_bd_out/12_abf_bd.rank.compressed.tif",
        bd_c="analyses/zonation/priocomp/14_abf_bd_cst/14_abf_bd_cst_out/14_abf_bd_cst.rank.compressed.tif"
    output:
        es_matched="analyses/zonation/priocomp/08_abf_es/08_abf_es_out/08_abf_es.rank_bd_matched.compressed.tif",
        es_c_matched="analyses/zonation/priocomp/10_abf_es_cst/10_abf_es_cst_out/10_abf_es_cst.rank_bd_matched.compressed.tif",
        bd_matched="analyses/zonation/priocomp/12_abf_bd/12_abf_bd_out/12_abf_bd.rank_es_matched.compressed.tif",
        bd_c_matched="analyses/zonation/priocomp/14_abf_bd_cst/14_abf_bd_cst_out/14_abf_bd_cst.rank_es_matched.compressed.tif"
    log:
        es="logs/match_zon_08_es_bd.log",
        es_c="logs/match_zon_10_es_bd_cost.log",
        bd="logs/match_zon_12_bd_es.log",
        bd_c="logs/match_zon_14_bd_es_cost.log"
    message:
        "Post-processing (matching) ZON results..."
    run:
        llogger = utils.get_local_logger("match_es_bd", log.es)
        coverage.expand_value_coverage(input.es, input.bd, output.es_matched,
                                       logger=llogger)
        llogger = utils.get_local_logger("match_es_bd_cost", log.es_c)
        coverage.expand_value_coverage(input.es_c, input.bd_c, output.es_c_matched,
                                       logger=llogger)

        llogger = utils.get_local_logger("match_bd_es", log.bd)
        coverage.expand_value_coverage(input.bd, input.es, output.bd_matched,
                                       logger=llogger)
        llogger = utils.get_local_logger("match_bd_es_cost", log.bd_c)
        coverage.expand_value_coverage(input.bd_c, input.es_c, output.bd_c_matched,
                                       logger=llogger)

# ILP ------------------------------------------------------------------------

rule prioritize_ilp_all:
    input:
        all=rules.harmonize_data.output.harmonized[:-1]+UDR_SRC_DATASETS,
        cost=rules.harmonize_data.output.harmonized[-1]
    output:
        all_w="analyses/ILP/ilp_all_weights.tif",
        all_w_c="analyses/ILP/ilp_all_weights_costs.tif"
    log:
        all_w="logs/prioritize_ilp_all_weights.log",
        all_w_c="logs/prioritize_ilp_all_weights_costs.log"
    message:
        "Optimizing ALL with Gurobi..."
    run:
        # Without weights
        #llogger = utils.get_local_logger("optimize_gurobi_all", log.all)
        #gurobi.prioritize_gurobi(input.all, output.all, logger=llogger,
        #                         ol_normalize=True, save_intermediate=False,
        #                         verbose=True)
        # With weights
        llogger = utils.get_local_logger("optimize_gurobi_all_weights",
                                         log.all_w)
        gurobi.prioritize_gurobi(input.all, output.all_w, logger=llogger,
                                 ol_normalize=True, weights=WEIGHTS,
                                 step=0.01, save_intermediate=True,
                                 verbose=True)
        # With weights and costs
        llogger = utils.get_local_logger("optimize_gurobi_all_weights_costs",
                                         log.all_w_c)
        gurobi.prioritize_gurobi(input.all, output.all_w_c, input.cost, logger=llogger,
                                 ol_normalize=True, weights=WEIGHTS,
                                 step=0.01, save_intermediate=True,
                                 verbose=True)

rule prioritize_ilp_es:
    input:
        # pop_density_v5.tif is the last item in harmonized datasets
        es=rules.harmonize_data.output.harmonized[:-1],
        cost=rules.harmonize_data.output.harmonized[-1]
    output:
        es="analyses/ILP/ilp_es.tif",
        es_c="analyses/ILP/ilp_es_costs.tif"
    log:
        es="logs/prioritize_ilp_es.log",
        es_c="logs/prioritize_ilp_es_costs.log"
    message:
        "Optimizing ES with Gurobi..."
    run:
        # Without costs
        llogger = utils.get_local_logger("optimize_gurobi_es", log.es)
        gurobi.prioritize_gurobi(input.es, output.es, logger=llogger,
                                 ol_normalize=True, step=0.01,
                                 save_intermediate=True, verbose=True)
        # With costs
        llogger = utils.get_local_logger("optimize_gurobi_es_c", log.es_c)
        gurobi.prioritize_gurobi(input.es, output.es_c, input.cost, logger=llogger,
                                 ol_normalize=True, step=0.01,
                                 save_intermediate=True, verbose=True)

rule prioritize_ilp_bd:
    input:
        bd=UDR_SRC_DATASETS,
        cost=rules.harmonize_data.output.harmonized[-1]
    output:
        bd="analyses/ILP/ilp_bd.tif",
        bd_c="analyses/ILP/ilp_bd_costs.tif"
    log:
        bd="logs/prioritize_ilp_bd.log",
        bd_c="logs/prioritize_ilp_bd_costs.log"
    message:
        "Optimizing BD with Gurobi..."
    run:
        # Without costs
        llogger = utils.get_local_logger("optimize_gurobi_bd", log.bd)
        gurobi.prioritize_gurobi(input.bd, output.bd, logger=llogger,
                                 ol_normalize=True, step=0.01,
                                 save_intermediate=True, verbose=True)

        # With costs
        llogger = utils.get_local_logger("optimize_gurobi_bd_c", log.bd_c)
        gurobi.prioritize_gurobi(input.bd, output.bd_c, input.cost, logger=llogger,
                                 ol_normalize=True, step=0.01,
                                 save_intermediate=True, verbose=True)

rule prioritize_ilp:
    input:
        rules.prioritize_ilp_all.output.all_w,
        rules.prioritize_ilp_es.output.es,
        rules.prioritize_ilp_bd.output.bd
    output:
        "analyses/ILP/ilp_done.txt"

rule postprocess_ilp:
    input:
        all_w=rules.prioritize_ilp_all.output.all_w,
        all_w_c=rules.prioritize_ilp_all.output.all_w_c,
        es=rules.prioritize_ilp_es.output.es,
        es_c=rules.prioritize_ilp_es.output.es_c,
        bd=rules.prioritize_ilp_bd.output.bd,
        bd_c=rules.prioritize_ilp_bd.output.bd_c,
        plu=utils.pick_from_list(rules.preprocess_nuts_level2_data.output.processed,
                                 ".shp")
    output:
        all_w="analyses/ILP/ilp_all_weights_stats.geojson",
        all_w_c="analyses/ILP/ilp_all_weights_costs_stats.geojson",
        es="analyses/ILP/ilp_es_stats.geojson",
        es_c="analyses/ILP/ilp_es_costs_stats.geojson",
        bd="analyses/ILP/ilp_bd_stats.geojson",
        bd_c="analyses/ILP/ilp_bd_costs_stats.geojson"
    log:
        all_w="logs/postprocess_ilp_all_weights.log",
        all_w_c="logs/postprocess_ilp_all_weights_costs.log",
        es="logs/postprocess_ilp_es.log",
        es_c="logs/postprocess_ilp_es_costs.log",
        bd="logs/postprocess_ilp_bd.log",
        bd_c="logs/postprocess_ilp_bd_costs.log"
    message:
        "Post-processing ILP results..."
    run:
        llogger = utils.get_local_logger("calculate_ilp_all_weights", log.all_w)
        llogger.info(" [1/6] Post-processing {}".format(input.all_w))
        shell("fio cat {input.plu} | rio zonalstats -r {input.all_w} > {output.all_w}")

        llogger = utils.get_local_logger("calculate_ilp_all_weights_costs", log.all_w_c)
        llogger.info(" [2/6] Post-processing {}".format(input.all_w_c))
        shell("fio cat {input.plu} | rio zonalstats -r {input.all_w_c} > {output.all_w_c}")

        llogger = utils.get_local_logger("calculate_ilp_es", log.es)
        llogger.info(" [3/6] Post-processing {}".format(input.es))
        shell("fio cat {input.plu} | rio zonalstats -r {input.es} > {output.es}")

        llogger = utils.get_local_logger("calculate_ilp_es_costs", log.es_c)
        llogger.info(" [4/6] Post-processing {}".format(input.es_c))
        shell("fio cat {input.plu} | rio zonalstats -r {input.es_c} > {output.es_c}")

        llogger = utils.get_local_logger("calculate_ilp_bd", log.bd)
        llogger.info(" [5/6] Post-processing {}".format(input.bd))
        shell("fio cat {input.plu} | rio zonalstats -r {input.bd} > {output.bd}")

        llogger = utils.get_local_logger("calculate_ilp_bd_costs", log.bd_c)
        llogger.info(" [6/6] Post-processing {}".format(input.bd_c))
        shell("fio cat {input.plu} | rio zonalstats -r {input.bd_c} > {output.bd_c}")

rule expand_ilp_coverage:
    input:
        template_w="analyses/zonation/priocomp/04_abf_all_wgt/04_abf_all_wgt_out/04_abf_all_wgt.rank.compressed.tif",
        template_w_c="analyses/zonation/priocomp/06_abf_all_wgt_cst/06_abf_all_wgt_cst_out/06_abf_all_wgt_cst.rank.compressed.tif",
        all_w=rules.prioritize_ilp_all.output.all_w,
        all_w_c=rules.prioritize_ilp_all.output.all_w_c,
    output:
        all_w="analyses/ILP/ilp_all_weights_expanded.tif",
        all_w_c="analyses/ILP/ilp_all_weights_costs_expanded.tif"
    log:
        all_w="logs/expand_ilp_all_w.log",
        all_w_c="logs/expand_ilp_all_w_c.log"
    message:
        "Post-processing (expanding) ILP results..."
    run:
        # With weights
        llogger = utils.get_local_logger("expand_ilp_w", log.all_w)
        coverage.expand_value_coverage(input.all_w, input.template_w,
                                       output.all_w, logger=llogger)
        # With weights and costs
        llogger = utils.get_local_logger("expand_ilp_w_c", log.all_w_c)
        coverage.expand_value_coverage(input.all_w_c, input.template_w,
                                       output.all_w_c, logger=llogger)

## Compare results ------------------------------------------------------------

rule compare_correlation:
    input:
        rules.prioritize_rwr.output.all_w,
        "analyses/zonation/priocomp/04_abf_all_wgt/04_abf_all_wgt_out/04_abf_all_wgt.rank.compressed.tif",
        rules.prioritize_ilp_all.output.all_w,
        rules.prioritize_rwr.output.all_w_c,
        "analyses/zonation/priocomp/06_abf_all_wgt_cst/06_abf_all_wgt_cst_out/06_abf_all_wgt_cst.rank.compressed.tif",
        rules.prioritize_ilp_all.output.all_w_c,
        rules.prioritize_rwr.output.es,
        "analyses/zonation/priocomp/08_abf_es/08_abf_es_out/08_abf_es.rank.compressed.tif",
        rules.prioritize_ilp_es.output.es,
        rules.prioritize_rwr.output.es_c,
        "analyses/zonation/priocomp/10_abf_es_cst/10_abf_es_cst_out/10_abf_es_cst.rank.compressed.tif",
        rules.prioritize_ilp_es.output.es_c,
        rules.prioritize_rwr.output.bd,
        "analyses/zonation/priocomp/12_abf_bd/12_abf_bd_out/12_abf_bd.rank.compressed.tif",
        rules.prioritize_ilp_bd.output.bd,
        rules.prioritize_rwr.output.bd_c,
        "analyses/zonation/priocomp/14_abf_bd_cst/14_abf_bd_cst_out/14_abf_bd_cst.rank.compressed.tif",
        rules.prioritize_ilp_bd.output.bd_c
    output:
        "analyses/comparison/cross_correlation.csv"
    log:
        "logs/compare_results_correlation.log"
    message:
        "Comparing results correlation using Kendall tau..."
    run:
        llogger = utils.get_local_logger("compare_correlation", log[0])
        correlations = similarity.cross_correlation(input, verbose=False,
                                                    logger=llogger)
        llogger.info("Saving results to {}".format(output[0]))
        correlations.to_csv(output[0], index=False)

rule compare_jaccard:
    input:
        rules.prioritize_rwr.output.all_w,
        "analyses/zonation/priocomp/04_abf_all_wgt/04_abf_all_wgt_out/04_abf_all_wgt.rank.compressed.tif",
        rules.prioritize_ilp_all.output.all_w,
        rules.prioritize_rwr.output.all_w_c,
        "analyses/zonation/priocomp/06_abf_all_wgt_cst/06_abf_all_wgt_cst_out/06_abf_all_wgt_cst.rank.compressed.tif",
        rules.prioritize_ilp_all.output.all_w_c,
        rules.prioritize_rwr.output.es,
        "analyses/zonation/priocomp/08_abf_es/08_abf_es_out/08_abf_es.rank.compressed.tif",
        rules.prioritize_ilp_es.output.es,
        rules.prioritize_rwr.output.es_c,
        "analyses/zonation/priocomp/10_abf_es_cst/10_abf_es_cst_out/10_abf_es_cst.rank.compressed.tif",
        rules.prioritize_ilp_es.output.es_c,
        rules.prioritize_rwr.output.bd,
        "analyses/zonation/priocomp/12_abf_bd/12_abf_bd_out/12_abf_bd.rank.compressed.tif",
        rules.prioritize_ilp_bd.output.bd,
        rules.prioritize_rwr.output.bd_c,
        "analyses/zonation/priocomp/14_abf_bd_cst/14_abf_bd_cst_out/14_abf_bd_cst.rank.compressed.tif",
        rules.prioritize_ilp_bd.output.bd_c
    output:
        "analyses/comparison/cross_jaccard.csv"
    log:
        "logs/compare_results_jaccard.log"
    message:
        "Comparing results overlap with Jaccard's index..."
    run:
        llogger = utils.get_local_logger("compare_jaccard", log[0])
        # Define thresholds for the top 10% and the low 10% for all rasters
        thresholds = [(0.0, 0.1, 0.0, 0.1), (0.9, 1.0, 0.9, 1.0)]
        jaccard_coefs = similarity.cross_jaccard(input, thresholds,
                                                 verbose=False, logger=llogger)
        llogger.info("Saving results to {}".format(output[0]))
        jaccard_coefs.to_csv(output[0], index=False)


rule compare_mcs:
    input:
        rules.postprocess_rwr.output.all_w,
        "analyses/zonation/priocomp/04_abf_all_wgt/04_abf_all_wgt_out/04_abf_all_wgt_nwout1.shp",
        rules.postprocess_ilp.output.all_w,
        rules.postprocess_rwr.output.all_w_c,
        "analyses/zonation/priocomp/06_abf_all_wgt_cst/06_abf_all_wgt_cst_out/06_abf_all_wgt_cst_nwout1.shp",
        rules.postprocess_ilp.output.all_w_c,
        rules.postprocess_rwr.output.es,
        "analyses/zonation/priocomp/08_abf_es/08_abf_es_out/08_abf_es_nwout1.shp",
        rules.postprocess_ilp.output.es,
        rules.postprocess_rwr.output.es_c,
        "analyses/zonation/priocomp/10_abf_es_cst/10_abf_es_cst_out/10_abf_es_cst_nwout1.shp",
        rules.postprocess_ilp.output.es_c,
        rules.postprocess_rwr.output.bd,
        "analyses/zonation/priocomp/12_abf_bd/12_abf_bd_out/12_abf_bd_nwout1.shp",
        rules.postprocess_ilp.output.bd,
        rules.postprocess_rwr.output.bd_c,
        "analyses/zonation/priocomp/14_abf_bd_cst/14_abf_bd_cst_out/14_abf_bd_cst_nwout1.shp",
        rules.postprocess_ilp.output.bd_c
    output:
        "analyses/comparison/cross_mcs.csv"
    log:
        "logs/compare_results_mcs.log"
    message:
        "Comparing NUTS2 mean ranks with MCS..."
    run:
        llogger = utils.get_local_logger("compare_mcs", log[0])
        value_fields = ['_mean', 'Men_rnk', '_mean', '_mean', 'Men_rnk', '_mean',
                        '_mean', 'Men_rnk', '_mean', '_mean', 'Men_rnk', '_mean',
                        '_mean', 'Men_rnk', '_mean', '_mean', 'Men_rnk', '_mean']
        mcs_scores_all = similarity.cross_mcs(input, value_fields,
                                              verbose=False, logger=llogger)
        llogger.info("Saving results to {}".format(output[0]))
        mcs_scores_all.to_csv(output[0], index=False)

rule compute_variation:
    input:
        no_costs=[rules.postprocess_rwr.output.all_w,
                  "analyses/zonation/priocomp/04_abf_all_wgt/04_abf_all_wgt_out/04_abf_all_wgt_nwout1.shp",
                  rules.postprocess_ilp.output.all_w,
                  rules.postprocess_rwr.output.es,
                  "analyses/zonation/priocomp/08_abf_es/08_abf_es_out/08_abf_es_nwout1.shp",
                  rules.postprocess_ilp.output.es,
                  rules.postprocess_rwr.output.bd,
                  "analyses/zonation/priocomp/12_abf_bd/12_abf_bd_out/12_abf_bd_nwout1.shp",
                  rules.postprocess_ilp.output.bd],
        costs=[rules.postprocess_rwr.output.all_w_c,
               "analyses/zonation/priocomp/06_abf_all_wgt_cst/06_abf_all_wgt_cst_out/06_abf_all_wgt_cst_nwout1.shp",
               rules.postprocess_ilp.output.all_w_c,
               rules.postprocess_rwr.output.es_c,
               "analyses/zonation/priocomp/10_abf_es_cst/10_abf_es_cst_out/10_abf_es_cst_nwout1.shp",
               rules.postprocess_ilp.output.es_c,
               rules.postprocess_rwr.output.bd_c,
               "analyses/zonation/priocomp/14_abf_bd_cst/14_abf_bd_cst_out/14_abf_bd_cst_nwout1.shp",
               rules.postprocess_ilp.output.bd_c]
    output:
        no_costs="analyses/comparison/nuts2_rank_variation.shp",
        costs="analyses/comparison/nuts2_rank_variation_costs.shp"
    log:
        no_costs="logs/compute_variation.log",
        costs="logs/compute_variation_costs.log"
    message:
        "Computing mean and variation stats of rank priorities for NUTS2 units..."
    run:
        # No costs
        llogger = utils.get_local_logger("compute_variation", log.no_costs)
        # FIXME: Codes are now hardcoded
        input_codes = ["rwr_all", "zon_all", "ilp_all",
                       "rwr_es", "zon_es", "ilp_es",
                       "rwr_bd", "zon_bd", "ilp_bd"]
        out_feature = similarity.plu_variation(input.no_costs, input_codes,
                                               logger=llogger)
        llogger.info("      Saving results to {}".format(output.no_costs))
        out_feature.to_file(output.no_costs)

        # Costs
        llogger = utils.get_local_logger("compute_variation_costs", log.costs)
        # FIXME: Codes are now hardcoded
        input_codes = ["rwr_all_costs", "zon_all_costs", "ilp_all_costs",
                       "rwr_es_costs", "zon_es_costs", "ilp_es_costs",
                       "rwr_bd_costs", "zon_bd_costs", "ilp_bd_costs"]
        out_feature = similarity.plu_variation(input.costs, input_codes,
                                               logger=llogger)
        llogger.info("      Saving results to {}".format(output.costs))
        out_feature.to_file(output.costs)

## Auxiliary operations -----------------------------------------------------

rule generate_range_data:
    input:
        "data/data_manifest.yml"
    output:
        csv="data/feature_ranges.csv"
    log:
        "logs/generate_range_data.log"
    message:
        "Generating feature ranges..."
    run:
        llogger = utils.get_local_logger("generate_range_data", log[0],
                                         debug=True)

        range_stats = pd.DataFrame({"feature": [], "count": [], "sum": [],
                                    "q25_ol": [], "mean_ol": [],
                                    "median_ol": [], "q75_ol": []})

        features = ES_DST_DATASETS + BD_DST_DATASETS
        for i, feature in enumerate(features):
            prefix = utils.get_iteration_prefix(i+1, len(features))

            llogger.info("{} Processing {}".format(prefix, feature))

            feature_stats = spatutils.get_range_size(feature, logger = llogger)
            range_stats = pd.concat([range_stats, feature_stats])

        llogger.info(" Saving results to {}".format(output.csv))
        range_stats.to_csv(output.csv, columns=["feature", "count", "sum",
                                                "q25_ol", "mean_ol",
                                                "median_ol", "q75_ol"],
                                                index=False)

rule generate_table_S1:
    input:
        "data/data_manifest.yml"
    output:
        md="reports/tables/02_table_S1.md",
        docx="reports/tables/02_table_S1.docx"
    log:
        "logs/generate_table_S1.log"
    message:
        "Generating table S1..."
    run:
        llogger = utils.get_local_logger("generate_table_S1", log[0])

        spp_table = dm.get_tabular(collection="european_tetrapods")
        with open(output.md, 'w') as outfile:
            outfile.write(spp_table)
        # You need to have pandox around for this to work
        output = pypandoc.convert_file(output.md, "docx",
                                       outputfile=output.docx)

## Tests -----------------------------------------------------

rule generate_test_data:
    input:
        features=UDR_SRC_DATASETS,
        clip_shp="tests/scratch/data/nl_clipper.shp"
    run:
        llogger = utils.get_local_logger("clip_test_data")
        nsteps = len(input.features)

        # Read in the clipping shapefile
        with fiona.open(input.clip_shp, "r") as shapefile:
            features = [feature["geometry"] for feature in shapefile]

            for i, s_raster in enumerate(input.features):
                prefix = utils.get_iteration_prefix(i+1, nsteps)

                llogger.info("{0} Clipping dataset {1}".format(prefix, s_raster))
                # Clip data. NOTE: UDR species rasters do not have a SRS defined,
                # but they are in EPSG:3035
                target_crs = rasterio.crs.CRS({'init':'epsg:3035'})
                with rasterio.open(s_raster) as src:
                    nodata_value = src.nodata
                    src_data = src.read(1, masked=True)

                    out_image, out_transform = rasterio.tools.mask.mask(src, features,
                                                                        crop=True)
                    # Replace NoData with 0s
                    out_image_src = out_image.copy()
                    out_image_src[out_image_src == nodata_value] = 0
                    if np.max(out_image_src) > 0:
                        out_meta = src.meta.copy()
                        out_meta.update({"driver": "GTiff",
                                         "height": out_image.shape[1],
                                         "width": out_image.shape[2],
                                         "transform": out_transform,
                                         "crs": target_crs})
                        # Check/create paths
                        input_path = os.path.dirname(s_raster)
                        input_name = os.path.basename(s_raster)
                        parent_dir = input_path.split(os.path.sep)[-1]
                        target_dir = os.path.join("tests/scratch/data/european_tetrapods_nl/",
                                                  parent_dir)
                        if not os.path.exists(target_dir):
                            os.makedirs(target_dir)
                        clipped_raster = os.path.join(target_dir, input_name)

                        llogger.debug("{0} Target dataset {1}".format(prefix, clipped_raster))
                        with rasterio.open(clipped_raster, "w", **out_meta) as dest:
                            dest.write(out_image)
                    else:
                        llogger.warning("{0} Empty dataset {1}".format(prefix, s_raster))

rule test_rwr:
    input:
        es=rules.generate_test_data.output
    output:
        es="tests/scratch/test_rwr_es.tif"
    message:
        "Testing RWR..."
    run:
        llogger = utils.get_local_logger("test_rwr_es")
        rwr.calculate_rwr(input.es, output.es, logger=llogger)

rule test_ilp:
    # Rule to test the result against the R implementation of Gurobi
    # maximum coverage problem
    input:
        spp_files=expand("/home/jlehtoma/dev/git-data/zonation-tutorial/data/species{ID}.tif",
                         ID=list(range(1, 8))),
        cost_file="/home/jlehtoma/dev/git-data/zonation-tutorial/data/cost.tif"
    output:
        "analyses/ILP/test_implementation/test_ilp_new_wcost.tif"
    log:
        "logs/test_ILP_new_wcost.log"
    message:
        "Runnning ILP tests..."
    run:
        llogger = utils.get_local_logger("test_optimize_gurobi", log[0])
        gurobi.prioritize_gurobi(input.spp_files, output[0], input.cost_file, logger=llogger,
                                 ol_normalize=True, save_intermediate=True,
                                 verbose=True)

rule test_ilp_hierarchy:
    # Rule to test if the solutions produced in the hierarchical ILP
    # optimization are indeed completely nested. In a set of
    # [f+x, ...,  f+nx] where f is 0,  n is a the number of consequtive
    # steps - 1, and x is 1 / n, make sure that the solution at f+(n-1) is
    # always a complete spatial subset of f+n.
    input:
        "analyses/ILP/ilp_es"
    output:
        "analyses/ILP/ilp_es/check.txt"
    log:
        "logs/test_ILP_hierarchy.log"
    message:
        "Testing ILP hierarchy..."
    run:
        llogger = utils.get_local_logger("test_ilp_hierarchy", log[0])
        # What's the file prefix?
        file_prefix = "budget_level_"
        # Define step used
        step = 0.05
        fractions = np.linspace(step, 1.0, 1/step)[:-1]
        # Genrate names for individual solutions
        solution_files = ["{0}{1}.tif".format(file_prefix, str(round(frac, 2)).replace(".", "_")) for frac in fractions]
        n_rasters = len(solution_files)

        previous_selected = None
        current_selected = None
        # Store names of non-subset rasters
        non_subsets = []

        for i, solution_file in enumerate(reversed(solution_files)):
            input_raster = os.path.join(input[0], solution_file)
            prefix = utils.get_iteration_prefix(i+1, n_rasters)
            llogger.info("{0} Processing raster {1}".format(prefix,
                                                            input_raster))
            with rasterio.open(input_raster) as in_raster:
                in_src = in_raster.read(1, masked=True)

                if i == 0:
                    # Store the selected elements. Original in_src has three
                    # values: [0, 1, --]. Fill masked array's masked values
                    # with 0's. Now 1s are selected unints, 0s either
                    # unselected units or NoData.
                    previous_selected = ma.filled(in_src, 0) == 1
                else:
                    # Get the current selected (see comment above)
                    current_selected = ma.filled(in_src, 0) == 1
                    # See if the current selected is a true subset of the
                    # previous
                    #import pdb; pdb.set_trace()
                    if not np.all(previous_selected[current_selected]):
                        llogger.warning("{0} is not a subset of {1}".format(solution_files[i],
                                                                            solution_files[i-1]))
                        non_subsets.append(solution_files[i])

                    previous_selected = current_selected

        with open(output[0], 'w') as file_handler:
            for item in non_subsets:
                file_handler.write("{}\n".format(item))

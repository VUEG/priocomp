import fiona
import gdal
import numpy as np
import os
import rasterio
from importlib.machinery import SourceFileLoader


rescale = SourceFileLoader("data_processing.rescale", "src/data_processing/rescale.py").load_module()
rwr = SourceFileLoader("rwr.calculate_rwr", "src/RWR/rwr.py").load_module()
utils = SourceFileLoader("src.utils", "src/utils.py").load_module()

## GLOBALS --------------------------------------------------------------------

# Analysis extent
PROJECT_EXTENT = {"bottom": 1000000.0, "left": 2000000.0, "right": 6526000.0,
                  "top": 5410000.0}
# EPSG for project
PROJECT_CRS = 3035

# Offset the bounds given in extent_yml. Values are in order
# (left, bottom, right, top) and interpreted in the CRS units. valuesq
# are added to bounds given by extent.yml
OFFSET = (100000, 100000, 0, 0)
# Which eurostat countries are included in the processed output? The
# following countries have been removed:
# "CY", "CH", "IS", "HR", "NO", "ME", "MT", "MK", "TR"
PROJECT_COUNTRIES = ["AT", "BE", "BG", "CZ", "DE", "DK", "ES", "EL", "EE",
                     "FR", "FI", "IT", "HU", "IE", "NL", "LU", "LI", "LT",
                     "LV", "PL", "SE", "RO", "PT", "SK", "SI", "UK"]

# Get datasets from the data manifest file.
data_manifest = utils.parse_data_manifest("data/data_manifest.yml")

external_data = "data/external"
beehub_url = "https://beehub.nl/environmental-geography-group"

# Define source and desination datasets. NOTE: data/Snakefile must be run
# before this Snakefile will work
DATADRYAD_SRC_DATASETS = [path for path in data_manifest['datadryad'] if path.endswith(".tif")]
DATADRYAD_SRC_DATASETS = [url.replace(beehub_url, external_data) for url in DATADRYAD_SRC_DATASETS]

EUROSTAT_SRC_DATASETS = [url.replace(beehub_url, external_data) for url in data_manifest['eurostat']]
# Define NUTS data separately. Shapefile constitute of several actual files on
# file system, so dealing with them properly is a bit tricky.
NUTS_LEVEL0_DATA = [path for path in EUROSTAT_SRC_DATASETS if "level0" in path and not path.endswith(".txt")]
NUTS_LEVEL2_DATA = [path for path in EUROSTAT_SRC_DATASETS if "level2" in path and not path.endswith(".txt")]

PROVIDE_SRC_DATASETS = [path for path in data_manifest['provide'] if path.endswith(".tif")]
PROVIDE_SRC_DATASETS = [url.replace(beehub_url, external_data) for url in PROVIDE_SRC_DATASETS]

# UDR collection "european_tetrapods" is already correctly formatted, place
# it directly to "processed/features"
UDR_SRC_DATASETS = [url.replace(beehub_url, "data/processed/features") for url in data_manifest["udr"]]

ALL_SRC_DATASETS = DATADRYAD_SRC_DATASETS + PROVIDE_SRC_DATASETS + UDR_SRC_DATASETS

# PROJECT RULES ----------------------------------------------------------------

rule all:
    input:
        "analyses/RWR/eu26_rwr.tif"

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
                # Insert a new field
                meta['schema']['properties']['ID'] = 'int'

                ID = 1
                with fiona.open(enhanced_shp, 'w', **meta) as sink:
                    # Loop over features
                    for f in source:
                        f['properties']['ID'] = ID
                        ID += 1
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
        enhanced=temp([path.replace("external", "interim/enhanced") for path in NUTS_LEVEL2_DATA]),
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
        # Reproject to EPSG:3035 from EPSG:4258
        input_shp = utils.pick_from_list(input.shp, ".shp")
        reprojected_shp = utils.pick_from_list(output.reprojected, ".shp")
        cmd_str = 'ogr2ogr {0} -t_srs "EPSG:{1}" {2}'.format(reprojected_shp, PROJECT_CRS, input_shp)
        shell(cmd_str)
        llogger.debug("Reprojected NUTS level 2 data from EPSG:4258 to EPSG:3035")
        llogger.debug(cmd_str)

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
        llogger.debug("Clipped NUTS level 2 data to analysis bounds: {}".format(bounds))
        llogger.debug("Selected only a subset of eurostat countries:")
        llogger.debug(" " + ", ".join(PROJECT_COUNTRIES))
        llogger.debug("Resulting file: {}".format(processed_shp))

rule rasterize_nuts_level0_data:
    input:
        rules.preprocess_nuts_level0_data.output.processed
    output:
        utils.pick_from_list(rules.preprocess_nuts_level0_data.output.processed, ".shp").replace(".shp", ".tif")
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
        # Rasterize
        cmd_str = "gdal_rasterize -l {} ".format(layer_shp) + \
                  "-a ID -tr 1000 1000 -te {} ".format(bounds) + \
                  "-ot Int16 -a_nodata -32768 -co COMPRESS=DEFLATE " + \
                  "{0} {1}".format(input_shp, output[0])
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

rule harmonize_data:
    input:
        external=DATADRYAD_SRC_DATASETS+PROVIDE_SRC_DATASETS,
        like_raster=[path for path in DATADRYAD_SRC_DATASETS if "woodprod_average" in path][0],
        clip_shp=utils.pick_from_list(rules.preprocess_nuts_level0_data.output.processed, ".shp")
    output:
        # NOTE: UDR_SRC_DATASETS do not need to processed
        warped=temp([path.replace("external", "interim/warped") for path in DATADRYAD_SRC_DATASETS+PROVIDE_SRC_DATASETS]),
        harmonized=[path.replace("external", "processed/features") for path in DATADRYAD_SRC_DATASETS+PROVIDE_SRC_DATASETS]
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
            prefix = utils.get_iteration_prexix(i+1, nsteps)
            if s_raster == input.like_raster:
                llogger.info("{0} Copying dataset {1}".format(prefix, s_raster))
                llogger.debug("{0} Target dataset {1}".format(prefix, warped_raster))
                ret = shell("cp {s_raster} {warped_raster}", read=True)
            else:
                llogger.info("{0} Warping dataset {1}".format(prefix, s_raster))
                llogger.debug("{0} Target dataset {1}".format(prefix, warped_raster))
                ret = shell("rio warp " + s_raster + " --like " + input.like_raster + \
                            " " + warped_raster + " --dst-crs " + str(PROJECT_CRS) + \
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

rule ol_normalize_data:
    input:
        rules.harmonize_data.output.harmonized#+UDR_SRC_DATASETS
    output:
        [path.replace("processed/features", "processed/features_ol_normalized") for path in rules.harmonize_data.output.harmonized]
    log:
        "logs/ol_normalize_data.log"
    message:
        "Normalizing data based on occurrence levels..."
    run:
        llogger = utils.get_local_logger("ol_normalize_data", log[0])
        for i, s_raster in enumerate(input):
            prefix = utils.get_iteration_prexix(i+1, len(input))
            llogger.info("{0} (OL) Normalizing dataset {1}".format(prefix, s_raster))
            # NOTE: looping over input and output only works if they have
            # exactly the same definition. Otherwise order may vary.
            rescale.rescale_raster(input[i], output[i], method="ol_normalize",
                                   logger=llogger)

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
            rescale.rescale_raster(input[i], output[i], method="normalize",
                                   verbose=False)

## Set up and run analyses -----------------------------------------------------

# RWR --------------------------------------------------------------------------

rule calculate_rwr:
    input:
        rules.harmonize_data.output.harmonized+UDR_SRC_DATASETS
    output:
        "analyses/RWR/eu26_rwr.tif"
    message:
        "Calculating RWR..."
    run:
        rwr.calculate_rwr(input, output[0], verbose=True)


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

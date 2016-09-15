import fiona
import gdal
import numpy as np
import os
import pdb
import rasterio
import rasterstats
from importlib.machinery import SourceFileLoader

## LOAD MODULES --------------------------------------------------------------
utils = SourceFileLoader("lib.utils", "src/00_lib/utils.py").load_module()
spatutils = SourceFileLoader("lib.spatutils", "src/00_lib/spatutils.py").load_module()
gurobi = SourceFileLoader("analysis.gurobi", "src/02_analysis/gurobi.py").load_module()
rwr = SourceFileLoader("analysis.rwr", "src/02_analysis/rwr.py").load_module()
similarity = SourceFileLoader("results.similarity", "src/03_post_processing/raster_similarity.py").load_module()


## GLOBALS --------------------------------------------------------------------

# Analysis extent
PROJECT_EXTENT = {"bottom": 1000000.0, "left": 2000000.0, "right": 6526000.0,
                  "top": 5410000.0}
# EPSG for project
PROJECT_CRS = 3035

# Offset the bounds given in extent_yml. Values are in order
# (left, bottom, right, top) and interpreted in the CRS units. values
# are added to bounds given by PROJECT_EXTENT
OFFSET = (100000, 100000, 0, 0)
# Which eurostat countries are included in the processed output? The
# following countries have been removed:
# "CY", "CH", "IS", "HR", "NO", "ME", "MT", "MK", "TR"
PROJECT_COUNTRIES = ["AT", "BE", "BG", "CZ", "DE", "DK", "ES", "EL", "EE",
                     "FR", "FI", "IT", "HU", "IE", "NL", "LU", "LI", "LT",
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

# Get specific NUTS collections from Eurostat
NUTS_LEVEL0_DATA = [url.replace(beehub_url, external_data) for url in dm.get_resources(collection="nuts_level0", full_path=True)]
NUTS_LEVEL2_DATA = [url.replace(beehub_url, external_data) for url in dm.get_resources(collection="nuts_level2", full_path=True)]

PROVIDE_SRC_DATASETS = [url.replace(beehub_url, external_data) for url in dm.get_resources(provider="provide", full_path=True)]

# UDR collection "european_tetrapods" is already correctly formatted, place
# it directly to "processed/features"
UDR_SRC_DATASETS = [url.replace(beehub_url, feature_data) for url in dm.get_resources(provider="udr", full_path=True)]

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

            # Finally, rescale (normalize) dataset if needed
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

rule calculate_rwr:
    input:
        all=rules.harmonize_data.output.harmonized+UDR_SRC_DATASETS,
        es=rules.harmonize_data.output.harmonized,
        bd=UDR_SRC_DATASETS
    output:
        all="analyses/RWR/rwr_eu26_all.tif",
        all_w="analyses/RWR/rwr_eu26_all_weights.tif",
        es="analyses/RWR/rwr_eu26_es.tif",
        bd="analyses/RWR/rwr_eu26_bd.tif"
    log:
        all="logs/calculate_rwr_eu26_all.log",
        all_w="logs/calculate_rwr_eu26_all_weights.log",
        es="logs/calculate_rwr_eu26_es.log",
        bd="logs/calculate_rwr_eu26_bd.log"
    message:
        "Calculating RWR..."
    run:
        # Without weights
        llogger = utils.get_local_logger("calculate_rwr_all", log.all)
        rwr.calculate_rwr(input.all, output.all, logger=llogger)
        # With weights
        llogger = utils.get_local_logger("calculate_rwr_all_weights", log.all_w)
        rwr.calculate_rwr(input.all, output.all_w, weights=WEIGHTS,
                          logger=llogger)

        llogger = utils.get_local_logger("calculate_rwr_es", log.es)
        rwr.calculate_rwr(input.es, output.es, logger=llogger)
        llogger = utils.get_local_logger("calculate_rwr_bd", log.bd)
        rwr.calculate_rwr(input.bd, output.bd, logger=llogger)

rule postprocess_rwr:
    input:
        all=rules.calculate_rwr.output.all,
        all_w=rules.calculate_rwr.output.all_w,
        es=rules.calculate_rwr.output.es,
        bd=rules.calculate_rwr.output.bd,
        plu=utils.pick_from_list(rules.preprocess_nuts_level0_data.output.processed,
                                 ".shp")
    output:
        all="analyses/RWR/rwr_eu26_all_stats.geojson",
        all_w="analyses/RWR/rwr_eu26_all_weights_stats.geojson",
        es="analyses/RWR/rwr_eu26_es_stats.geojson",
        bd="analyses/RWR/rwr_eu26_bd_stats.geojson"
    log:
        all="logs/postprocess_rwr_eu26_all.log",
        all_w="logs/postprocess_rwr_eu26_all_weights.log",
        es="logs/postprocess_rwr_eu26_es.log",
        bd="logs/postprocess_rwr_eu26_bd.log"
    message:
        "Post-processing RWR..."
    run:
        llogger = utils.get_local_logger("calculate_rwr_all", log.all)
        llogger.info(" [1/4] Post-processing {}".format(input.all))
        shell("fio cat {input.plu} | rio zonalstats -r {input.all} > {output.all}")
        llogger = utils.get_local_logger("calculate_rwr_all_weights", log.all_w)
        llogger.info(" [2/4] Post-processing {}".format(input.all_w))
        shell("fio cat {input.plu} | rio zonalstats -r {input.all_w} > {output.all_w}")
        llogger = utils.get_local_logger("calculate_rwr_es", log.es)
        llogger.info(" [3/4] Post-processing {}".format(input.es))
        shell("fio cat {input.plu} | rio zonalstats -r {input.es} > {output.es}")
        llogger = utils.get_local_logger("calculate_rwr_bd", log.bd)
        llogger.info(" [4/4] Post-processing {}".format(input.bd))
        shell("fio cat {input.plu} | rio zonalstats -r {input.bd} > {output.bd}")





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

# Gurobi ----------------------------------------------------------------------

rule prioritize_gurobi:
    input:
        all=rules.harmonize_data.output.harmonized+UDR_SRC_DATASETS,
        es=rules.harmonize_data.output.harmonized,
        bd=UDR_SRC_DATASETS
    output:
        #all="analyses/ILP/ilp_eu26_all.tif",
        all_w="analyses/ILP/ilp_eu26_all_weights.tif"
        #es="analyses/ILP/ilp_eu26_es.tif",
        #bd="analyses/ILP/ilp_eu26_bd.tif"
    log:
        all="logs/prioritize_ilp_eu26_all.log",
        all_w="logs/prioritize_ilp_eu26_all_weights.log",
        es="logs/prioritize_ilp_eu26_es.log",
        bd="logs/prioritize_ilp_eu26_bd.log"
    message:
        "Optimizing with Gurobi..."
    run:
        # Without weights
        #llogger = utils.get_local_logger("optimize_gurobi_all", log.all)
        #gurobi.prioritize_gurobi(input.all, output.all, logger=llogger,
        #                         ol_normalize=True, save_intermediate=True,
        #                         verbose=True)
        # With weights
        llogger = utils.get_local_logger("optimize_gurobi_all_weights",
                                         log.all_w)
        gurobi.prioritize_gurobi(input.all, output.all_w, logger=llogger,
                                 ol_normalize=True, weights=WEIGHTS,
                                 save_intermediate=False, verbose=True)

        #llogger = utils.get_local_logger("optimize_gurobi_es", log.es)
        #gurobi.prioritize_gurobi(input.es, output.es, logger=llogger,
        #                         ol_normalize=True, save_intermediate=False,
        #                         verbose=True)
        #llogger = utils.get_local_logger("optimize_gurobi_bd", log.bd)
        #gurobi.prioritize_gurobi(input.bd, output.bd, logger=llogger,
        #                         ol_normalize=True, save_intermediate=False,
        #                         verbose=True)

## Compare results ------------------------------------------------------------

rule compare_correlation:
    input:
        rwr=rules.calculate_rwr.output.all_w,
        zon="analyses/zonation/priocomp/04_abf_wgt/04_abf_wgt_out/04_abf_wgt.rank.compressed.tif",
        ilp=rules.prioritize_gurobi.output.all_w
    output:
        "analyses/comparison/cross_correlation.csv"
    log:
        all="logs/compare_results_correlation.log"
    message:
        "Comparing results correlation using Kendal tau..."
    run:
        llogger = utils.get_local_logger("compare_correlation", log[0])
        correlations = similarity.cross_correlation(input, verbose=False,
                                                    logger=llogger)
        llogger.info("Saving results to {}".format(output[0]))
        correlations.to_csv(output[0], index=False)

rule compare_jaccard:
    input:
        rwr=rules.calculate_rwr.output.all_w,
        zon="analyses/zonation/priocomp/04_abf_wgt/04_abf_wgt_out/04_abf_wgt.rank.compressed.tif",
        ilp=rules.prioritize_gurobi.output.all_w
    output:
        "analyses/comparison/cross_jaccard.csv"
    log:
        all="logs/compare_results_jaccard.log"
    message:
        "Comparing results overlap with Jaccard's index..."
    run:
        llogger = utils.get_local_logger("compare_jaccard", log[0])
        # Define thresholds every 10% of the rank priority map
        thresholds = np.arange(0.05, 1.0, 0.05)
        jaccard_coefs = similarity.cross_jaccard(input, thresholds,
                                                 verbose=False, logger=llogger)
        llogger.info("Saving results to {}".format(output[0]))
        jaccard_coefs.to_csv(output[0], index=False)

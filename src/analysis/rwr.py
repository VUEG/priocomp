#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Functions and utilities for calculating the rarity-weighted richness score.

Module can be used alone or as part of Snakemake workflow.
"""

import click
import glob
import logging
import numpy as np
import numpy.ma as ma
import pdb
import rasterio
import time

from importlib.machinery import SourceFileLoader
from scipy.stats.mstats import rankdata
from timeit import default_timer as timer

spatutils = SourceFileLoader("data_processing.spatutils",
                             "src/data_processing/spatutils.py").load_module()
utils = SourceFileLoader("src.utils", "src/utils.py").load_module()


def calculate_rwr(input_rasters, output_raster, compress='DEFLATE',
                  verbose=False, logger=None):
    """ Calculate rarity-weighted richness (RWR) based on input rasters.

    Function does the following steps:
    1. Read in data from a raster and occurrence level normalize the data
       using ol_normalize().
    2. Add the the OL normalized data to a dynamically increasing array (i.e.
       sum all the values).
    3. Rank elements of the array.
    4. Rescale ranks in range [0, 1].

    Data is read in as ma.MaskedArrays where the mask represents raster NoData.
    NoData is not propagated, i.e. value + NoData is always value. While bounds
    and resolution are checked, summation will fail if all the input rasters do
    not have the same extent and resolution (and hence implicitly CRS).

    NOTE: 0 is used as a NoData value internally. Having acutal 0s as values in
    the input data may cause unintended consequences.

    Ranking of the values in the summed array is done using "average" strategy,
    i.e. the average of the ranks that would have been assigned to all the tied
    values is assigned to each value.

    :param input_rasters: List of String paths of input rasters.
    :param output_raster: String path to raster file to be created.
    :param compress: String compression level used for the output raster.
    :param verbose Boolean indicating how much information is printed out.
    :param log_file String path to log file used.
    """
    # 1. Setup  --------------------------------------------------------------

    all_start = timer()
    load_start = timer()

    if not logger:
        logging.basicConfig()
        llogger = logging.getLogger('calculate_rwr')
        llogger.setLevel(logging.DEBUG if verbose else logging.INFO)
    else:
        llogger = logger

    # Check inputs
    assert len(input_rasters) > 0, "Input rasters list cannot be empty"
    assert len(output_raster) != "", "Output raster path cannot be empty"

    # 2. Sum value  -----------------------------------------------------

    llogger.info(" [** SUMMING RASTERS **]")

    # Create a sum array, i.e. sum all (occurrence level normalized) raster
    # values in input_rasters together.
    sum_array = spatutils.sum_raster(input_rasters, olnormalize=True,
                                     logger=llogger)

    load_end = timer()
    load_elapsed = round(load_end - load_start, 2)
    llogger.info(" [TIME] Pre-processing took {} sec".format(load_elapsed))

    # 3. Rank RWR data --------------------------------------------------------

    post_start = timer()
    llogger.info(" [** POST-PROCESSING **]")
    llogger.info(" [1/3] Ranking values (this can take a while...)")

    nodata_value = -3.4e+38

    # Use 0s from summation as a mask
    sum_array = ma.masked_values(sum_array, 0.0)
    rank_array = rankdata(sum_array)
    # rankdata() will mark all masked values with 0, make those NoData again.
    np.place(rank_array, rank_array == 0, nodata_value)
    # Create a masked array. Convert remaining zeros to NoData and use that as
    # a mask
    # NOTE: use ma.masked_values() because we're replacing with a float value
    rank_array = ma.masked_values(rank_array, nodata_value)

    # 4. Recale data into range [0, 1] ----------------------------------------
    llogger.info(" [2/3] Rescaling ranks")
    rank_array = spatutils.normalize(rank_array)
    rank_array = rank_array.astype(np.float32)

    # 5. Write out the data
    llogger.info(" [3/3] Writing output to {}".format(output_raster))

    # Get the raster profile from the input raster files
    profile = spatutils.get_profile(input_rasters, logger=llogger)

    # Rescaled data is always float32, and we have only 1 band. Remember
    # to set NoData-value correctly.
    profile.update(dtype=rasterio.float32, compress=compress,
                   nodata=-3.4e+38)

    with rasterio.open(output_raster, 'w', **profile) as dst:
        dst.write_mask(ma.getmask(rank_array))
        dst.write(rank_array, 1)

    post_end = timer()
    post_elapsed = round(post_end - post_start, 2)
    llogger.info(" [TIME] Post-processing took {} sec".format(post_elapsed))

    all_end = timer()
    all_elapsed = round(all_end - all_start, 2)
    llogger.info(" [TIME] All processing took {} sec".format(all_elapsed))


@click.command()
@click.option('-v', '--verbose', is_flag=True)
@click.option('-e', '--extension', nargs=1, default=".tif",
              help='Raster file extension')
@click.argument('indir', nargs=1, type=click.Path(exists=True))
@click.argument('outfile', nargs=1)
def cli(indir, outfile, extension, verbose):
    """ Command-line interface."""
    start_time = time.time()
    # List files in input directorys
    file_iterator = glob.iglob('{0}/**/*{1}'.format(indir, extension),
                               recursive=True)
    input_rasters = [item for item in file_iterator]
    if verbose:
        click.echo(click.style(" Found {} rasters".format(len(input_rasters)),
                               fg='green'))
    calculate_rwr(input_rasters, outfile, verbose=verbose)
    exc_time = time.time() - start_time
    click.echo(click.style('Done in {} seconds!'.format(exc_time), fg='green'))

if __name__ == '__main__':
    cli()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Functions and utilities for calculating the rarity-weighted richness score.

Module can be used alone or as part of Snakemake workflow.
"""

import click
import glob
import numpy as np
import numpy.ma as ma
import os
import rasterio
import time

from importlib.machinery import SourceFileLoader
from scipy.stats.mstats import rankdata

rescale = SourceFileLoader("data_processing.rescale",
                           "src/data_processing/rescale.py").load_module()
utils = SourceFileLoader("src.utils", "src/utils.py").load_module()


def calculate_rwr(input_rasters, output_raster, compress='DEFLATE',
                  verbose=False, log_file=None):
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
    """
    assert len(input_rasters) > 0, "Input rasters list cannot be empty"

    # Set up logging
    llogger = utils.get_local_logger(calculate_rwr.__name__, log_file,
                                     debug=verbose)

    # Initiate the data structure for holding the summed values.
    sum_array = None
    # Placeholder for raster metadata
    profile = None

    for i, input_raster in enumerate(input_rasters):
        no_raster = i + 1
        n_rasters = len(input_rasters)

        if not os.path.exists(input_raster):
            raise OSError("Input raster {} not found".format(input_raster))

        with rasterio.open(input_raster) as in_src:
            llogger.info(" [{0}/{1}] Processing".format(no_raster, n_rasters) +
                         " raster {0}".format(input_raster))
            llogger.debug(" [{0}/{1} step 1] Reading".format(no_raster,
                                                             n_rasters) +
                          " in data and OL normalizing")

            # Read the first band, use zeros for NoData
            src_data = in_src.read(1, masked=True)
            src_data = ma.filled(src_data, 0)

            # 1. Occurrence level normalize data ------------------------------
            src_data = rescale.ol_normalize(src_data)

            # 2. Sum OL normalized data ---------------------------------------
            llogger.debug(" [{0}/{1} step 2] Summing values".format(no_raster,
                                                                    n_rasters))
            # If this is the first raster, use its dimensions to build an array
            # that holds the summed values.
            if i == 0:
                # Set up an array of zeros that has the correct dimensions
                sum_array = np.zeros_like(src_data, dtype=np.float32)
                # Also store the necessarry information for the output raster
                # Write the output product
                profile = in_src.profile
            sum_array += src_data

    # 3. Rank RWR data --------------------------------------------------------
    llogger.info(" Ranking values (this can take a while...)")

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
    llogger.debug(" Rescaling ranks")
    rank_array = rescale.normalize(rank_array)

    rank_array = rank_array.astype(np.float32)

    # 5. Write out the data
    llogger.info(" Writing output to {}".format(output_raster))

    # Rescaled data is always float32, and we have only 1 band. Remember
    # to set NoData-value correctly.
    profile.update(dtype=rasterio.float32, compress=compress,
                   nodata=-3.4e+38)

    with rasterio.open(output_raster, 'w', **profile) as dst:
        dst.write_mask(ma.getmask(rank_array))
        dst.write(rank_array, 1)


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

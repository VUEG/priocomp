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
import os
import rasterio
import time

from gurobipy import *
from importlib.machinery import SourceFileLoader

rescale = SourceFileLoader("data_processing.rescale",
                           "../data_processing/rescale.py").load_module()
utils = SourceFileLoader("src.utils", "../utils.py").load_module()


def optimize_gurobi(input_rasters, output_raster, budget,
                    compress='DEFLATE', verbose=False, logger=None):
    """ Solve maximum coverage problem for a set of input rasters and
    target budget using Gurobi solver.

    :param input_rasters: List of String paths of input rasters.
    :param output_raster: String path to raster file to be created.
    :param budget: numeric indicating the budget for the maximum coverage
                   problem.
    :param compress: String compression level used for the output raster.
    :param verbose Boolean indicating how much information is printed out.
    :param logger logger object to be used.
    """
    assert len(input_rasters) > 0, "Input rasters list cannot be empty"

    # Set up logging
    if not logger:
        logging.basicConfig()
        llogger = logging.getLogger('optimize_gurobi')
        llogger.setLevel(logging.DEBUG if verbose else logging.INFO)
    else:
        llogger = logger

    # Initiate the data structure for holding the summed values.
    sum_array = None
    # Placeholder for raster metadata
    profile = None
    n_rasters = len(input_rasters)

    for i, input_raster in enumerate(input_rasters):
        no_raster = i + 1

        if not os.path.exists(input_raster):
            raise OSError("Input raster {} not found".format(input_raster))

        with rasterio.open(input_raster) as in_src:
            prefix = utils.get_iteration_prexix(no_raster, n_rasters)
            llogger.info("{0} Processing raster {1}".format(prefix,
                                                            input_raster))
            llogger.debug("{0} Reading in data and OL ".format(prefix) +
                          "normalizing")

            # Read the first band, use zeros for NoData
            src_data = in_src.read(1, masked=True)
            src_data = ma.filled(src_data, 0)

            # 1. Occurrence level normalize data ------------------------------
            src_data = rescale.ol_normalize(src_data)

            # 2. Sum OL normalized data ---------------------------------------
            llogger.debug("{0} Summing values".format(prefix))
            # If this is the first raster, use its dimensions to build an array
            # that holds the summed values.
            if i == 0:
                # Set up an array of zeros that has the correct dimensions
                sum_array = np.zeros_like(src_data, dtype=np.float32)
                # Also store the necessarry information for the output raster
                # Write the output product
                profile = in_src.profile
            sum_array += src_data

    # Flatten the sum array
    sum_array = sum_array.flatten()
    # Create equal cost array
    cost = np.ones(sum_array.size)
    budget = budget * cost.size

    # 3. Optimize  --------------------------------------------------------
    llogger.info(" Optimizing...")

    try:
        # Create a new model
        m = Model("mip1")

        # Create variables
        for x in np.nditer(sum_array):
            m.addVar(vtype=GRB.BINARY)

        # Integrate new variables
        m.update()
        vars = m.getVars()

        # Set objective and constraint
        obj = LinExpr()
        expr = LinExpr()
        for i in range(sum_array.size):
            obj.addTerms(sum_array[i], vars[i])
            expr.addTerms(cost[i], vars[i])

        m.setObjective(obj, sense=GRB.MAXIMIZE)
        print(type(expr))
        m.addConstr(lhs=expr, sense=GRB.LESS_EQUAL, rhs=budget)

        m.update()

        #m.params.timeLimit = 100.0
        #m.params.mipgap = 0.001
        #m.params.solutionLimit = 1
        #m.params.presolve = -1

        m.optimize()

        print('Obj: %g' % m.objVal)

    except GurobiError:
        raise

    nodata_value = -3.4e+38

    # Use 0s from summation as a mask
    #sum_array = ma.masked_values(sum_array, 0.0)

    #np.place(rank_array, rank_array == 0, nodata_value)
    # Create a masked array. Convert remaining zeros to NoData and use that as
    # a mask
    # NOTE: use ma.masked_values() because we're replacing with a float value
    #rank_array = ma.masked_values(rank_array, nodata_value)

    # 4. Write out the data
    #llogger.info(" Writing output to {}".format(output_raster))

    # Rescaled data is always float32, and we have only 1 band. Remember
    # to set NoData-value correctly.
    #profile.update(dtype=rasterio.float32, compress=compress,
    #               nodata=-3.4e+38)

    #with rasterio.open(output_raster, 'w', **profile) as dst:
    #    dst.write_mask(ma.getmask(rank_array))
    #    dst.write(rank_array, 1)


@click.command()
@click.option('-v', '--verbose', is_flag=True)
@click.option('-e', '--extension', nargs=1, default=".tif",
              help='Raster file extension')
@click.option('-b', '--budget', required=True,
              help='Budget for the maximum coverage problem')
@click.argument('indir', nargs=1, type=click.Path(exists=True))
@click.argument('outfile', nargs=1)
def cli(indir, outfile, budget, extension, verbose):
    """ Command-line interface."""
    start_time = time.time()
    # List files in input directorys
    file_iterator = glob.iglob('{0}/**/*{1}'.format(indir, extension),
                               recursive=True)
    input_rasters = [item for item in file_iterator]
    if verbose:
        click.echo(click.style(" Found {} rasters".format(len(input_rasters)),
                               fg='green'))
    optimize_gurobi(input_rasters, outfile, budget, verbose=verbose)
    exc_time = time.time() - start_time
    click.echo(click.style('Done in {} seconds!'.format(exc_time), fg='green'))

if __name__ == '__main__':
    cli()

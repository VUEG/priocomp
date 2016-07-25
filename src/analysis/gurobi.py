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
import ipdb
import rasterio
import sys

from gurobipy import *
from importlib.machinery import SourceFileLoader
from scipy.stats.mstats import rankdata
from timeit import default_timer as timer

spatutils = SourceFileLoader("data_processing.spatutils",
                           "src/data_processing/spatutils.py").load_module()
utils = SourceFileLoader("src.utils", "src/utils.py").load_module()


def optimize_maxcover(x, budget, rij, verbose=False, logger=None):
    """ Solve maximum coverage problem using Gurobi.

    :param x: Numpy ndarray of planning unit costs. In this case the spatial
              representation of the planning units is not provided.
    :param budget: numeric indicating the budget for the maximum coverage
                   problem.
    :param rij: Numpy ndarray; Column sums of the matrix of representation
                levels of conservation features (rows) within planning units
                (columns).
    :param verbose: Boolean indicating how much information is printed out.
    :param logger: logger object to be used.
    :return
    """
    try:
        # Create a new model
        m = Model("mip1")

        # Create variables
        for element in np.nditer(rij):
            m.addVar(vtype=GRB.BINARY)

        # Integrate new variables
        m.update()
        vars = m.getVars()

        # Set objective and constraint
        obj = LinExpr()
        expr = LinExpr()

        for i in range(rij.size):
            obj.addTerms(rij[i], vars[i])
            expr.addTerms(x[i], vars[i])

        m.setObjective(obj, sense=GRB.MAXIMIZE)
        m.addConstr(lhs=expr, sense=GRB.LESS_EQUAL, rhs=budget)
        m.update()

        # The the log file directory from the logger, but log into a separate
        # file. NOTE: handler location is hard coded.
        log_dir = os.path.dirname(logger.handlers[1].baseFilename)
        log_file = os.path.join(log_dir, "gurobi.log")
        if os.path.exists(log_dir):
            logger.debug("See {} for Gurobi log".format(log_file))
            m.params.logToConsole = 0
            m.params.logFile = log_file
        else:
            logger.debug("Log dir {} not found".format(log_dir))
            m.params.logToConsole = 0
        # m.params.mipgap = 0.001
        # m.params.solutionLimit = 1
        # m.params.presolve = -1

        m.optimize()

    except GurobiError:
        raise

    # Construct a result array and return that
    res = np.asarray([var.x for var in m.getVars()], dtype=np.uint8)
    return res


def prioritize_gurobi(input_rasters, output_rank_raster, step=0.05,
                      save_intermediate=False, compress='DEFLATE',
                      ol_normalize=False, verbose=False, logger=None):
    """ Solve (multiple) maximum coverage problems for a set of input rasters.

    Create a hierarchical spatial prioritization using Gurobi solver to solve
    multiple optimization problems with different budget levels. Each budget
    level must be in range [0.0, 1.0] and corresponds to the area of interest.
    In other words, value of 0.1 would correspond to best 10% of the landscape.
    Budget levels will be automatically created using step value defined by the
    step arguments. Gurobi solver will then solve each problem using the
    representation levels in input_rasters. All the (binary) results are then
    summed together forming a selection frequency. Finally, the selection
    frequency is rescaled into range [0, 1] forming a rank priority raster.

    :param input_rasters: List of String paths of input rasters.
    :param output_raster: String path to the rank raster file to be created.
    :param step: numeric value in (0, 1) defining the step length for budget
                 levels.
    :param compress: String compression level used for the output raster.
    :param compress: Boolean setting OL (Occurrence Level) normalization.
    :param verbose Boolean indicating how much information is printed out.
    :param logger logger object to be used.
    """
    # Set up logging

    all_start = timer()
    load_start = timer()

    if not logger:
        logging.basicConfig()
        llogger = logging.getLogger('optimize_gurobi')
        llogger.setLevel(logging.DEBUG if verbose else logging.INFO)
    else:
        llogger = logger

    if len(input_rasters) < 1:
        llogger.error("Input rasters list cannot be empty")
        sys.exit(1)

    # Check inputs
    try:
        step = float(step)
    except ValueError:
        llogger.error("'step' must be coercible to float")
        sys.exit(1)
    assert step > 0.0 and step < 1.0, "Step argument must be in range (0, 1)"
    # Construct budget levels based on the step provided. 0.0 (nothing) and
    # 1.0 (everyhting) are not needed.
    budget_levels = np.linspace(0.0+step, 1.0, 1/step)[:-1]

    # Initiate the data structure for holding the summed values.
    sum_array = None
    # Placeholder for raster metadata
    profile = None
    n_rasters = len(input_rasters)

    llogger.info(" [** PRE-PROCESSING **]")

    if ol_normalize:
        llogger.info(" [NOTE] Using occurrence level normalization")

    for i, input_raster in enumerate(input_rasters):
        no_raster = i + 1

        if not os.path.exists(input_raster):
            raise OSError("Input raster {} not found".format(input_raster))

        with rasterio.open(input_raster) as in_src:
            prefix = utils.get_iteration_prexix(no_raster, n_rasters)
            llogger.info("{0} Processing raster {1}".format(prefix,
                                                            input_raster))
            llogger.debug("{0} Reading in data".format(prefix))

            # Read the first band, use zeros for NoData
            src_data = in_src.read(1, masked=True)
            src_data = ma.filled(src_data, 0)

            if ol_normalize:
                # 1. Occurrence level normalize data --------------------------
                llogger.debug("{0} OL normalizing data".format(prefix))
                src_data = spatutils.ol_normalize(src_data)

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
    # Get rid of zeros for now, but first get the mask
    mask = ma.getmask(ma.masked_values(sum_array, 0.0))
    sum_array = sum_array[~mask]
    mask = mask.reshape((profile["height"], profile["width"]))
    # Create equal cost array
    cost = np.ones(sum_array.size)

    load_end = timer()
    load_elapsed = round(load_end - load_start, 2)
    llogger.info(" [TIME] Pre-processing took {} sec".format(load_elapsed))

    # 3. Optimize  -----------------------------------------------------------

    opt_start = timer()
    llogger.info(" [** OPTIMIZING **]")
    blevels_str = ", ".join([str(level) for level in budget_levels])
    llogger.info(" [NOTE] Target budget levels: {}".format(blevels_str))

    # Construct a ndarray (matrix) that will hold the selection frequency
    sel_freq = np.zeros_like(src_data, dtype=np.float32)

    # Define budget and optimize_maxcover
    for i, blevel in enumerate(budget_levels):
        no_blevel = i + 1
        prefix = utils.get_iteration_prexix(no_blevel, len(budget_levels))

        budget = blevel * cost.size
        llogger.info("{} Optimizing with ".format(prefix) +
                     "budget level {}...".format(blevel))
        x = optimize_maxcover(cost, budget, sum_array, verbose=verbose,
                              logger=llogger)
        # Create a full (filled with 0s) raster template
        x_selection = np.full((profile["height"], profile["width"]), 0.0)
        # Place the values of result array (it's binary = {0, 1}) into template
        # elements that are False in the original mask
        x_selection[~mask] = x
        # Add the selected elements (planning units) into the selection
        # frequency matrix
        sel_freq += x_selection

        if save_intermediate:
            # Replace the real nodata values with a proper NoData value
            nodata_value = 255
            x_selection[mask] = nodata_value
            # Create a masked array
            output_x = ma.masked_values(x_selection, nodata_value)
            # Construct the output raster name
            btoken = str(blevel).replace('.', '_')
            # Construct a subdir name based on the basename
            output_bname = os.path.basename(output_rank_raster).split('.')[0]
            output_subir = os.path.join(os.path.dirname(output_rank_raster),
                                        output_bname)
            if not os.path.exists(output_subir):
                os.makedirs(output_subir)
            output_raster = os.path.join(output_subir,
                                         "budget_level_{}.tif".format(btoken))
            # Write out the data
            llogger.debug(" Writing intermediate output to " +
                          "{}".format(output_raster))
            profile.update(dtype=rasterio.uint8, compress=compress,
                           nodata=nodata_value)

            with rasterio.open(output_raster, 'w', **profile) as dst:
                dst.write_mask(mask)
                dst.write(output_x.astype(np.uint8), 1)

    opt_end = timer()
    opt_elapsed = round(opt_end - opt_start, 2)
    llogger.info(" [TIME] Optimization took {} sec".format(opt_elapsed))

    # 4. Rank values ---------------------------------------------------------

    post_start = timer()
    llogger.info(" [** POST-PROCESSING **]")

    llogger.info(" [1/3] Ranking values (this can take a while...)")
    # Use 0s from summation as a mask
    rank_array = ma.masked_values(sel_freq, 0.0)
    rank_array = rankdata(rank_array)

    # 5. Recale data into range [0, 1] ---------------------------------------

    llogger.info(" [2/3] Rescaling ranks")
    rank_array = spatutils.normalize(rank_array)
    # Force float32
    rank_array = rank_array.astype(np.float32)

    # 6. Prepare and write output --------------------------------------------

    llogger.info(" [3/3] Writing output to {}".format(output_rank_raster))
    # Replace the real nodata values with a proper NoData value
    nodata_value = -3.4e+38
    rank_array[mask] = nodata_value
    # Create a masked array
    rank_array = ma.masked_values(rank_array, nodata_value)
    profile.update(dtype=rasterio.float32, compress=compress,
                   nodata=nodata_value)

    with rasterio.open(output_rank_raster, 'w', **profile) as dst:
        dst.write_mask(mask)
        dst.write(rank_array.astype(np.float32), 1)

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
@click.option('-b', '--budget', required=True,
              help='Budget for the maximum coverage problem')
@click.argument('indir', nargs=1, type=click.Path(exists=True))
@click.argument('outfile', nargs=1)
def cli(indir, outfile, budget, extension, verbose):
    """ Command-line interface."""
    # List files in input directorys
    file_iterator = glob.iglob('{0}/**/*{1}'.format(indir, extension),
                               recursive=True)
    input_rasters = [item for item in file_iterator]
    if verbose:
        click.echo(click.style(" Found {} rasters".format(len(input_rasters)),
                               fg='green'))
    optimize_gurobi(input_rasters, outfile, budget, verbose=verbose)


if __name__ == '__main__':
    cli()

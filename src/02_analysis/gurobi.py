#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Functions and utilities for calculating the rarity-weighted richness score.

Module can be used alone or as part of Snakemake workflow.
"""

import click
import logging
import numpy as np
import numpy.ma as ma
import pdb
import rasterio

from gurobipy import *
from importlib.machinery import SourceFileLoader
from scipy.stats.mstats import rankdata
from timeit import default_timer as timer

utils = SourceFileLoader("lib.utils", "src/00_lib/utils.py").load_module()
spatutils = SourceFileLoader("lib.spatutils", "src/00_lib/spatutils.py").load_module()


def optimize_maxcover(cost, fraction, rij, normalize=False, verbose=False,
                      logger=None):
    """ Solve maximum coverage problem using Gurobi.

    :param cost: Numpy ndarray of planning unit costs. In this case the spatial
                 representation of the planning units is not provided.
    :param fraction: numeric indicating the fraction of the landscape that's
                     the target for the maximum coverage problem.
    :param rij: Numpy ndarray; Column sums of the matrix of representation
                levels of conservation features (rows) within planning units
                (columns).
    :param normalize: Boolean defining whehter cost and rij arrays are
                      normalized in range [0, 1] before optimization. This may
                      be necessary since the multiobjective optimization uses
                      blending approach.
    :param verbose: Boolean indicating how much information is printed out.
    :param logger: logger object to be used.
    :return
    """

    # cost array contains the actual cost values per pixel. Create another
    # array that
    assert rij.size == cost.size, 'Representation and cost array sizes must match'

    # Fractions of elements out of the total number of elements
    elem_fractions = np.full(rij.size, 1 / rij.size, dtype=np.float)
    # Calculate cumsum to define what is the closest fraction multiple to
    # desired fraction
    cs_elem_fractions = np.cumsum(elem_fractions)
    # Find the index of closest cumsum value to the requested fraction
    fraction_approx = utils.find_nearest(cs_elem_fractions, fraction)

    # NO COSTS USED
    if np.max(cost) == np.min(cost):
        try:
            logger.debug(" [DEBUG] Starting single-objective optimization")
            # Create a new model
            model = Model("singleobj")

            # Create variables
            for element in np.nditer(rij):
                model.addVar(vtype=GRB.BINARY)

            # Integrate new variables
            model.update()
            vars = model.getVars()

            # Set objective and constraint
            obj = LinExpr()
            expr = LinExpr()

            for i in range(rij.size):
                obj.addTerms(rij[i], vars[i])
                expr.addTerms(x[i], vars[i])

            model.setObjective(obj, sense=GRB.MAXIMIZE)
            model.addConstr(lhs=expr, sense=GRB.LESS_EQUAL,
                            rhs=fraction_approx)
            model.update()

            # The the log file directory from the logger, but log into a
            # separate file. NOTE: handler location is hard coded.
            log_dir = os.path.dirname(logger.handlers[1].baseFilename)
            log_file = os.path.join(log_dir, "gurobi.log")
            if os.path.exists(log_dir):
                logger.debug("See {} for Gurobi log".format(log_file))
                model.params.logToConsole = 0
                model.params.logFile = log_file
            else:
                logger.debug("Log dir {} not found".format(log_dir))
                model.params.logToConsole = 0
            # model.params.mipgap = 0.001
            # model.params.solutionLimit = 1
            # model.params.presolve = -1
            model.optimize()

        except GurobiError:
            raise
    # COSTS USED
    else:
        if normalize:
            logger.debug(" [DEBUG] Normalizing cost and rij arrays")
            rij = spatutils.normalize(rij)
            cost = spatutils.normalize(cost)
        try:
            logger.debug(" [DEBUG] Starting multi-objective optimization")
            # Create initial model
            model = Model('multiobj')

            # The the log file directory from the logger, but log into a
            # separate file. NOTE: handler location is hard coded.
            log_dir = os.path.dirname(logger.handlers[1].baseFilename)
            log_file = os.path.join(log_dir, "gurobi.log")
            if os.path.exists(log_dir):
                logger.debug("See {} for Gurobi log".format(log_file))
                model.params.logToConsole = 0
                model.params.logFile = log_file
            else:
                logger.debug("Log dir {} not found".format(log_dir))
                model.params.logToConsole = 0
                # model.params.mipgap = 0.001
                # model.params.solutionLimit = 1
                # model.params.presolve = -1

            # Initialize decision variables for ground set:
            for element in np.nditer(rij):
                model.addVar(vtype=GRB.BINARY)
            # Integrate new variables
            model.update()
            vars = model.getVars()

            # Constraint: limit total number of elements to be picked to be at
            # most fraction of total number of elements
            frac_constr = LinExpr()
            frac_constr.addTerms(elem_fractions.tolist(), vars)
            model.addConstr(lhs=frac_constr, sense=GRB.EQUAL,
                            rhs=fraction_approx)

            # Set global sense for ALL objectives
            model.ModelSense = GRB.MAXIMIZE

            # Limit how many solutions to collect
            model.setParam(GRB.Param.PoolSolutions, 100)

            # Set number of objectives
            model.NumObj = 2

            # Objective 1: maximize values
            model.setParam(GRB.Param.ObjNumber, 0)
            model.ObjNWeight = 1
            model.ObjNName = "Values"
            # model.ObjNRelTol = 0.01
            # model.ObjNAbsTol = 1.0
            model.setAttr(GRB.Attr.ObjN, vars, rij.tolist())

            # Objective 2: minimize values
            model.setParam(GRB.Param.ObjNumber, 1)
            # Minimizing an objective function is equivalent to maximizing the
            # negation of that function -> use ObjNWeight = -1.0
            model.ObjNWeight = -1.0
            model.ObjNName = "Costs"
            # model.ObjNRelTol = 0.01
            # model.ObjNAbsTol = 2.0
            model.setAttr(GRB.Attr.ObjN, vars, cost.tolist())

            model.optimize()

        except GurobiError:
            raise

    # Check the number of solutions
    n_solutions = model.SolCount
    if n_solutions > 1:
        logger.warning( " [WARN] {} solutions found, using the first (best) solution".format(n_solutions))

    # Construct a result array and return that
    res = np.asarray([var.x for var in model.getVars()], dtype=np.uint8)
    return res


def prioritize_gurobi(input_rasters, output_rank_raster, cost_raster=None,
                      step=0.05, save_intermediate=False, compress='DEFLATE',
                      ol_normalize=False, weights=None, verbose=False,
                      logger=None):
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

    It is possible to provide a list (vector) of weights for each features.
    These values are used as simple multipliers for each feature when summing
    the values over all features. If provided, the list must match the number
    of input rasters exactly.

    :param input_rasters: List of String paths of input rasters.
    :param output_raster: String path to the rank raster file to be created.
    :param cost_rasters: String path to a raster defining the cost surface. If None
                 (default), an equal cost is used for all cells.
    :param step: numeric value in (0, 1) defining the step length for budget
                 levels.
    :param save_intermediate: should intermediate optiomization results be
                              saved?
    :param compress: String compression level used for the output raster.
    :param ol_normalize: Boolean setting OL (Occurrence Level) normalization.
    :param weights: list of weights. Length must match the number of input
                    rasters.
    :param verbose Boolean indicating how much information is printed out.
    :param logger logger object to be used.
    """
    # 1. Setup  --------------------------------------------------------------

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
    assert len(input_rasters) > 0, "Input rasters list cannot be empty"
    assert len(output_rank_raster) != "", "Output raster path cannot be empty"
    try:
        step = float(step)
    except ValueError:
        llogger.error("'step' must be coercible to float")
        sys.exit(1)
    assert step > 0.0 and step < 1.0, "Step argument must be in range (0, 1)"
    # Construct fraction levels based on the step provided. 0.0 (nothing) and
    # 1.0 (everything) are not needed.
    fraction_levels = np.linspace(0.0+step, 1.0, 1/step)

    # 2. Pre-processing  -----------------------------------------------------

    llogger.info(" [** PRE-PROCESSING **]")

    # Create a sum array, i.e. sum all (occurrence level normalized) raster
    # values in input_rasters together.  NOTE: sum_raster() returns a masked
    # array.
    sum_array_masked = spatutils.sum_raster(input_rasters, olnormalize=True,
                                            weights=weights, logger=llogger)
    # To speed up things, do 2 things: 1) save mask (NoData) and get rid of
    # NoData cells for now, 2) flatten the array.
    (height, width) = sum_array_masked.shape
    mask = ma.getmask(sum_array_masked)

    # If no cost feature is provided, use equal cost for everything
    if cost_raster is not None:
        cost_array_masked = rasterio.open(cost_raster, 'r').read(1, masked=True)
        # Get the mask to see if it's different to sum_array
        cost_mask = ma.getmask(cost_array_masked)
        # Check that cost value is available for all sum_array cells
        if not np.all(cost_mask == mask):
            mask_diff = np.sum(cost_mask) - np.sum(mask)
            llogger.warning(" [WARNING] Different masks for cost and feature sum rasters, diff: {}".format(mask_diff))
            llogger.warning(" [WARNING] Using data only from cells that have both cost and occurrence data")
            # Get an intersection of the cost and sum array masks
            intersect_mask = np.logical_and(~cost_mask, ~mask)
            # Set the intersection mask to both cost and sum array data. NOTE: no
            # need to set masked values since using intersection.
            sum_array_masked.mask = ~intersect_mask
            cost_array_masked.mask = ~intersect_mask
            mask = ~intersect_mask

        sum_array = ma.compressed(sum_array_masked)
        cost_array = ma.compressed(cost_array_masked)
    else:
        # Get all the non-masked data as a 1-D array
        sum_array = ma.compressed(sum_array_masked)
        cost_array = np.ones(sum_array.size)

    load_end = timer()
    load_elapsed = round(load_end - load_start, 2)
    llogger.info(" [TIME] Pre-processing took {} sec".format(load_elapsed))

    # 3. Optimize  -----------------------------------------------------------

    opt_start = timer()
    llogger.info(" [** OPTIMIZING **]")
    flevels_str = ", ".join([str(level) for level in fraction_levels])
    llogger.info(" [NOTE] Target budget levels: {}".format(flevels_str))

    # Construct a ndarray (matrix) that will hold the selection frequency.
    # Populate it with 0s
    sel_freq = np.full((height, width), 0, dtype=rasterio.float64)

    # Define budget and optimize_maxcover
    for i, flevel in enumerate(fraction_levels):
        no_flevel = i + 1
        prefix = utils.get_iteration_prefix(no_flevel, len(fraction_levels))

        llogger.info("{} Optimizing with ".format(prefix) +
                     "fraction level {}...".format(np.round(flevel, 2)))
        x = optimize_maxcover(cost=cost_array, fraction=flevel, rij=sum_array,
                              normalize=True, verbose=verbose, logger=llogger)
        # Create a full (filled with 0s) raster template
        x_selection = np.full((height, width), 0.0, dtype=rasterio.float64)
        # Place the values of result array (it's binary = {0, 1}) into template
        # elements that are False in the original mask
        x_selection[~mask] = x
        # Add the selected elements (planning units) into the selection
        # frequency matrix
        sel_freq += x_selection
        # Get the raster profile from the input raster files
        profile = spatutils.get_profile(input_rasters, logger=llogger)

        if save_intermediate:
            # Replace the real nodata values with a proper NoData value
            nodata_value = 255
            x_selection[mask] = nodata_value
            # Create a masked array
            output_x = ma.masked_values(x_selection, nodata_value)
            # Construct the output raster name
            ftoken = str(flevel).replace('.', '_')
            # Construct a subdir name based on the basename
            output_bname = os.path.basename(output_rank_raster).split('.')[0]
            output_subir = os.path.join(os.path.dirname(output_rank_raster),
                                        output_bname)
            if not os.path.exists(output_subir):
                os.makedirs(output_subir)
            output_raster = os.path.join(output_subir,
                                         "budget_level_{}.tif".format(ftoken))
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

    llogger.info(" [1/3] Ranking selection frequencies")
    # Use 0s from summation as a mask
    rank_array = ma.masked_values(sel_freq, 0.0)
    rank_array = rankdata(rank_array)

    # 5. Recale data into range [0, 1] ---------------------------------------

    llogger.info(" [2/3] Rescaling ranks")
    rank_array = spatutils.normalize(rank_array)

    # 6. Prepare and write output --------------------------------------------

    llogger.info(" [3/3] Writing output to {}".format(output_rank_raster))
    # Replace the real nodata values with a proper NoData value
    nodata_value = -3.4e+38
    rank_array[mask] = nodata_value
    # Create a masked array
    rank_array = ma.masked_values(rank_array, nodata_value)
    profile.update(dtype=rasterio.float64, compress=compress,
                   nodata=nodata_value)

    with rasterio.open(output_rank_raster, 'w', **profile) as dst:
        dst.write_mask(mask)
        dst.write(rank_array.astype(np.float64), 1)

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

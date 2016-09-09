#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Functions and utilities comparing raster similarities.

Module can be used alone or as part of Snakemake workflow.
"""
import rasterio
import pandas as pd
import numpy as np

from scipy.spatial.distance import jaccard


def jaccard(x, y, x_min=0.0, x_max=1.0, y_min=0.0, y_max=1.0,
            warn_uneven=False, limit_tolerance=4, disable_checks=FALSE):
    """Calculate the Jaccard coefficient.

    The Jaccard coefficient measures similarity between sample sets, and is
    defined as the size of the intersection divided by the size of the union of
    the sample sets. The Jaccard coefficient can be calculated for a subset of
    rasters provided by using the threshold argument.

    Min and max values must be provided for both RasterLayer objects x
    and y. Method can be used with RasterLayers of any value range, but
    the defaults [0.0, 1.0] are geared towards comparing Zonation rank priority
    rasters. Limits provided are inclusive.

    :param x ndarray object.
    :param y ndarray object.
    :param x_min Numeric minimum threshold value for x to be used
                 (default 0.0).
    :param x_max Numeric maximum threshold value for x to be used
                 (default 1.0).
    :param y_min Numeric minimum threshold value for y to be used
                 (default 0.0).
    :param y_max Numeric maximum threshold value for y to be used
                 (default 1.0).
    :param warn_uneven Boolean indicating whether a warning is raised if the
                       compared raster coverages are very (>20x) uneven.
    :param limit_tolerance integer values that defines to which precision x and
                           y limits are rounded to. This helps e.g. with values
                           that close to 0 but not quite 0 (default: 4, i.e.
                           round(x, 4)).
    :param disable_checks boolean indicating if the input limit values are
                          checked against the actual raster values in x and y.

    :return numeric value in [0, 1].
    """

    if not disable_checks:
        assert x_min < np.round(np.min(x), limit_tolerance), "Min threshold smaller than computed min of x"
        assert x_max < np.round(np.max(x), limit_tolerance), "Max threshold smaller than computed max of x"
        assert x_min >= x_max, "Min threshold for x larger or equal to max threshold"
        assert y_min < np.round(np.min(y), limit_tolerance), "Min threshold smaller than computed min of y"
        assert y_max < np.round(np.max(y), limit_tolerance), "Max threshold smaller than computed max of y"
        assert y_min >= y_max, "Min threshold for y larger or equal to max threshold"

  # [fixme] - using cellStats(X, "sum") should be safe as we're dealing with
  # binary 0/1 rasters. count() would be preferable, but apparently raster
  # (>= 2.2 at least) doesn't support it anymore.

  # Get the values according to the limits provided
  x_bin <- (x >= x.min & x <=x.max)
  y_bin <- (y >= y.min & y <=y.max)

  if (warn.uneven) {
    x.size <- cellStats(x.bin, "sum")
    y.size <- cellStats(y.bin, "sum")
    # Sort from smaller to larger
    sizes <- sort(c(x.size, y.size))
    if (sizes[2] / sizes[1] > 20) {
      warning("The extents of raster values above the threshhold differ more",
              "than 20-fold: Jaccard coefficient may not be informative.")
    }
  }

  # Calculate the intersection of the two rasters, this is given by adding
  # the binary rasters together -> 2 indicates intersection
  combination <- x.bin + y.bin
  intersection <- combination == 2

  # Union is all the area covered by the both rasters
  union <- combination >= 1

  return(cellStats(intersection, "sum") / cellStats(union, "sum"))
}


def cross_jaccard(input_rasters, thresholds, verbose=False, logger=None,
                  *args, **kwargs):
    """ Calculate Jaccard coefficients bewteen all the inpur rasters.

    This is a utility function that is intented to be used to compare
    top-fractions of the landscape. Thus, x_max and y_max for
    jaccard are fixed to 1.0.

    :param input_rasters list of input raster paths.
    :param thresholds Numeric vector values of thresholds.
    :param verbose: Boolean indicating how much information is printed out.
    :param logger: logger object to be used.
    :param ... additional arguments passed on to jaccard().

    :return dict of Pandas Dataframes with Jaccard coefficients between all
            rasters.
    """
    # 1. Setup  --------------------------------------------------------------

    all_start = timer()
    load_start = timer()

    if not logger:
        logging.basicConfig()
        llogger = logging.getLogger('cross_jaccard')
        llogger.setLevel(logging.DEBUG if verbose else logging.INFO)
    else:
        llogger = logger

    # Check the inputs
    assert len(input_rasters) > 1, "More than one input rasters are needed"

    # 2. Calculations --------------------------------------------------------

    all_jaccards = {}
    n_rasters = len(input_rasters)

    for threshold in thresholds:
        # Initialize a matrix to hold the jaccard coefficients and populate it
        # with -1s.
        jaccards = np.empty([n_rasters, n_rasters])
        jaccards[:] = -1.0

        for i in range(0, n_rasters):
            raster1 = input_rasters[i]
            raster1_src = rasterio.open(raster1).read(1)
            for j in range(0, n_rasters):
                raster2 = input_rasters[j]
                raster2_src = rasterio.open(raster2).read(1)
                # Jaccard coefficient is always 1.0 on the diagonal
                if i == j:
                    jaccards[i, j] = 1.0
                else:
                    # See the complement, if it's not NA then the pair has
                    # already been compared
                    if jaccards[j, i] != -1.0:
                        logger.info(("Calculating Jaccard index for [{}".format(threshold) +
                                 ", 1.0, ] between {0} and {1}".format(raster1,
                                                                       raster2)))

                        jaccards[i, j] = jaccard(raster1_src, raster2_src,
                                                 x_min=threshold, x_max=1.0,
                                                 y_min=threshold, y_max=1.0,
                                                 args, kwargs)
                    else:
                        jaccards[i, j] = jaccards[j, i]

        jaccards = pd.DataFrame(jaccards)
        jaccards.columns = input_rasters
        jaccards.index = input_rasters
        all_jaccards[threshold] = jaccards

    return all_jaccards

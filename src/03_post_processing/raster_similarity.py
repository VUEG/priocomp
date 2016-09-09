#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Functions and utilities comparing raster similarities.

Module can be used alone or as part of Snakemake workflow.
"""
import rasterio
import pandas as pd
import numpy as np

from scipy.spatial.distance import jaccard


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

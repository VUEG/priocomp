#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import rasterio


def normalize(x):
    """ Rescale all numeric values in range [0, 1].

        Input must be a numpy ndarray, no coercion is tried.

        :param x: numpy ndarray to be rescaled.
        :return: numpy ndarray with rescaled values.
        """
    if type(x) is not np.ndarray:
        raise TypeError("x must be a numpy.ndarray")

    # true_divide() is needed to handle both floats and integers
    np.true_divide(x, np.max(np.abs(x), axis=0), out=x, casting='unsafe')
    return x


def standardize(x):
    pass


def rescale_raster(input_raster, output_raster, method, compression='DEFLATE'):
    """ Rescale all numeric values of a raster according ot a given method.

    Currently two methods are implemented:
        1. 'normalize'
        2. 'standardize'

    :param input_raster: String path to raster file to be normalized.
    :param output_raster: String path to raster file to be created.
    :param compression: String compression level used for the output raster
    :return: Boolean indicating whether the operation was successful
    """
    return None


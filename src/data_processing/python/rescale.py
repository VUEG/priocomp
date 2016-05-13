#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import os
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
    np.true_divide(x, np.max(np.abs(x)), out=x, casting='unsafe')
    return x


def standardize(x):
    pass


def rescale_raster(input_raster, output_raster, method, compress='DEFLATE'):
    """ Rescale all numeric values of a raster according ot a given method.

    Currently two methods are implemented:
        1. 'normalize'
        2. 'standardize'

    :param input_raster: String path to raster file to be normalized.
    :param output_raster: String path to raster file to be created.
    :param compress: String compression level used for the output raster
    """
    if not os.path.exists(input_raster):
        raise OSError("Input raster {} not found".format(input_raster))

    with rasterio.open(input_raster) as in_src:
        # Read the first band
        src_data = in_src.read(1)
        if method == 'normalize':
            rescaled_data = normalize(src_data)
        elif method == 'standardize':
            rescaled_data = standardize(src_data)
        else:
            raise TypeError("Method {} not implemented".format(method))

        # Write the product.
        profile = rasterio.default_gtiff_profile()
        # Rescaled data is always float32, and we have only 1 band
        profile.update(dtype=rasterio.float32, count=1, compress=compress, height=rescaled_data.shape[0],
                       width=rescaled_data.shape[1])

        with rasterio.open(output_raster, 'w', **profile) as dst:
            dst.write(rescaled_data.astype(rasterio.float32), 1)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Utility functions for spatial processing."""
import click
import logging
import numpy as np
import numpy.ma as ma
import os
import rasterio
import scipy.stats

from importlib.machinery import SourceFileLoader

utils = SourceFileLoader("src.utils", "src/utils.py").load_module()


def get_profile(input_rasters, warn_inconistent=True, pick="first",
                logger=None):
    """ Construct a raster profile base on a list of input raster files.

    Function will loop over a list of input raster files, extract raster
    profile (metadata), possibly check for inconsistencies and return
    a profile. The profile returned is always that of the raster last in the
    list.

    :param input_rasters: String list of raster paths.
    :param warn_inconistent: Boolean should inconsistent profiles cause a
                             warning?
    :param pick: String of either value "first" or "last". Former picks the
                 profile from the first raster in the list, "last" the last.
    :param logger: Logger object.
    :return: Dict profile (raster metadata).
    """
    # Set up logging
    if not logger:
        logging.basicConfig()
        llogger = logging.getLogger('get_profile')
        llogger.setLevel(logging.INFO)
    else:
        llogger = logger

    if pick not in ["first", "last"]:
        raise ValueError("Value for 'pick' must be either 'first' or 'last'")

    current_profile = None
    for i, raster in enumerate(input_rasters):
        with rasterio.open(raster) as src:
            profile = src.profile
            if i == 0:
                current_profile = profile
            else:
                if warn_inconistent and profile != current_profile:
                    llogger.warning(" [WARN] Inconsistent raster profiles " +
                                    "detected")
                    warn_inconistent = False
                if pick == 'last':
                    current_profile = profile
    return current_profile


def normalize(x):
    """ Rescale all numeric values in range [0, 1].

    Input must be a numpy ndarray, no coercion is tried.

    :param x: numpy ndarray to be rescaled.
    :return: numpy ndarray with rescaled values.
    """
    if type(x) is not np.ndarray and type(x) is not ma.core.MaskedArray:
        raise TypeError("x must be a numpy.ndarray or numpy.ma.MaskedArray")

    # NOTE: the approach commented out below would be more memory efficient,
    # but doesn't work as such with masked arrays
    # np.true_divide(x, np.max(np.abs(x)), out=x, casting='unsafe')
    # Data may have negative values, thus first add the abs(min(x)) to
    # everything.
    x_min = ma.min(x)
    x_max = ma.max(x)
    return (x - x_min) / (x_max - x_min)


def ol_normalize(x):
    """ Normalize layer based on occurrence levels in the array.

    The value of each element is divided by the sum off all elements. Input
    must be a numpy ndarray, no coercion is tried.

    :param x: numpy ndarray to be rescaled.
    :return: numpy ndarray with transformed values.
    """
    if type(x) is not np.ndarray and type(x) is not ma.core.MaskedArray:
        raise TypeError("x must be a numpy.ndarray or numpy.ma.MaskedArray")

    min_val = ma.min(x)

    return (x - min_val) / (ma.sum(x - min_val))


def rescale_raster(input_raster, output_raster, method, fill_w_zeros=False,
                   only_positive=False, compress='DEFLATE', verbose=False,
                   logger=None):
    """ Rescale all numeric values of a raster according ot a given method.

    Currently two methods are implemented:
        1. 'normalize'
        2. 'normalize_ol'
        3. 'standardize'

    :param input_raster: String path to raster file to be normalized.
    :param output_raster: String path to raster file to be created.
    :param method: String method to use.
    :param fill_w_zeros:  Boolean indicating whether NoData is encoded as real
                          NoData or filled with zeros.
    :param only_positive: Boolean indicating if negative values (in the
                          original raster) should be ignored.
    :param compress: String compression level used for the output raster.
    :param verbose: Boolean indicating how much information is printed out.
    :param logger: Logger object.
    :return Boolean True if success, False otherwise
    """
    # Set up logging
    if not logger:
        logging.basicConfig()
        llogger = logging.getLogger('rescale_raster')
        llogger.setLevel(logging.DEBUG if verbose else logging.INFO)
    else:
        llogger = logger
    NODATA_VALUE = -3.4e+38
    llogger.debug("Using internal NoData value {}".format(NODATA_VALUE))

    if not os.path.exists(input_raster):
        raise OSError("Input raster {} not found".format(input_raster))

    with rasterio.open(input_raster) as in_src:
        # Read the first band
        src_data = in_src.read(1, masked=True)
        # Cast data to float32 as that's what's used anyways
        src_data = src_data.astype(np.float32)
        # Check min and max
        decimals = 12
        src_min = np.round(src_data.min(), decimals)
        src_max = np.round(src_data.max(), decimals)
        if src_min == 0.0 and src_max == 0.0:
            llogger.warning('{} has all zero values, '.format(input_raster +
                            'skipping'))
            return False

        # Get information on the distribution of values before rescaling
        src_q75, src_q25 = np.nanpercentile(ma.filled(src_data, np.nan),
                                            (75, 25))
        src_q25 = np.round(src_q25, decimals)
        src_q75 = np.round(src_q75, decimals)
        src_mean = np.round(src_data.mean(), decimals)
        llogger.debug(" Before rescaling (min, q25, mean, q75, max): " +
                      "{0}, {1}, {2}, {3}, {4}".format(src_min, src_q25,
                                                       src_mean, src_q75,
                                                       src_max))
        # Only positive values?
        if only_positive:
            if src_min < 0:
                llogger.info("Ignoring negative values")
                src_data = src_data.clip(min=0)
            else:
                llogger.warning("Ignore negative values requested but " +
                                " there are none")

        # Do the actual rescaling
        if method == 'normalize':
            rescaled_data = normalize(src_data)
        elif method == 'ol_normalize':
            rescaled_data = ol_normalize(src_data)
        elif method == 'standardize':
            pass
        else:
            raise TypeError("Method {} not implemented".format(method))

        dst_q75, dst_q25 = np.nanpercentile(ma.filled(rescaled_data, np.nan),
                                            (75, 25))
        dst_q25 = np.round(dst_q25, decimals)
        dst_q75 = np.round(dst_q75, decimals)
        dst_min = np.round(rescaled_data.min(), decimals)
        dst_mean = np.round(rescaled_data.mean(), decimals)
        dst_max = np.round(rescaled_data.max(), decimals)
        llogger.debug("  After rescaling (min, q25, mean, q75, max): " +
                      "{0}, {1}, {2}, {3}, {4}".format(dst_min, dst_q25,
                                                       dst_mean, dst_q75,
                                                       dst_max))

        # Write the product.
        profile = in_src.profile
        # Rescaled data is always float32, and we have only 1 band.
        profile.update(dtype=rasterio.float32, compress=compress,
                       nodata=NODATA_VALUE)
        # Replace the NoData values and set the fill value for NoData, unless
        # the result is to be filled with zeros
        if fill_w_zeros:
            rescaled_data = ma.filled(rescaled_data, 0)
            llogger.debug("Filling NoData with zeros")
        else:
            rescaled_data = ma.filled(rescaled_data, NODATA_VALUE)

        with rasterio.open(output_raster, 'w', **profile) as dst:
            if not fill_w_zeros:
                # Set the NoData mask
                dst.write_mask(in_src.read_masks(1))
            dst.write(rescaled_data, 1)
            llogger.debug("Wrote raster {}".format(output_raster))
        return True


def sum_raster(input_rasters, olnormalize=False, verbose=False, logger=None):
    """ Sum a group of input rasters.

    Read in a group of input rasters and sum the values in each cell. All
    rasters must have the same dimensions, no broadcasting allowed.

    Optionally, values in all rasters can occurrence level (OL) normalized
    before summation.

    :param input_rasters: String list of input raster paths.
    :param ol_normalize: Boolean indicating wether OL normalization is done.
    :param verbose: Boolean indicating how much information is printed out.
    :param logger: Logger object.
    :return: numpy ndarray of summed values.
    """
    # Set up logging
    if not logger:
        logging.basicConfig()
        llogger = logging.getLogger('sum_raster')
        llogger.setLevel(logging.DEBUG if verbose else logging.INFO)
    else:
        llogger = logger

    if ol_normalize:
        llogger.info(" [NOTE] Using occurrence level normalization")

    # Initiate the data structure for holding the summed values.
    sum_array = None
    n_rasters = len(input_rasters)

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

            if olnormalize:
                # 1. Occurrence level normalize data --------------------------
                llogger.debug("{0} OL normalizing data".format(prefix))
                src_data = ol_normalize(src_data)

            # Sum OL normalized data ---------------------------------------
            llogger.debug("{0} Summing values".format(prefix))
            # If this is the first raster, use its dimensions to build an array
            # that holds the summed values.
            if i == 0:
                # Set up an array of zeros that has the correct dimensions
                sum_array = np.zeros_like(src_data, dtype=np.float32)
            sum_array += src_data

    return sum_array


def winsorize(x, limits=0.05):
    """ Winsorize values in an array.

    The distribution of many statistics can be heavily influenced by outliers.
    A typical strategy is to set all outliers to a specified percentile of the
    data; for example, a 90% winsorization would see all data below the 5th
    percentile set to the 5th percentile, and data above the 95th percentile
    set to the 95th percentile. Winsorized estimators are usually more robust
    to outliers than their more standard forms.

    See https://en.wikipedia.org/wiki/Winsorizing

    :param x: numpy ndarray to be rescaled.
    :param limits: Tuple of the percentages to cut on each side of the array,
                   with respect to the number of unmasked data, as floats
                   between 0. and 1. Noting n the number of unmasked data
                   before trimming, the (n*limits[0])th smallest data and the
                   (n*limits[1])th largest data are masked, and the total
                   number of unmasked data after trimming is n*(1.-sum(limits))
                   The value of one limit can be set to None to indicate an
                   open interval.
    :return: numpy ndarray with transformed values.
    """
    if type(x) is not np.ndarray and type(x) is not ma.core.MaskedArray:
        raise TypeError("x must be a numpy.ndarray or numpy.ma.MaskedArray")

    return scipy.stats.mstats.winsorize(x, limits=limits)


@click.command()
@click.option('-m', '--method', default='normalize',
              help='Rescaling method used.')
@click.option('-v', '--verbose', is_flag=True)
@click.argument('infile', nargs=1, type=click.Path(exists=True))
@click.argument('outfile', nargs=1)
def cli(infile, outfile, method, verbose):
    """ Command-line interface."""
    click.echo(click.style('Rescaling (method: {0}) file {1}'.format(method,
                                                                     infile),
                           fg='green'))
    success = rescale_raster(infile, outfile, method=method, verbose=verbose)
    if success:
        click.echo(click.style('Done!', fg='green'))
    else:
        click.echo(click.style('Normalization failed', fg='red'))


if __name__ == '__main__':
    cli()

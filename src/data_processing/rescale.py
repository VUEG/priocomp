#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Rescaling and normalization functions."""
import click
import numpy as np
import numpy.ma as ma
import os
import rasterio

from importlib.machinery import SourceFileLoader

utils = SourceFileLoader("src.utils", "src/utils.py").load_module()


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
    x = x + ma.abs(ma.min(x))
    x = x / (ma.max(ma.abs(x)))
    return x


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


def rescale_raster(input_raster, output_raster, method, compress='DEFLATE',
                   verbose=False, log_file=None):
    """ Rescale all numeric values of a raster according ot a given method.

    Currently two methods are implemented:
        1. 'normalize'
        2. 'normalize_ol'
        3. 'standardize'

    :param input_raster: String path to raster file to be normalized.
    :param output_raster: String path to raster file to be created.
    :param compress: String compression level used for the output raster.
    :param verbose Boolean indicating how much information is printed out.
    :param log_file String path to log file used.
    :return Boolean True if success, False otherwise
    """
    # Set up logging
    llogger = utils.get_local_logger(rescale_raster.__name__, log_file,
                                     debug=verbose)
    NODATA_VALUE = -3.4e+38
    llogger.debug("Using internal NoData value {}".NODATA_VALUE)

    if not os.path.exists(input_raster):
        raise OSError("Input raster {} not found".format(input_raster))

    with rasterio.open(input_raster) as in_src:
        # Read the first band
        src_data = in_src.read(1, masked=True)
        # Cast data to float32 as that's what's used anyways
        src_data = src_data.astype(np.float32)
        # Check min and max
        src_min = src_data.min()
        src_max = src_data.max()
        if src_min == 0.0 and src_max == 0.0:
            llogger.warning('{} has all zero values, '.format(input_raster +
                            'skipping'))
            return False
        if verbose:
            q75, q25 = np.nanpercentile(ma.filled(src_data, np.nan), (75, 25))
            mean = src_data.mean()
            llogger.debug(' min before rescaling:     {}'.format(src_min))
            llogger.debug(' lower Q before rescaling: {}'.format(q25))
            llogger.debug(' mean before rescaling:    {}'.format(mean))
            llogger.debug(' upper Q before rescaling: {}'.format(q75))
            llogger.debug(' max before rescaling:     {}'.format(src_max))
        if method == 'normalize':
            rescaled_data = normalize(src_data)
        elif method == 'ol_normalize':
            rescaled_data = ol_normalize(src_data)
        elif method == 'standardize':
            pass
        else:
            raise TypeError("Method {} not implemented".format(method))

        if verbose:
            q75, q25 = np.nanpercentile(ma.filled(rescaled_data, np.nan),
                                        (75, 25))
            min = rescaled_data.min()
            mean = rescaled_data.mean()
            max = rescaled_data.max()
            llogger.debug(' min after rescaling:     {}'.format(min))
            llogger.debug(' lower Q after rescaling: {}'.format(q25))
            llogger.debug(' mean after rescaling:    {}'.format(mean))
            llogger.debug(' upper Q after rescaling: {}'.format(q75))
            llogger.debug(' max after rescaling:     {}'.format(max))
            llogger.debug(' saving data to:          {}'.format(output_raster))

        # Write the product.
        profile = in_src.profile
        # Rescaled data is always float32, and we have only 1 band.
        profile.update(dtype=rasterio.float32, compress=compress,
                       nodata=NODATA_VALUE)
        # Also remember to replace the NoData values and set the fill value for
        # NoData
        rescaled_data = ma.filled(rescaled_data, NODATA_VALUE)
        with rasterio.open(output_raster, 'w', **profile) as dst:
            # Set the NoData mask
            dst.write_mask(in_src.read_masks(1))
            dst.write(rescaled_data, 1)
            llogger.debug("Wrote raster {}".format(output_raster))
        return True


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

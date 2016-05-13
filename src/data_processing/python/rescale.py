#!/usr/bin/python3
# -*- coding: utf-8 -*-
import click
import numpy as np
import numpy.ma as ma
import os
import rasterio


def normalize(x):
    """ Rescale all numeric values in range [0, 1].

        Input must be a numpy ndarray, no coercion is tried.

        :param x: numpy ndarray to be rescaled.
        :return: numpy ndarray with rescaled values.
        """
    if type(x) is not np.ndarray and type(x) is not ma.core.MaskedArray:
        raise TypeError("x must be a numpy.ndarray or numpy.ma.MaskedArray")

    # NOTE: the approach commented out below would be more memory efficient, but
    # doesn't work as such with masked arrays
    # np.true_divide(x, np.max(np.abs(x)), out=x, casting='unsafe')
    # Data may have negative values, thus first add the abs(min(x)) to everything.
    x = x + ma.abs(ma.min(x))
    x = x / (ma.max(ma.abs(x)))
    return x


def standardize(x):
    pass


def rescale_raster(input_raster, output_raster, method, compress='DEFLATE', verbose=False):
    """ Rescale all numeric values of a raster according ot a given method.

    Currently two methods are implemented:
        1. 'normalize'
        2. 'standardize'

    :param input_raster: String path to raster file to be normalized.
    :param output_raster: String path to raster file to be created.
    :param compress: String compression level used for the output raster.
    :param verbose Boolean indicating how much information is printed out.
    """
    if not os.path.exists(input_raster):
        raise OSError("Input raster {} not found".format(input_raster))

    with rasterio.open(input_raster) as in_src:
        # Read the first band
        src_data = in_src.read(1, masked=True)
        # Cast data to float32 as that's what's used anyways
        src_data = src_data.astype(np.float32)
        if verbose:
            q75, q25 = np.nanpercentile(ma.filled(src_data, np.nan), (75, 25))
            click.echo(click.style(' min before rescaling:     {}'.format(src_data.min()), fg='green'))
            click.echo(click.style(' lower Q before rescaling: {}'.format(q25), fg='green'))
            click.echo(click.style(' mean before rescaling:    {}'.format(src_data.mean()), fg='green'))
            click.echo(click.style(' upper Q before rescaling: {}'.format(q75), fg='green'))
            click.echo(click.style(' max before rescaling:     {}'.format(src_data.max()), fg='green'))
        if method == 'normalize':
            rescaled_data = normalize(src_data)
        elif method == 'standardize':
            rescaled_data = standardize(src_data)
        else:
            raise TypeError("Method {} not implemented".format(method))

        if verbose:
            q75, q25 = np.nanpercentile(ma.filled(rescaled_data, np.nan), (75, 25))
            click.echo(click.style('\n min after rescaling:     {}'.format(rescaled_data.min()), fg='green'))
            click.echo(click.style(' lower Q after rescaling: {}'.format(q25), fg='green'))
            click.echo(click.style(' mean after rescaling:    {}'.format(rescaled_data.mean()), fg='green'))
            click.echo(click.style(' upper Q after rescaling: {}'.format(q75), fg='green'))
            click.echo(click.style(' max after rescaling:     {}'.format(rescaled_data.max()), fg='green'))

        # Write the product.
        profile = in_src.profile
        # Rescaled data is always float32, and we have only 1 band
        profile.update(dtype=rasterio.float32, compress=compress)

        with rasterio.open(output_raster, 'w', **profile) as dst:
            dst.write(rescaled_data, 1)

@click.command()
@click.option('-m', '--method', default='normalize', help='Rescaling method used.')
@click.option('-v', '--verbose', is_flag=True)
@click.argument('infile', nargs=1, type=click.Path(exists=True))
@click.argument('outfile', nargs=1)
def cli(infile, outfile, method, verbose):
    click.echo(click.style('Rescaling file {}'.format(infile), fg='green'))
    rescale_raster(infile, outfile, method=method, verbose=verbose)
    click.echo(click.style('Done!', fg='green'))


if __name__ == '__main__':
    cli()
#!/usr/bin/python3
# -*- coding: utf-8 -*-
import click
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
        profile = in_src.profile
        # Rescaled data is always float32, and we have only 1 band
        profile.update(dtype=rasterio.float32, compress=compress)

        with rasterio.open(output_raster, 'w', **profile) as dst:
            dst.write(rescaled_data.astype(rasterio.float32), 1)

@click.command()
@click.option('-m', '--method', default='normalize', help='Rescaling method used.')
@click.argument('infile', nargs=1, type=click.Path(exists=True))
@click.argument('outfile', nargs=1)
def cli(infile, outfile, method):
    click.echo(click.style('Rescaling file {}'.format(infile), fg='green'))
    rescale_raster(infile, outfile, method)
    click.echo(click.style('Done!', fg='green'))


if __name__ == '__main__':
    cli()
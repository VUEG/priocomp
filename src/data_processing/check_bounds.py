#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import click
import fnmatch
import logging
import os
import rasterio

# Set up logging
log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.DEBUG, format=log_fmt)
logger = logging.getLogger(__name__)

def find_files(path, extension):
    """ Simple utility script to find files with given file extension in path.

    """
    matches = []
    for root, dirnames, filenames in os.walk(path):
        for filename in fnmatch.filter(filenames, extension):
            matches.append(os.path.join(root, filename))
    return matches

def process_rasters(input_rasters):

    for raster in input_rasters:
        # Register GDAL format drivers and configuration options with a
        # context manager.
        with rasterio.drivers():
            # Read raster bands directly to Numpy arrays.
            with rasterio.open(raster) as src:
                logger.info("Bounds for {0}:\n   {1}".format(raster, src.bounds))
                print("\n")

@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.option('-f', '--filetype', default="*.tif", help='Raster file type')
def main(input_filepath, filetype="*.tif"):

    logger.info('Harmonising rasters in {}'.format(input_filepath))
    # Get a list of raster files
    raster_files = find_files(input_filepath, filetype)
    process_rasters(raster_files)

if __name__ == '__main__':

    # not used in this stub but often useful for finding various files
    project_dir = os.path.join(os.path.dirname(__file__), os.pardir, os.pardir)

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    #load_dotenv(find_dotenv())

    main()

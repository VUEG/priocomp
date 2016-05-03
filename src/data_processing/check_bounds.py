#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import click
import fnmatch
import logging
import os
import rasterio

# Set up logging
log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(format=log_fmt)
logger = logging.getLogger(__name__)

def find_files(path, extension):
    """ Simple utility script to find files with given file extension in path.

    """
    matches = []
    for root, dirnames, filenames in os.walk(path):
        for filename in fnmatch.filter(filenames, extension):
            matches.append(os.path.join(root, filename))
    return matches

def process_rasters(input_rasters, verbose=False):

    # Store min/max extents. NOTE: this makes sense only if the rasters are in
    # the same CRS or have the same resolution and map unit.

    # Store CRSs
    crss = {}
    dimensions = {}

    for raster in input_rasters:
        # Register GDAL format drivers and configuration options with a
        # context manager.
        with rasterio.drivers():
            # Read raster bands directly to Numpy arrays.
            with rasterio.open(raster) as src:
                if verbose:
                    click.echo(click.style("Bounds for {0}:".format(raster), fg="green"))
                    click.echo(click.style(" <left> = {0}".format(src.bounds[0]), fg="green"))
                    click.echo(click.style(" <bottom> = {0}".format(src.bounds[1]), fg="green"))
                    click.echo(click.style(" <right> = {0}".format(src.bounds[2]), fg="green"))
                    click.echo(click.style(" <top> = {0}\n".format(src.bounds[3]), fg="green"))

                # Get the CRS; for items > first test if the CRS is the same
                crs_string = str(src.crs)
                if crs_string in crss.keys():
                    crss[crs_string].append(raster)

                    if src.bounds[0] < dimensions[crs_string]['min_left']['value']:
                        dimensions[crs_string]['min_left'] = {'raster': raster, 'value': src.bounds[0]}
                    if src.bounds[1] < dimensions[crs_string]['min_bottom']['value']:
                        dimensions[crs_string]['min_bottom'] = {'raster': raster, 'value': src.bounds[1]}
                    if src.bounds[2] > dimensions[crs_string]['max_right']['value']:
                        dimensions[crs_string]['max_right'] = {'raster': raster, 'value': src.bounds[2]}
                    if src.bounds[3] > dimensions[crs_string]['max_top']['value']:
                        dimensions[crs_string]['max_top'] = {'raster': raster, 'value': src.bounds[3]}
                else:
                    crss[crs_string] = [raster]
                    # Check max dimensions
                    dimensions[crs_string] = {'min_left': {'raster': raster, 'value': src.bounds[0]},
                                              'min_bottom': {'raster': raster, 'value': src.bounds[1]},
                                              'max_right': {'raster': raster, 'value': src.bounds[2]},
                                              'max_top': {'raster': raster, 'value': src.bounds[3]}}


    if len(crss.keys()) > 1:
        click.echo(click.style("** Found {0} CRSs *****\n".format(len(crss.keys())), fg="red"))
    else:
        click.echo(click.style("** CRS *****\n", fg="blue"))

    for i, crs in enumerate(crss.keys()):
        click.echo(click.style("## CRS {0} ###### \n{1}".format(i+1, crs), fg="blue"))
        for raster in crss[crs]:
            click.echo(click.style(" {0}".format(raster), fg="green"))
        crs_dimensions = dimensions[crs]
        click.echo(click.style("\n Min/max extent values \n", fg="blue"))
        click.echo(click.style(" <Min left value> {0} in raster {1}".format(crs_dimensions['min_left']['value'],crs_dimensions['min_left']['raster']), fg="green"))
        click.echo(click.style(" <Min bottom value> {0} in raster {1}".format(crs_dimensions['min_bottom']['value'],crs_dimensions['min_bottom']['raster']), fg="green"))
        click.echo(click.style(" <Max right value> {0} in raster {1}".format(crs_dimensions['max_right']['value'],crs_dimensions['max_right']['raster']), fg="green"))
        click.echo(click.style(" <Max top value> {0} in raster {1}\n".format(crs_dimensions['max_top']['value'],crs_dimensions['max_top']['raster']), fg="green"))



@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.option('-f', '--filetype', default="*.tif", help='Raster file type')
@click.option('-v', '--verbose', is_flag=True, help='Verbose output')
def main(input_filepath, filetype, verbose):
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    logger.info('Harmonising rasters in {}'.format(input_filepath))
    # Get a list of raster files
    raster_files = find_files(input_filepath, filetype)
    process_rasters(raster_files, verbose=verbose)

if __name__ == '__main__':

    # not used in this stub but often useful for finding various files
    project_dir = os.path.join(os.path.dirname(__file__), os.pardir, os.pardir)

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    #load_dotenv(find_dotenv())

    main()

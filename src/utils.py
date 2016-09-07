#!/usr/bin/python3
# -*- coding: utf-8 -*-
""" Utility module for the workflow."""
import logging
import os
import pdb
import yaml
from colorlog import ColoredFormatter

# Classes ---------------------------------------------------------------------


class DataManifest(object):

    """ Class representing a data manifest file.

    Methods of this class can be used to access the content of the data
    manifest file in various ways.
    """

    def __init__(self, infile):
        if os.path.exists(infile):
            self._infile = infile
        else:
            raise IOError("Input data manifest file does not exist")

        # Read and parse the data manifest file
        self._data = self.parse_data_manifest()

    def parse_data_manifest(self):
        """ Parse datasets from a data manifest file.

        Current implementation can work with the following hierarchy:
            [URI]: dict
                "provider": str
                "collections": list
                    [NAME]: dict
                        "category": str
                        "metadata": list
                        "resources": list
        """

        def parse_collection(collection, collection_name):
            collection_data = {}
            collection_item = collection[collection_name]
            # Get the values for checking
            collection_keys = list(collection_item.keys())
            # Category is optional
            if 'category' in collection_keys:
                collection_data['collection_category'] = collection_item['category']
            # Metadata is optional
            if 'metadata' in collection_keys:
                collection_data['collection_metadata'] = collection_item['metadata']
            # Resources is required
            if 'resources' in collection_keys:
                collection_data['collection_resorces'] = collection_item['resources']
            else:
                raise ValueError("Collection {} contains no resources".format(collection_name))

            return collection_data

        items = []

        data_manifest = yaml.safe_load(open(self._infile, 'r'))
        for item in data_manifest:
            # Item should only have one key, which is the URI
            uri = list(item.keys())[0]
            # Loop over the provider content. A single provider can have
            # multiple collections.
            for provider_content in item[uri]:
                # Get the textual name of the provider
                provider_name = provider_content['provider']
                # Loop over the collections for this provider
                for collection in provider_content['collections']:
                    # Collection name is the only key
                    collection_name = list(collection.keys())[0]
                    # Check whether collection item is a dict or a list
                    # (indicating subgroups)
                    if isinstance(collection[collection_name], dict):
                        collection_data = parse_collection(collection,
                                                           collection_name)
                        collection_data['uri'] = uri
                        collection_data['provider_name'] = provider_name
                        collection_data['collection_name'] = collection_name
                        items.append(collection_data)
                    elif isinstance(collection[collection_name], list):
                        for subcollection in collection[collection_name]:
                            subcollection_name = list(subcollection.keys())[0]
                            subcollection_data = parse_collection(subcollection,
                                                                  subcollection_name)
                            subcollection_data['uri'] = uri
                            subcollection_data['provider_name'] = provider_name
                            subcollection_data['collection_name'] = collection_name
                            subcollection_data['subcollection_name'] = subcollection_name
                            items.append(subcollection_data)
                    else:
                        raise ValueError("Invalid collection type")

        return items

# Functions -------------------------------------------------------------------


def get_iteration_prefix(i, total):
    """ Return a String prefix for itarative task phases.

    :param i int current step.
    :param total int total steps.
    """
    return " [{0}/{1}]".format(i, total)


def get_local_logger(name, log_file=None, debug=False):
    """ Return a local logger."""
    date_format = "%Y-%m-%d %H:%M:%S"
    colFormatter = ColoredFormatter("%(log_color)s %(message)s%(reset)s",
                                    datefmt=date_format,
                                    reset=True,
                                    log_colors={'DEBUG':    'cyan',
                                                'INFO':     'green',
                                                'WARNING':  'yellow',
                                                'ERROR':    'red',
                                                'CRITICAL': 'red',
                                                })
    llogger = logging.getLogger(name)
    llogger.setLevel(logging.DEBUG)
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(colFormatter)
    if not debug:
        consoleHandler.setLevel(logging.INFO)
    llogger.addHandler(consoleHandler)

    if log_file is not None:
        fileFormatter = logging.Formatter("%(asctime)s [%(name)-10s] " +
                                          "[%(levelname)-5.5s] %(message)s",
                                          datefmt=date_format)

        fileHandler = logging.FileHandler(log_file, mode='w')

        fileHandler.setFormatter(fileFormatter)
        llogger.addHandler(fileHandler)

    return llogger


def pick_from_list(items, suffix):
    """ Pick an element from list ending with suffix.

    If no match is found, return None.

    :param items list of items to be searched.
    :suffix String suffix defining the match.
    """
    match = None
    for item in items:
        if item.endswith(suffix):
            match = item
    return match


def process_stdout(x, prefix=""):
    """Process stdout string returned by a shell command.

    If x is None, return an empty list.
    """
    if x is None or x == "":
        return []
    x = x.decode('utf-8')
    return [prefix + " " + item for item in x.split("\n") if item != ""]

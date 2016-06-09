#!/usr/bin/python3
# -*- coding: utf-8 -*-
""" Utility module for the workflow."""
import logging
import os
import yaml
from colorlog import ColoredFormatter


def get_iteration_prexix(i, total):
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


def parse_data_manifest(infile):
    """ Parse datasets from a data manifest file.

    Current implementation can work with the following hierarchy:
        uri:
            provider: str
            collections: list | dict

    If "collections" value is a list, then there are 3 hierarchy levels. If
    it is a dict, then there are 4 hierarchy levels. The key of the dict is a
    subcollection name and value is a list of datasets.

    :param infile String path to input YAML data manifest file.
    :return a dictionary of files with provider name as key and downloadable
            URIs as values
    """
    provider_datasets = {}
    if not os.path.exists(infile):
        raise IOError("Input data manifest file does not exist")

    data_manifest = yaml.safe_load(open(infile, 'r'))
    for item in data_manifest:
        # Item should only have one key, which is the URIq
        uri = list(item.keys())[0]
        # Loop over the collections
        for provider_collection in item[uri]:
            provider = provider_collection['provider']
            provider_datasets[provider] = []
            collection = provider_collection['collections']
            # Loop over collections and datasets
            for collection_name, datasets in collection.items():
                # 3 levels of hierarchy
                if isinstance(datasets, list):
                    for dataset in datasets:
                        dataset_uri = "/".join([uri, provider, collection_name,
                                               dataset])
                        provider_datasets[provider].append(dataset_uri)
                elif isinstance(datasets, dict):
                    for subcollection_name, subdatasets in datasets.items():
                        for subdataset in subdatasets:
                            dataset_uri = "/".join([uri, provider,
                                                    collection_name,
                                                    subcollection_name,
                                                    subdataset])
                            provider_datasets[provider].append(dataset_uri)
                else:
                    raise TypeError("Unknown collections type")

    return provider_datasets


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
    return [prefix + " " + item for item in x.split("\n")]

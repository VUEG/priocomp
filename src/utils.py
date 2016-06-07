#!/usr/bin/python3
# -*- coding: utf-8 -*-
import os
import yaml


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

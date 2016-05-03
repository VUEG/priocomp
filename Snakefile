import os
import requests
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

#################################################################################
# GLOBALS                                                                       #
#################################################################################

import dotenv
dotenv_path = '.env'
dotenv.load_dotenv(dotenv_path)


#################################################################################
# PROJECT RULES                                                                 #
#################################################################################

PROVIDE_DATASETS = ['agrobiodiversity_species_richness',
                    'carbon_sequestration',
                    'cultural_landscape_index_agro',
                    'cultural_landscape_index_forest',
                    'erosion_prevention',
                    'flood_regulation_supply',
                    'floodregulation',
                    'megafauna',
                    'nature_tourism',
                    'pollination_flows',
                    'pollination_visitprob',
                    'species_richness_farmland_birds_original1',
                    'species_richness_vascular_plants']

rule all:
    input:
        expand(["data/external/{dataset}/datapackage.json", "data/external/{dataset}/{dataset}.tif"], dataset=PROVIDE_DATASETS)

rule get_beehub_data:
    input:
        HTTP.remote(expand(["beehub.nl/environmental-geography-group/provide/{dataset}/datapackage.json",
                            "beehub.nl/environmental-geography-group/provide/{dataset}/{dataset}.tif"], dataset=PROVIDE_DATASETS),
                    username=os.environ.get("BEEHUB_USERNAME"), password=os.environ.get("BEEHUB_PASSWORD"), keep_local=False)
    output:
        expand(["data/external/{dataset}/datapackage.json", "data/external/{dataset}/{dataset}.tif"], dataset=PROVIDE_DATASETS)
    message:
        "Downloaded data from BeeHub"
    run:
        for i in range(0, len(input)):
            shell("mv {0} {1}".format(input[i], output[i]))

import requests
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

#################################################################################
# GLOBALS                                                                       #
#################################################################################


#################################################################################
# PROJECT RULES                                                                 #
#################################################################################

configfile: "data/data_manifest.yml"

rule get_beehub_data:
    input:
        HTTP.remote("beehub.nl/environmental-geography-group/provide/agrobiodiversity_species_richness/datapackage.json",
                    auth=requests.auth.HTTPDigestAuth("", ""), keep_local=True)
    run:
        shell("mv {input} datapackage.json")

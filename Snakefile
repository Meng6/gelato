from snakemake.utils import validate
configfile: "config.yaml"
validate(config, "tools/config.schema.yaml")
include: "rules/common.smk"
include: "rules/preprocessing.smk"
include: "rules/query.smk"

import itertools

files_to_compute = []
if len(config["GRAPHS"]) == 0:
    raise ValueError("Add at least one graph to GRAPHS in config.yaml. Remember to put the related data in data/external")

for graph in config["GRAPHS"]:
    if len(config[graph]["QUERIES"]["RUN"]) == 0:
        raise ValueError("Add at least one query to [{graph}][QUERIES][RUN] in config.yaml.".format(graph))
    if len(config[graph]["LANGUAGES_AND_TOOLS"]) == 0:
        raise ValueError("Add at least one language or tool to [{graph}][LANGUAGES_AND_TOOLS] in config.yaml.".format(graph))
    
    # Preprocess
    #files_to_compute.extend(expand("data/raw/{graph}/{lat}/{graph_data}.tsv", graph=graph.lower(), lat=map(str.lower, config[graph]["LANGUAGES_AND_TOOLS"]), graph_data=["advised", "dissertation", "person"]))

    # Query
    if "SEARCH_DESCENDANTS" in config[graph]["QUERIES"]["RUN"]:
        files_to_compute.extend(expand("data/query/mgdb/{lat}/output_search_descendants_for_{pid}.txt", lat=map(str.lower, config[graph]["LANGUAGES_AND_TOOLS"]), pid=config[graph]["QUERIES"]["SEARCH_DESCENDANTS"]["PIDS"]))
        files_to_compute.extend(expand("data/query/mgdb/{lat}/benchmark_search_descendants_for_{pid}.txt", lat=map(str.lower, config[graph]["LANGUAGES_AND_TOOLS"]), pid=config[graph]["QUERIES"]["SEARCH_DESCENDANTS"]["PIDS"]))

rule all:
    input:
        files_to_compute

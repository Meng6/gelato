from snakemake.utils import validate
configfile: "config.yaml"
validate(config, "tools/config.schema.yaml")
include: "rules/common.smk"
include: "rules/preprocessing.smk"
include: "rules/query.smk"
include: "rules/report.smk"

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
    if graph == "MGDB":
        files_to_compute.extend(expand("data/raw/{graph}/tsv/{graph_data}.tsv", graph=graph.lower(), graph_data=["advised", "dissertation", "person"]))
    if "CLINGO" in config[graph]["LANGUAGES_AND_TOOLS"]:
        files_to_compute.extend(expand("data/raw/{graph}/clingo/facts.tsv", graph=graph.lower()))
    if (graph == "MGDB") and ("BASHLOG" in config[graph]["LANGUAGES_AND_TOOLS"]):
        files_to_compute.extend(expand("data/raw/{graph}/bashlog/{graph_data}.tsv", graph=graph.lower(), graph_data=["advised", "dissertation", "person"]))

    # Query
    if "SEARCH_ANCESTORS" in config[graph]["QUERIES"]["RUN"]:
        if "BASHLOG" in config[graph]["LANGUAGES_AND_TOOLS"]:
            files_to_compute.extend(expand("data/query/{graph}/bashlog/search_ancestors_for_{pid}.dlog", graph=graph.lower(), pid=config[graph]["QUERIES"]["SEARCH_ANCESTORS"]["PIDS"]))
            files_to_compute.extend(expand("data/query/{graph}/bashlog/search_ancestors_for_{pid}.sh", graph=graph.lower(), pid=config[graph]["QUERIES"]["SEARCH_ANCESTORS"]["PIDS"]))
        files_to_compute.extend(expand("data/query/{graph}/{lat}/output_search_ancestors_for_{pid}.txt", graph=graph.lower(), lat=map(str.lower, config[graph]["LANGUAGES_AND_TOOLS"]), pid=config[graph]["QUERIES"]["SEARCH_ANCESTORS"]["PIDS"]))
        files_to_compute.extend(expand("data/query/{graph}/{lat}/benchmark_search_ancestors_for_{pid}.txt", graph=graph.lower(), lat=map(str.lower, config[graph]["LANGUAGES_AND_TOOLS"]), pid=config[graph]["QUERIES"]["SEARCH_ANCESTORS"]["PIDS"]))
    if "SEARCH_DESCENDANTS" in config[graph]["QUERIES"]["RUN"]:
        if "BASHLOG" in config[graph]["LANGUAGES_AND_TOOLS"]:
            files_to_compute.extend(expand("data/query/{graph}/bashlog/search_descendants_for_{pid}.dlog", graph=graph.lower(), pid=config[graph]["QUERIES"]["SEARCH_DESCENDANTS"]["PIDS"]))
            files_to_compute.extend(expand("data/query/{graph}/bashlog/search_descendants_for_{pid}.sh", graph=graph.lower(), pid=config[graph]["QUERIES"]["SEARCH_DESCENDANTS"]["PIDS"]))
        files_to_compute.extend(expand("data/query/{graph}/{lat}/output_search_descendants_for_{pid}.txt", graph=graph.lower(), lat=map(str.lower, config[graph]["LANGUAGES_AND_TOOLS"]), pid=config[graph]["QUERIES"]["SEARCH_DESCENDANTS"]["PIDS"]))
        files_to_compute.extend(expand("data/query/{graph}/{lat}/benchmark_search_descendants_for_{pid}.txt", graph=graph.lower(), lat=map(str.lower, config[graph]["LANGUAGES_AND_TOOLS"]), pid=config[graph]["QUERIES"]["SEARCH_DESCENDANTS"]["PIDS"]))

    # Report
    if config["BENCHMARK"]["REPORT"]:
        files_to_compute.extend(expand("data/report/{graph}_benchmarks.csv", graph=graph.lower()))

rule all:
    input:
        files_to_compute

import glob

def get_mgdb_external_data(wildcards):
    for file_path in glob.glob(config["MGDB"]["DATA_SOURCE"]["FOLDER_PATH"] + "/*"):
        file_name = file_path.split("/")[-1].split(".")[0]
        if file_name == wildcards.graph_data:
            return file_path
    raise ValueError("MGDB: external file {graph_data} does not exists".format(wildcards.graph_data))

def get_mgdb_merge_benchmarks_input(wildcards):
    input = []
    if "SEARCH_DESCENDANTS" in config["MGDB"]["QUERIES"]["RUN"]:
        input.extend(expand("data/query/mgdb/{lat}/benchmark_search_descendants_for_{pid}.txt", lat=map(str.lower, config[graph]["LANGUAGES_AND_TOOLS"]), pid=config[graph]["QUERIES"]["SEARCH_DESCENDANTS"]["PIDS"]))
    if "SEARCH_ANCESTORS" in config["MGDB"]["QUERIES"]["RUN"]:
        input.extend(expand("data/query/mgdb/{lat}/benchmark_search_ancestors_for_{pid}.txt", lat=map(str.lower, config[graph]["LANGUAGES_AND_TOOLS"]), pid=config[graph]["QUERIES"]["SEARCH_ANCESTORS"]["PIDS"]))
    return input
    
import glob

def get_mgdb_external_data(wildcards):
    for file_path in glob.glob(config["MGDB"]["DATA_SOURCE"]["FOLDER_PATH"] + "/*"):
        file_name = file_path.split("/")[-1].split(".")[0]
        if file_name == wildcards.graph_data:
            return file_path
    raise ValueError("MGDB: external file {graph_data} does not exists".format(wildcards.graph_data))

def get_mgdb_merge_benchmarks_input(wildcards):
    input = []
    for query in config["MGDB"]["QUERIES"]["RUN"]:
        if query in ["UNARY_SEARCH_ANCESTORS", "BINARY_SEARCH_ANCESTORS", "UNARY_SEARCH_DESCENDANTS", "BINARY_SEARCH_DESCENDANTS"]:
            input.extend(expand("data/query/mgdb/{lat}/benchmark_{query}_for_{pid}.txt", lat=map(str.lower, config[graph]["LANGUAGES_AND_TOOLS"]), query=query.lower(), pid=config[graph]["QUERIES"][query]["PIDS"]))
    return input

def optional_mgdb_cypher_input(wildcards):
    if config["MGDB"]["DATA_SOURCE"]["CYPHER"]["LOAD_DATA"]:
        return "data/raw/mgdb/cypher/load_mgdb_data_to_neo4j.done"
    return []

def optional_mgdb_sparql_input(wildcards):
    if config["MGDB"]["DATA_SOURCE"]["SPARQL"]["LOAD_DATA"]:
        return "data/raw/mgdb/sparql/load_mgdb_data_to_blazegraph.done"
    return []

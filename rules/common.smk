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
            input.extend(expand("data/query/mgdb/{lat}/benchmark_{query}_for_{pid}.txt", lat=map(str.lower, config["MGDB"]["LANGUAGES_AND_TOOLS"]), query=query.lower(), pid=config["MGDB"]["QUERIES"][query]["PIDS"]))
        if query == "LOWEST_COMMON_ANCESTORS":
            input.extend(expand("data/query/mgdb/{lat}/benchmark_{query}_of_{pid1}_and_{pid2}.txt", lat=map(str.lower, config["MGDB"]["LANGUAGES_AND_TOOLS"]), query=query.lower(), pid1=config["MGDB"]["QUERIES"][query]["PID1"], pid2=config["MGDB"]["QUERIES"][query]["PID2"]))
    return input

def optional_mgdb_cypher_input(wildcards):
    if config["MGDB"]["DATA_SOURCE"]["CYPHER"]["LOAD_DATA"]:
        return "data/raw/mgdb/cypher/load_mgdb_data_to_neo4j.done"
    return []

def optional_mgdb_sparql_input(wildcards):
    if config["MGDB"]["DATA_SOURCE"]["SPARQL"]["LOAD_DATA"]:
        return "data/raw/mgdb/sparql/load_mgdb_data_to_blazegraph.done"
    return []

def get_mgdb_prepare_bashscript_for_bashlog_2_input(wildcards):
    if wildcards.query.startswith("interim"):
        input = ["data/raw/mgdb/bashlog/advised.tsv", "data/raw/mgdb/bashlog/dissertation.tsv", "data/raw/mgdb/bashlog/person.tsv"]
    else:
        input = ["data/query/mgdb/bashlog/interim_advise.tsv", "data/query/mgdb/bashlog/interim_common_ancestors_of_{pid1}_and_{pid2}.tsv"]
    return input

def get_fish_merge_benchmarks_input(wildcards):
    input = []
    for query in config["FISH"]["QUERIES"]["RUN"]:
        if query in ["UNARY_SEARCH_ANCESTORS", "BINARY_SEARCH_ANCESTORS", "UNARY_SEARCH_DESCENDANTS", "BINARY_SEARCH_DESCENDANTS"]:
            input.extend(expand("data/query/fish/{max_hamming_number}/{lat}/benchmark_{query}_for_{pid}.txt", max_hamming_number=wildcards.max_hamming_number, lat=map(str.lower, config["FISH"]["LANGUAGES_AND_TOOLS"]), query=query.lower(), pid=config["FISH"]["QUERIES"][query]["PIDS"]))
        if query == "LOWEST_COMMON_ANCESTORS":
            input.extend(expand("data/query/fish/{max_hamming_number}/{lat}/benchmark_{query}_of_{pid1}_and_{pid2}.txt", max_hamming_number=wildcards.max_hamming_number, lat=map(str.lower, config["FISH"]["LANGUAGES_AND_TOOLS"]), query=query.lower(), pid1=config["FISH"]["QUERIES"][query]["PID1"], pid2=config["FISH"]["QUERIES"][query]["PID2"]))
    return input

def get_sail_merge_benchmarks_input(wildcards):
    input = []
    for query in config["SAIL"]["QUERIES"]["RUN"]:
        if query in ["UNARY_SEARCH_ANCESTORS", "BINARY_SEARCH_ANCESTORS", "UNARY_SEARCH_DESCENDANTS", "BINARY_SEARCH_DESCENDANTS"]:
            input.extend(expand("data/query/sail/{max_hamming_number}/{lat}/benchmark_{query}_for_{pid}.txt", max_hamming_number=wildcards.max_hamming_number, lat=map(str.lower, config["SAIL"]["LANGUAGES_AND_TOOLS"]), query=query.lower(), pid=config["SAIL"]["QUERIES"][query]["PIDS"]))
        if query == "LOWEST_COMMON_ANCESTORS":
            input.extend(expand("data/query/sail/{max_hamming_number}/{lat}/benchmark_{query}_of_{pid1}_and_{pid2}.txt", max_hamming_number=wildcards.max_hamming_number, lat=map(str.lower, config["SAIL"]["LANGUAGES_AND_TOOLS"]), query=query.lower(), pid1=config["SAIL"]["QUERIES"][query]["PID1"], pid2=config["SAIL"]["QUERIES"][query]["PID2"]))
    return input

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
    
    if graph == "MGDB":
        # Preprocess
        files_to_compute.extend(expand("data/raw/{graph}/tsv/{graph_data}.tsv", graph=graph.lower(), graph_data=["advised", "dissertation", "person"]))
        if "CLINGO" in config[graph]["LANGUAGES_AND_TOOLS"]:
            files_to_compute.extend(expand("data/raw/{graph}/clingo/facts.tsv", graph=graph.lower()))
        if "BASHLOG" in config[graph]["LANGUAGES_AND_TOOLS"]:
            files_to_compute.extend(expand("data/raw/{graph}/bashlog/{graph_data}.tsv", graph=graph.lower(), graph_data=["advised", "dissertation", "person"]))
        # Load data into SQLite or Neo4j or Blazegraph DB
        files_to_compute.extend(expand("data/external/{graph}/{graph}.db", graph=graph.lower()))        
        if config[graph]["DATA_SOURCE"]["CYPHER"]["LOAD_DATA"]:
            files_to_compute.extend(expand("data/raw/{graph}/cypher/load_{graph}_data_to_neo4j.done", graph=graph.lower()))
        if config[graph]["DATA_SOURCE"]["SPARQL"]["LOAD_DATA"]:
            files_to_compute.extend(expand("data/raw/{graph}/sparql/{graph}-ttl.txt", graph=graph.lower()))
            files_to_compute.extend(expand("data/raw/{graph}/sparql/load_{graph}_data_to_blazegraph.done", graph=graph.lower()))
        
        # Query
        for query in config[graph]["QUERIES"]["RUN"]:
            if query in ["UNARY_SEARCH_ANCESTORS", "BINARY_SEARCH_ANCESTORS", "UNARY_SEARCH_DESCENDANTS", "BINARY_SEARCH_DESCENDANTS"]:
                files_to_compute.extend(expand("data/query/{graph}/{lat}/output_{query}_for_{pid}.csv", graph=graph.lower(), lat=map(str.lower, config[graph]["LANGUAGES_AND_TOOLS"]), query=query.lower(), pid=config[graph]["QUERIES"][query]["PIDS"]))
                files_to_compute.extend(expand("data/query/{graph}/{lat}/benchmark_{query}_for_{pid}.txt", graph=graph.lower(), lat=map(str.lower, config[graph]["LANGUAGES_AND_TOOLS"]), query=query.lower(), pid=config[graph]["QUERIES"][query]["PIDS"]))
            if query in ["LOWEST_COMMON_ANCESTORS", "LOWEST_COMMON_ANCESTORS_PATH"]:
                files_to_compute.extend(expand("data/query/{graph}/{lat}/output_{query}_of_{pid1}_and_{pid2}.csv", graph=graph.lower(), lat=map(str.lower, config[graph]["LANGUAGES_AND_TOOLS"]), query=query.lower(), pid1=config[graph]["QUERIES"][query]["PID1"], pid2=config[graph]["QUERIES"][query]["PID2"]))
                files_to_compute.extend(expand("data/query/{graph}/{lat}/benchmark_{query}_of_{pid1}_and_{pid2}.txt", graph=graph.lower(), lat=map(str.lower, config[graph]["LANGUAGES_AND_TOOLS"]), query=query.lower(), pid1=config[graph]["QUERIES"][query]["PID1"], pid2=config[graph]["QUERIES"][query]["PID2"]))

        # Report
        if config["BENCHMARK"]["REPORT"]:
            files_to_compute.extend(expand("data/report/{graph}/benchmarks.csv", graph=graph.lower()))
        if config["OUTPUT_PER_QUERY"]["REPORT"]:
            for query in config[graph]["QUERIES"]["RUN"]:
                if query in ["UNARY_SEARCH_ANCESTORS", "BINARY_SEARCH_ANCESTORS", "UNARY_SEARCH_DESCENDANTS", "BINARY_SEARCH_DESCENDANTS"]:
                    files_to_compute.extend(expand("data/report/{graph}/report_{query}_for_{pid}_{lat}.html", graph=graph.lower(), lat=map(str.lower, config[graph]["LANGUAGES_AND_TOOLS"]), query=query.lower(), pid=config[graph]["QUERIES"][query]["PIDS"]))
                if query in ["LOWEST_COMMON_ANCESTORS", "LOWEST_COMMON_ANCESTORS_PATH"]:
                    files_to_compute.extend(expand("data/report/{graph}/report_{query}_of_{pid1}_and_{pid2}_{lat}.html", graph=graph.lower(), lat=map(str.lower, config[graph]["LANGUAGES_AND_TOOLS"]), query=query.lower(), pid1=config[graph]["QUERIES"][query]["PID1"], pid2=config[graph]["QUERIES"][query]["PID2"]))


    if graph == "FISH" or graph == "SAIL":
        # Preprocess
        files_to_compute.extend(expand("data/raw/{graph}/tsv/{graph}_{max_hamming_number}.tsv", graph=graph.lower(), max_hamming_number=config[graph]["DATA_SOURCE"]["MAX_HAMMING_NUMBER"]))
        if "CLINGO" in config[graph]["LANGUAGES_AND_TOOLS"]:
            files_to_compute.extend(expand("data/raw/{graph}/clingo/facts_{max_hamming_number}.tsv", graph=graph.lower(), max_hamming_number=config[graph]["DATA_SOURCE"]["MAX_HAMMING_NUMBER"]))
        if "BASHLOG" in config[graph]["LANGUAGES_AND_TOOLS"]:
            files_to_compute.extend(expand("data/raw/{graph}/bashlog/{graph}_{max_hamming_number}.tsv", graph=graph.lower(), max_hamming_number=config[graph]["DATA_SOURCE"]["MAX_HAMMING_NUMBER"]))
        # Load data into SQLite or Neo4j or Blazegraph DB
        files_to_compute.extend(expand("data/external/{graph}/{graph}_{max_hamming_number}.db", graph=graph.lower(), max_hamming_number=config[graph]["DATA_SOURCE"]["MAX_HAMMING_NUMBER"]))
    
        # Query
        for query in config[graph]["QUERIES"]["RUN"]:
            if query == "LOWEST_COMMON_ANCESTORS":
                files_to_compute.extend(expand("data/query/{graph}/{max_hamming_number}/{lat}/output_{query}_of_{pid1}_and_{pid2}.csv", graph=graph.lower(), max_hamming_number=config[graph]["DATA_SOURCE"]["MAX_HAMMING_NUMBER"], lat=map(str.lower, config[graph]["LANGUAGES_AND_TOOLS"]), query=query.lower(), pid1=config[graph]["QUERIES"][query]["PID1"], pid2=config[graph]["QUERIES"][query]["PID2"]))
                files_to_compute.extend(expand("data/query/{graph}/{max_hamming_number}/{lat}/benchmark_{query}_of_{pid1}_and_{pid2}.txt", graph=graph.lower(), max_hamming_number=config[graph]["DATA_SOURCE"]["MAX_HAMMING_NUMBER"], lat=map(str.lower, config[graph]["LANGUAGES_AND_TOOLS"]), query=query.lower(), pid1=config[graph]["QUERIES"][query]["PID1"], pid2=config[graph]["QUERIES"][query]["PID2"]))

        # Report
        if config["BENCHMARK"]["REPORT"]:
            files_to_compute.extend(expand("data/report/{graph}/{max_hamming_number}/benchmarks.csv", graph=graph.lower(), max_hamming_number=config[graph]["DATA_SOURCE"]["MAX_HAMMING_NUMBER"]))
        if config["OUTPUT_PER_QUERY"]["REPORT"]:
            for query in config[graph]["QUERIES"]["RUN"]:
                if query in ["UNARY_SEARCH_ANCESTORS", "BINARY_SEARCH_ANCESTORS", "UNARY_SEARCH_DESCENDANTS", "BINARY_SEARCH_DESCENDANTS"]:
                    files_to_compute.extend(expand("data/report/{graph}/{max_hamming_number}/report_{query}_for_{pid}_{lat}.html", graph=graph.lower(), max_hamming_number=config[graph]["DATA_SOURCE"]["MAX_HAMMING_NUMBER"], lat=map(str.lower, config[graph]["LANGUAGES_AND_TOOLS"]), query=query.lower(), pid=config[graph]["QUERIES"][query]["PIDS"]))
                if query == "LOWEST_COMMON_ANCESTORS":
                    files_to_compute.extend(expand("data/report/{graph}/{max_hamming_number}/report_{query}_of_{pid1}_and_{pid2}_{lat}.html", graph=graph.lower(), max_hamming_number=config[graph]["DATA_SOURCE"]["MAX_HAMMING_NUMBER"], lat=map(str.lower, config[graph]["LANGUAGES_AND_TOOLS"]), query=query.lower(), pid1=config[graph]["QUERIES"][query]["PID1"], pid2=config[graph]["QUERIES"][query]["PID2"]))

# Stats of the graphs
files_to_compute.append("data/report/overall_stats_of_graphs.csv")

rule all:
    input:
        files_to_compute

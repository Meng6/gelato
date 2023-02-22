# BENCHMARK
rule mgdb_stats:
    input:
        advised = "data/raw/mgdb/tsv/advised.tsv",
        dissertation = "data/raw/mgdb/tsv/dissertation.tsv",
        person = "data/raw/mgdb/tsv/person.tsv"
    output:
        "data/report/mgdb/stats.csv"
    script:
        "../src/report/mgdb_stats.py"

rule fish_stats:
    input:
        expand("data/raw/fish/tsv/fish_{max_hamming_number}.tsv", max_hamming_number=config["FISH"]["DATA_SOURCE"]["MAX_HAMMING_NUMBER"])
    output:
        "data/report/fish/stats.csv"
    script:
        "../src/report/fish_stats.py"

rule sail_stats:
    input:
        expand("data/raw/sail/tsv/sail_{max_hamming_number}.tsv", max_hamming_number=config["SAIL"]["DATA_SOURCE"]["MAX_HAMMING_NUMBER"])
    output:
        "data/report/sail/stats.csv"
    script:
        "../src/report/sail_stats.py"

rule merge_stats_of_graphs:
    input:
        expand("data/report/{graph}/stats.csv", graph=map(str.lower, config["GRAPHS"]))
    params:
        graphs = config["GRAPHS"]
    output:
        "data/report/overall_stats_of_graphs.csv"
    script:
        "../src/report/merge_stats_of_graphs.py"

rule mgdb_merge_benchmarks:
    input:
        get_mgdb_merge_benchmarks_input
    output:
        "data/report/mgdb/benchmarks.csv"
    script:
        "../src/report/merge_benchmarks.py"

rule fish_merge_benchmarks:
    input:
        get_fish_merge_benchmarks_input
    params:
        max_hamming_number = "{max_hamming_number}"
    output:
        "data/report/fish/{max_hamming_number}/benchmarks.csv"
    script:
        "../src/report/merge_benchmarks.py"

rule sail_merge_benchmarks:
    input:
        get_sail_merge_benchmarks_input
    params:
        max_hamming_number = "{max_hamming_number}"
    output:
        "data/report/sail/{max_hamming_number}/benchmarks.csv"
    script:
        "../src/report/merge_benchmarks.py"

# OUTPUT_PER_QUERY
# Search ancestors
rule report_mgdb_unary_search_ancestors:
    input:
        formated_data = "data/query/mgdb/{lat}/output_unary_search_ancestors_for_{pid}.txt"
    params:
        query = "unary_search_ancestors",
        edge = [],
        include_interactive_table = config["OUTPUT_PER_QUERY"]["INCLUDE_INTERACTIVE_TABLE"],
        include_network_graph = config["OUTPUT_PER_QUERY"]["INCLUDE_NETWORK_GRAPH"]
    output:
        "data/report/mgdb/report_unary_search_ancestors_for_{pid}_{lat}.html"
    script:
        "../src/report/report_output_per_query.py"

rule report_mgdb_binary_search_ancestors:
    input:
        formated_data = "data/query/mgdb/{lat}/output_binary_search_ancestors_for_{pid}.txt"
    params:
        query = "binary_search_ancestors",
        edge = [["student_name", "advisor_name"]],
        include_interactive_table = config["OUTPUT_PER_QUERY"]["INCLUDE_INTERACTIVE_TABLE"],
        include_network_graph = config["OUTPUT_PER_QUERY"]["INCLUDE_NETWORK_GRAPH"]
    output:
        "data/report/mgdb/report_binary_search_ancestors_for_{pid}_{lat}.html"
    script:
        "../src/report/report_output_per_query.py"

# Search descendants
rule report_mgdb_unary_search_descendants:
    input:
        formated_data = "data/query/mgdb/{lat}/output_unary_search_descendants_for_{pid}.txt"
    params:
        query = "unary_search_descendants",
        edge = [],
        include_interactive_table = config["OUTPUT_PER_QUERY"]["INCLUDE_INTERACTIVE_TABLE"],
        include_network_graph = config["OUTPUT_PER_QUERY"]["INCLUDE_NETWORK_GRAPH"]
    output:
        "data/report/mgdb/report_unary_search_descendants_for_{pid}_{lat}.html"
    script:
        "../src/report/report_output_per_query.py"

rule report_mgdb_binary_search_descendants:
    input:
        formated_data = "data/query/mgdb/{lat}/output_binary_search_descendants_for_{pid}.txt"
    params:
        query = "binary_search_descendants",
        edge = [["student_name", "advisor_name"]],
        include_interactive_table = config["OUTPUT_PER_QUERY"]["INCLUDE_INTERACTIVE_TABLE"],
        include_network_graph = config["OUTPUT_PER_QUERY"]["INCLUDE_NETWORK_GRAPH"]
    output:
        "data/report/mgdb/report_binary_search_descendants_for_{pid}_{lat}.html"
    script:
        "../src/report/report_output_per_query.py"

# Lowest common ancestors
rule report_mgdb_lowest_common_ancestors:
    input:
        formated_data = "data/query/mgdb/{lat}/output_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt"
    params:
        query = "lowest_common_ancestors",
        edge = [],
        include_interactive_table = config["OUTPUT_PER_QUERY"]["INCLUDE_INTERACTIVE_TABLE"],
        include_network_graph = config["OUTPUT_PER_QUERY"]["INCLUDE_NETWORK_GRAPH"]
    output:
        "data/report/mgdb/report_lowest_common_ancestors_of_{pid1}_and_{pid2}_{lat}.html"
    script:
        "../src/report/report_output_per_query.py"

rule report_fish_lowest_common_ancestors:
    input:
        formated_data = "data/query/fish/{max_hamming_number}/{lat}/output_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt"
    params:
        query = "lowest_common_ancestors",
        edge = [],
        include_interactive_table = config["OUTPUT_PER_QUERY"]["INCLUDE_INTERACTIVE_TABLE"],
        include_network_graph = config["OUTPUT_PER_QUERY"]["INCLUDE_NETWORK_GRAPH"]
    output:
        "data/report/fish/{max_hamming_number}/report_lowest_common_ancestors_of_{pid1}_and_{pid2}_{lat}.html"
    script:
        "../src/report/report_output_per_query.py"

rule report_sail_lowest_common_ancestors:
    input:
        formated_data = "data/query/sail/{max_hamming_number}/{lat}/output_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt"
    params:
        query = "lowest_common_ancestors",
        edge = [],
        include_interactive_table = config["OUTPUT_PER_QUERY"]["INCLUDE_INTERACTIVE_TABLE"],
        include_network_graph = config["OUTPUT_PER_QUERY"]["INCLUDE_NETWORK_GRAPH"]
    output:
        "data/report/sail/{max_hamming_number}/report_lowest_common_ancestors_of_{pid1}_and_{pid2}_{lat}.html"
    script:
        "../src/report/report_output_per_query.py"
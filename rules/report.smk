rule mgdb_stats:
    input:
        advised = "data/raw/mgdb/tsv/advised.tsv",
        dissertation = "data/raw/mgdb/tsv/dissertation.tsv",
        person = "data/raw/mgdb/tsv/person.tsv"
    output:
        "data/report/mgdb_stats.csv"
    script:
        "../src/report/mgdb_stats.py"

rule fish_stats:
    input:
        expand("data/raw/fish/tsv/fish_{max_hamming_number}.tsv", max_hamming_number=config["FISH"]["DATA_SOURCE"]["MAX_HAMMING_NUMBER"])
    output:
        "data/report/fish_stats.csv"
    script:
        "../src/report/fish_stats.py"

rule sail_stats:
    input:
        expand("data/raw/sail/tsv/sail_{max_hamming_number}.tsv", max_hamming_number=config["SAIL"]["DATA_SOURCE"]["MAX_HAMMING_NUMBER"])
    output:
        "data/report/sail_stats.csv"
    script:
        "../src/report/sail_stats.py"

rule merge_stats_of_graphs:
    input:
        expand("data/report/{graph}_stats.csv", graph=map(str.lower, config["GRAPHS"]))
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
        "data/report/mgdb_benchmarks.csv"
    script:
        "../src/report/merge_benchmarks.py"

rule fish_merge_benchmarks:
    input:
        get_fish_merge_benchmarks_input
    params:
        max_hamming_number = "{max_hamming_number}"
    output:
        "data/report/fish_{max_hamming_number}_benchmarks.csv"
    script:
        "../src/report/merge_benchmarks.py"

rule sail_merge_benchmarks:
    input:
        get_sail_merge_benchmarks_input
    params:
        max_hamming_number = "{max_hamming_number}"
    output:
        "data/report/sail_{max_hamming_number}_benchmarks.csv"
    script:
        "../src/report/merge_benchmarks.py"
rule mgdb_merge_benchmarks:
    input:
        get_mgdb_merge_benchmarks_input
    output:
        "data/report/mgdb_benchmarks.csv"
    script:
        "../src/report/merge_benchmarks.py"
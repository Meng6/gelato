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
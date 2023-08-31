# Query
# Search ancestors
rule mgdb_unary_search_ancestors_with_python:
    input:
        advised = "data/raw/mgdb/tsv/advised.tsv",
        dissertation = "data/raw/mgdb/tsv/dissertation.tsv",
        person = "data/raw/mgdb/tsv/person.tsv"
    params:
        pid = "{pid}",
        lat = "python",
        query = "unary_search_ancestors"
    output:
        output = "data/query/mgdb/python/output_unary_search_ancestors_for_{pid}.csv"
    log:
        "data/query/mgdb/python/unary_search_ancestors_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/python/benchmark_unary_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_binary_search_ancestors_with_python:
    input:
        advised = "data/raw/mgdb/tsv/advised.tsv",
        dissertation = "data/raw/mgdb/tsv/dissertation.tsv",
        person = "data/raw/mgdb/tsv/person.tsv"
    params:
        pid = "{pid}",
        lat = "python",
        query = "binary_search_ancestors"
    output:
        "data/query/mgdb/python/output_binary_search_ancestors_for_{pid}.csv"
    log:
        "data/query/mgdb/python/binary_search_ancestors_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/python/benchmark_binary_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_unary_search_ancestors_with_sql:
    input:
        database = "data/external/mgdb/mgdb.db"
    params:
        pid = "{pid}",
        lat = "sql",
        query = "unary_search_ancestors"
    output:
        "data/query/mgdb/sql/output_unary_search_ancestors_for_{pid}.csv"
    log:
        "data/query/mgdb/sql/unary_search_ancestors_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/sql/benchmark_unary_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_binary_search_ancestors_with_sql:
    input:
        database = "data/external/mgdb/mgdb.db"
    params:
        pid = "{pid}",
        lat = "sql",
        query = "binary_search_ancestors"
    output:
        "data/query/mgdb/sql/output_binary_search_ancestors_for_{pid}.csv"
    log:
        "data/query/mgdb/sql/binary_search_ancestors_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/sql/benchmark_binary_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_unary_search_ancestors_with_cypher:
    input:
        optional_mgdb_cypher_input
    params:
        pid = "{pid}",
        database_group = config["MGDB"]["DATA_SOURCE"]["CYPHER"]["DATABASE_GROUP"],
        lat = "cypher",
        query = "unary_search_ancestors"
    output:
        "data/query/mgdb/cypher/output_unary_search_ancestors_for_{pid}.csv"
    log:
        "data/query/mgdb/cypher/unary_search_ancestors_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/cypher/benchmark_unary_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_binary_search_ancestors_with_cypher:
    input:
        optional_mgdb_cypher_input
    params:
        pid = "{pid}",
        database_group = config["MGDB"]["DATA_SOURCE"]["CYPHER"]["DATABASE_GROUP"],
        lat = "cypher",
        query = "binary_search_ancestors"
    output:
        "data/query/mgdb/cypher/output_binary_search_ancestors_for_{pid}.csv"
    log:
        "data/query/mgdb/cypher/binary_search_ancestors_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/cypher/benchmark_binary_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_unary_search_ancestors_with_sparql:
    input:
        optional_mgdb_sparql_input
    params:
        pid = "{pid}",
        namespace="MGDB",
        database_group = config["MGDB"]["DATA_SOURCE"]["SPARQL"]["DATABASE_GROUP"],
        lat = "sparql",
        query = "unary_search_ancestors"
    output:
        "data/query/mgdb/sparql/output_unary_search_ancestors_for_{pid}.csv"
    log:
        "data/query/mgdb/sparql/unary_search_ancestors_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/sparql/benchmark_unary_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_binary_search_ancestors_with_sparql:
    input:
        optional_mgdb_sparql_input
    params:
        pid = "{pid}",
        namespace="MGDB",
        database_group = config["MGDB"]["DATA_SOURCE"]["SPARQL"]["DATABASE_GROUP"],
        lat = "sparql",
        query = "binary_search_ancestors"
    output:
        "data/query/mgdb/sparql/output_binary_search_ancestors_for_{pid}.csv"
    log:
        "data/query/mgdb/sparql/binary_search_ancestors_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/sparql/benchmark_binary_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_unary_search_ancestors_with_clingo:
    input:
        facts = "data/raw/mgdb/clingo/facts.tsv"
    params:
        pid = "{pid}",
        lat = "clingo",
        query = "unary_search_ancestors"
    output:
        "data/query/mgdb/clingo/output_unary_search_ancestors_for_{pid}.csv"
    log:
        "data/query/mgdb/clingo/unary_search_ancestors_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/clingo/benchmark_unary_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_binary_search_ancestors_with_clingo:
    input:
        facts = "data/raw/mgdb/clingo/facts.tsv"
    params:
        pid = "{pid}",
        lat = "clingo",
        query = "binary_search_ancestors"
    output:
        "data/query/mgdb/clingo/output_binary_search_ancestors_for_{pid}.csv"
    log:
        "data/query/mgdb/clingo/binary_search_ancestors_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/clingo/benchmark_binary_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_unary_search_ancestors_with_bashlog:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv"
    params:
        pid = "{pid}",
        lat = "bashlog",
        query = "unary_search_ancestors"
    output:
        "data/query/mgdb/bashlog/output_unary_search_ancestors_for_{pid}.csv"
    log:
        "data/query/mgdb/bashlog/unary_search_ancestors_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/bashlog/benchmark_unary_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_binary_search_ancestors_with_bashlog:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv"
    params:
        pid = "{pid}",
        lat = "bashlog",
        query = "binary_search_ancestors"
    output:
        "data/query/mgdb/bashlog/output_binary_search_ancestors_for_{pid}.csv"
    log:
        "data/query/mgdb/bashlog/binary_search_ancestors_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/bashlog/benchmark_binary_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_unary_search_ancestors_with_logica:
    input:
        database = "data/external/mgdb/mgdb.db"
    params:
        pid = "{pid}",
        lat = "logica",
        query = "unary_search_ancestors"
    output:
        "data/query/mgdb/logica/output_unary_search_ancestors_for_{pid}.csv"
    log:
        "data/query/mgdb/logica/unary_search_ancestors_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/logica/benchmark_unary_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_binary_search_ancestors_with_logica:
    input:
        database = "data/external/mgdb/mgdb.db"
    params:
        pid = "{pid}",
        lat = "logica",
        query = "binary_search_ancestors"
    output:
        "data/query/mgdb/logica/output_binary_search_ancestors_for_{pid}.csv"
    log:
        "data/query/mgdb/logica/binary_search_ancestors_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/logica/benchmark_binary_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

# Search descendants
rule mgdb_unary_search_descendants_with_python:
    input:
        advised = "data/raw/mgdb/tsv/advised.tsv",
        dissertation = "data/raw/mgdb/tsv/dissertation.tsv",
        person = "data/raw/mgdb/tsv/person.tsv"
    params:
        pid = "{pid}",
        lat = "python",
        query = "unary_search_descendants"
    output:
        "data/query/mgdb/python/output_unary_search_descendants_for_{pid}.csv"
    log:
        "data/query/mgdb/python/unary_search_descendants_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/python/benchmark_unary_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_binary_search_descendants_with_python:
    input:
        advised = "data/raw/mgdb/tsv/advised.tsv",
        dissertation = "data/raw/mgdb/tsv/dissertation.tsv",
        person = "data/raw/mgdb/tsv/person.tsv"
    params:
        pid = "{pid}",
        lat = "python",
        query = "binary_search_descendants"
    output:
        "data/query/mgdb/python/output_binary_search_descendants_for_{pid}.csv"
    log:
        "data/query/mgdb/python/binary_search_descendants_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/python/benchmark_binary_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_unary_search_descendants_with_sql:
    input:
        database = "data/external/mgdb/mgdb.db"
    params:
        pid = "{pid}",
        lat = "sql",
        query = "unary_search_descendants"
    output:
        "data/query/mgdb/sql/output_unary_search_descendants_for_{pid}.csv"
    log:
        "data/query/mgdb/sql/unary_search_descendants_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/sql/benchmark_unary_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_binary_search_descendants_with_sql:
    input:
        database = "data/external/mgdb/mgdb.db"
    params:
        pid = "{pid}",
        lat = "sql",
        query = "binary_search_descendants"
    output:
        "data/query/mgdb/sql/output_binary_search_descendants_for_{pid}.csv"
    log:
        "data/query/mgdb/sql/binary_search_descendants_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/sql/benchmark_binary_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_unary_search_descendants_with_cypher:
    input:
        optional_mgdb_cypher_input
    params:
        pid = "{pid}",
        database_group = config["MGDB"]["DATA_SOURCE"]["CYPHER"]["DATABASE_GROUP"],
        lat = "cypher",
        query = "unary_search_descendants"
    output:
        "data/query/mgdb/cypher/output_unary_search_descendants_for_{pid}.csv"
    log:
        "data/query/mgdb/cypher/unary_search_descendants_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/cypher/benchmark_unary_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_binary_search_descendants_with_cypher:
    input:
        optional_mgdb_cypher_input
    params:
        pid = "{pid}",
        database_group = config["MGDB"]["DATA_SOURCE"]["CYPHER"]["DATABASE_GROUP"],
        lat = "cypher",
        query = "binary_search_descendants"
    output:
        "data/query/mgdb/cypher/output_binary_search_descendants_for_{pid}.csv"
    log:
        "data/query/mgdb/cypher/binary_search_descendants_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/cypher/benchmark_binary_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_unary_search_descendants_with_sparql:
    input:
        optional_mgdb_sparql_input
    params:
        pid = "{pid}",
        namespace="MGDB",
        database_group = config["MGDB"]["DATA_SOURCE"]["SPARQL"]["DATABASE_GROUP"],
        lat = "sparql",
        query = "unary_search_descendants"
    output:
        "data/query/mgdb/sparql/output_unary_search_descendants_for_{pid}.csv"
    log:
        "data/query/mgdb/sparql/unary_search_descendants_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/sparql/benchmark_unary_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_binary_search_descendants_with_sparql:
    input:
        optional_mgdb_sparql_input
    params:
        pid = "{pid}",
        namespace="MGDB",
        database_group = config["MGDB"]["DATA_SOURCE"]["SPARQL"]["DATABASE_GROUP"],
        lat = "sparql",
        query = "binary_search_descendants"
    output:
        "data/query/mgdb/sparql/output_binary_search_descendants_for_{pid}.csv"
    log:
        "data/query/mgdb/sparql/binary_search_descendants_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/sparql/benchmark_binary_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_unary_search_descendants_with_clingo:
    input:
        facts = "data/raw/mgdb/clingo/facts.tsv"
    params:
        pid = "{pid}",
        lat = "clingo",
        query = "unary_search_descendants"
    output:
        "data/query/mgdb/clingo/output_unary_search_descendants_for_{pid}.csv"
    log:
        "data/query/mgdb/clingo/unary_search_descendants_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/clingo/benchmark_unary_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_binary_search_descendants_with_clingo:
    input:
        facts = "data/raw/mgdb/clingo/facts.tsv"
    params:
        pid = "{pid}",
        lat = "clingo",
        query = "binary_search_descendants"
    output:
        "data/query/mgdb/clingo/output_binary_search_descendants_for_{pid}.csv"
    log:
        "data/query/mgdb/clingo/binary_search_descendants_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/clingo/benchmark_binary_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_unary_search_descendants_with_bashlog:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv"
    params:
        pid = "{pid}",
        lat = "bashlog",
        query = "unary_search_descendants"
    output:
        "data/query/mgdb/bashlog/output_unary_search_descendants_for_{pid}.csv"
    log:
        "data/query/mgdb/bashlog/unary_search_descendants_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/bashlog/benchmark_unary_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_binary_search_descendants_with_bashlog:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv"
    params:
        pid = "{pid}",
        lat = "bashlog",
        query = "binary_search_descendants"
    output:
        "data/query/mgdb/bashlog/output_binary_search_descendants_for_{pid}.csv"
    log:
        "data/query/mgdb/bashlog/binary_search_descendants_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/bashlog/benchmark_binary_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_unary_search_descendants_with_logica:
    input:
        database = "data/external/mgdb/mgdb.db"
    params:
        pid = "{pid}",
        lat = "logica",
        query = "unary_search_descendants"
    output:
        "data/query/mgdb/logica/output_unary_search_descendants_for_{pid}.csv"
    log:
        "data/query/mgdb/logica/unary_search_descendants_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/logica/benchmark_unary_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_binary_search_descendants_with_logica:
    input:
        database = "data/external/mgdb/mgdb.db"
    params:
        pid = "{pid}",
        lat = "logica",
        query = "binary_search_descendants"
    output:
        "data/query/mgdb/logica/output_binary_search_descendants_for_{pid}.csv"
    log:
        "data/query/mgdb/logica/binary_search_descendants_for_{pid}.log"
    benchmark:
        repeat("data/query/mgdb/logica/benchmark_binary_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

# Lowest common ancestors
rule mgdb_lowest_common_ancestors_with_python:
    input:
        advised = "data/raw/mgdb/tsv/advised.tsv",
        dissertation = "data/raw/mgdb/tsv/dissertation.tsv",
        person = "data/raw/mgdb/tsv/person.tsv"
    params:
        pid1 = "{pid1}",
        pid2 = "{pid2}",
        lat = "python",
        query = "lowest_common_ancestors"
    output:
        "data/query/mgdb/python/output_lowest_common_ancestors_of_{pid1}_and_{pid2}.csv"
    log:
        "data/query/mgdb/python/lowest_common_ancestors_of_{pid1}_and_{pid2}.log"
    benchmark:
        repeat("data/query/mgdb/python/benchmark_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_lowest_common_ancestors_with_sql:
    input:
        database = "data/external/mgdb/mgdb.db"
    params:
        pid1 = "{pid1}",
        pid2 = "{pid2}",
        lat = "sql",
        query = "lowest_common_ancestors"
    output:
        "data/query/mgdb/sql/output_lowest_common_ancestors_of_{pid1}_and_{pid2}.csv"
    log:
        "data/query/mgdb/sql/lowest_common_ancestors_of_{pid1}_and_{pid2}.log"
    benchmark:
        repeat("data/query/mgdb/sql/benchmark_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_lowest_common_ancestors_with_cypher:
    input:
        optional_mgdb_cypher_input
    params:
        pid1 = "{pid1}",
        pid2 = "{pid2}",
        database_group = config["MGDB"]["DATA_SOURCE"]["CYPHER"]["DATABASE_GROUP"],
        lat = "cypher",
        query = "lowest_common_ancestors"
    output:
        "data/query/mgdb/cypher/output_lowest_common_ancestors_of_{pid1}_and_{pid2}.csv"
    log:
        "data/query/mgdb/cypher/lowest_common_ancestors_of_{pid1}_and_{pid2}.log"
    benchmark:
        repeat("data/query/mgdb/cypher/benchmark_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_lowest_common_ancestors_with_sparql:
    input:
        optional_mgdb_sparql_input
    params:
        pid1 = "{pid1}",
        pid2 = "{pid2}",
        namespace="MGDB",
        database_group = config["MGDB"]["DATA_SOURCE"]["SPARQL"]["DATABASE_GROUP"],
        lat = "sparql",
        query = "lowest_common_ancestors"
    output:
        "data/query/mgdb/sparql/output_lowest_common_ancestors_of_{pid1}_and_{pid2}.csv"
    log:
        "data/query/mgdb/sparql/lowest_common_ancestors_of_{pid1}_and_{pid2}.log"
    benchmark:
        repeat("data/query/mgdb/sparql/benchmark_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_lowest_common_ancestors_with_clingo:
    input:
        facts = "data/raw/mgdb/clingo/facts.tsv"
    params:
        pid1 = "{pid1}",
        pid2 = "{pid2}",
        lat = "clingo",
        query = "lowest_common_ancestors"
    output:
        "data/query/mgdb/clingo/output_lowest_common_ancestors_of_{pid1}_and_{pid2}.csv"
    log:
        "data/query/mgdb/clingo/lowest_common_ancestors_of_{pid1}_and_{pid2}.log"
    benchmark:
        repeat("data/query/mgdb/clingo/benchmark_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_lowest_common_ancestors_with_logica:
    input:
        database = "data/external/mgdb/mgdb.db"
    params:
        pid1 = "{pid1}",
        pid2 = "{pid2}",
        lat = "logica",
        query = "lowest_common_ancestors"
    output:
        "data/query/mgdb/logica/output_lowest_common_ancestors_of_{pid1}_and_{pid2}.csv"
    log:
        "data/query/mgdb/logica/lowest_common_ancestors_of_{pid1}_and_{pid2}.log"
    benchmark:
        repeat("data/query/mgdb/logica/benchmark_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_lowest_common_ancestors_with_bashlog:
    input:
        person = "data/raw/mgdb/bashlog/person.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        advised = "data/raw/mgdb/bashlog/advised.tsv"
    params:
        pid1 = "{pid1}",
        pid2 = "{pid2}",
        lat = "bashlog",
        query = "lowest_common_ancestors"
    output:
        "data/query/mgdb/bashlog/output_lowest_common_ancestors_of_{pid1}_and_{pid2}.csv"
    log:
        "data/query/mgdb/bashlog/lowest_common_ancestors_of_{pid1}_and_{pid2}.log"
    benchmark:
        repeat("data/query/mgdb/bashlog/benchmark_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule fish_lowest_common_ancestors_with_clingo:
    input:
        facts = "data/raw/fish/clingo/facts_{max_hamming_number}.tsv"
    params:
        pid1 = "{pid1}",
        pid2 = "{pid2}",
        lat = "clingo",
        query = "lowest_common_ancestors"
    output:
        "data/query/fish/{max_hamming_number}/clingo/output_lowest_common_ancestors_of_{pid1}_and_{pid2}.csv"
    log:
        "data/query/fish/{max_hamming_number}/clingo/lowest_common_ancestors_of_{pid1}_and_{pid2}.log"
    benchmark:
        repeat("data/query/fish/{max_hamming_number}/clingo/benchmark_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/fish_entry.py"

rule fish_lowest_common_ancestors_with_bashlog:
    input:
        fish = "data/raw/fish/bashlog/fish_{max_hamming_number}.tsv"
    params:
        pid1 = "{pid1}",
        pid2 = "{pid2}",
        lat = "bashlog",
        query = "lowest_common_ancestors"
    output:
        "data/query/fish/{max_hamming_number}/bashlog/output_lowest_common_ancestors_of_{pid1}_and_{pid2}.csv"
    log:
        "data/query/fish/{max_hamming_number}/bashlog/lowest_common_ancestors_of_{pid1}_and_{pid2}.log"
    benchmark:
        repeat("data/query/fish/{max_hamming_number}/bashlog/benchmark_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/fish_entry.py"

rule sail_lowest_common_ancestors_with_clingo:
    input:
        facts = "data/raw/sail/clingo/facts_{max_hamming_number}.tsv"
    params:
        pid1 = "{pid1}",
        pid2 = "{pid2}",
        lat = "clingo",
        query = "lowest_common_ancestors"
    output:
        "data/query/sail/{max_hamming_number}/clingo/output_lowest_common_ancestors_of_{pid1}_and_{pid2}.csv"
    log:
        "data/query/sail/{max_hamming_number}/clingo/lowest_common_ancestors_of_{pid1}_and_{pid2}.log"
    benchmark:
        repeat("data/query/sail/{max_hamming_number}/clingo/benchmark_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/sail_entry.py"

rule sail_lowest_common_ancestors_with_bashlog:
    input:
        sail = "data/raw/sail/bashlog/sail_{max_hamming_number}.tsv"
    params:
        pid1 = "{pid1}",
        pid2 = "{pid2}",
        lat = "bashlog",
        query = "lowest_common_ancestors"
    output:
        "data/query/sail/{max_hamming_number}/bashlog/output_lowest_common_ancestors_of_{pid1}_and_{pid2}.csv"
    log:
        "data/query/sail/{max_hamming_number}/bashlog/lowest_common_ancestors_of_{pid1}_and_{pid2}.log"
    benchmark:
        repeat("data/query/sail/{max_hamming_number}/bashlog/benchmark_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/sail_entry.py"

# Lowest common ancestors path
rule mgdb_lowest_common_ancestors_path_with_cypher:
    input:
        optional_mgdb_cypher_input
    params:
        pid1 = "{pid1}",
        pid2 = "{pid2}",
        database_group = config["MGDB"]["DATA_SOURCE"]["CYPHER"]["DATABASE_GROUP"],
        lat = "cypher",
        query = "lowest_common_ancestors_path"
    output:
        "data/query/mgdb/cypher/output_lowest_common_ancestors_path_of_{pid1}_and_{pid2}.csv"
    log:
        "data/query/mgdb/cypher/lowest_common_ancestors_path_of_{pid1}_and_{pid2}.log"
    benchmark:
        repeat("data/query/mgdb/cypher/benchmark_lowest_common_ancestors_path_of_{pid1}_and_{pid2}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_lowest_common_ancestors_path_with_clingo:
    input:
        facts = "data/raw/mgdb/clingo/facts.tsv"
    params:
        pid1 = "{pid1}",
        pid2 = "{pid2}",
        lat = "clingo",
        query = "lowest_common_ancestors_path"
    output:
        "data/query/mgdb/clingo/output_lowest_common_ancestors_path_of_{pid1}_and_{pid2}.csv"
    log:
        "data/query/mgdb/clingo/lowest_common_ancestors_path_of_{pid1}_and_{pid2}.log"
    benchmark:
        repeat("data/query/mgdb/clingo/benchmark_lowest_common_ancestors_path_of_{pid1}_and_{pid2}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

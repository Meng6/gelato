# Preprare bashlog
rule mgdb_prepare_unary_search_ancestors_datalog_for_bashlog:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv"
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/bashlog/unary_search_ancestors_for_{pid}.dlog"
    script:
        "../src/query/mgdb/prepare_unary_search_ancestors_datalog_for_bashlog.py"

rule mgdb_prepare_binary_search_ancestors_datalog_for_bashlog:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv"
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/bashlog/binary_search_ancestors_for_{pid}.dlog"
    threads: workflow.cores
    script:
        "../src/query/mgdb/prepare_binary_search_ancestors_datalog_for_bashlog.py"

rule mgdb_prepare_unary_search_descendants_datalog_for_bashlog:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv"
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/bashlog/unary_search_descendants_for_{pid}.dlog"
    script:
        "../src/query/mgdb/prepare_unary_search_descendants_datalog_for_bashlog.py"

rule mgdb_prepare_binary_search_descendants_datalog_for_bashlog:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv"
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/bashlog/binary_search_descendants_for_{pid}.dlog"
    script:
        "../src/query/mgdb/prepare_binary_search_descendants_datalog_for_bashlog.py"

rule mgdb_prepare_lowest_common_ancestors_datalog_for_bashlog:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv"
    params:
        pid1 = "{pid1}",
        pid2 = "{pid2}"
    output:
        "data/query/mgdb/bashlog/lowest_common_ancestors_of_{pid1}_and_{pid2}.dlog"
    script:
        "../src/query/mgdb/prepare_lowest_common_ancestors_datalog_for_bashlog.py"

rule mgdb_prepare_bashscript_for_bashlog_1:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv",
        datalog = "data/query/mgdb/bashlog/{query}_for_{pid}.dlog"
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/bashlog/{query}_for_{pid}.sh"
    threads: workflow.cores
    script:
        "../src/query/mgdb/prepare_bashscript_for_bashlog.sh"

rule mgdb_prepare_bashscript_for_bashlog_2:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv",
        datalog = "data/query/mgdb/bashlog/{query}_of_{pid1}_and_{pid2}.dlog"
    params:
        pid1 = "{pid1}",
        pid2 = "{pid2}"
    output:
        "data/query/mgdb/bashlog/{query}_of_{pid1}_and_{pid2}.sh"
    threads: workflow.cores
    script:
        "../src/query/mgdb/prepare_bashscript_for_bashlog.sh"

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
        "data/query/mgdb/python/output_unary_search_ancestors_for_{pid}.txt"
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
        "data/query/mgdb/python/output_binary_search_ancestors_for_{pid}.txt"
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
        "data/query/mgdb/sql/output_unary_search_ancestors_for_{pid}.txt"
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
        "data/query/mgdb/sql/output_binary_search_ancestors_for_{pid}.txt"
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
        "data/query/mgdb/cypher/output_unary_search_ancestors_for_{pid}.txt"
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
        "data/query/mgdb/cypher/output_binary_search_ancestors_for_{pid}.txt"
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
        "data/query/mgdb/sparql/output_unary_search_ancestors_for_{pid}.txt"
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
        "data/query/mgdb/sparql/output_binary_search_ancestors_for_{pid}.txt"
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
        "data/query/mgdb/clingo/output_unary_search_ancestors_for_{pid}.txt"
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
        "data/query/mgdb/clingo/output_binary_search_ancestors_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/clingo/benchmark_binary_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_unary_search_ancestors_with_bashlog:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv",
        bashlog = "data/query/mgdb/bashlog/unary_search_ancestors_for_{pid}.sh"
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/bashlog/output_unary_search_ancestors_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/bashlog/benchmark_unary_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    threads: workflow.cores
    shell:
        "bash data/query/mgdb/bashlog/unary_search_ancestors_for_{wildcards.pid}.sh > {output}"

rule mgdb_binary_search_ancestors_with_bashlog:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv",
        bashlog = "data/query/mgdb/bashlog/binary_search_ancestors_for_{pid}.sh"
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/bashlog/output_binary_search_ancestors_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/bashlog/benchmark_binary_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    threads: workflow.cores
    shell:
        "bash data/query/mgdb/bashlog/binary_search_ancestors_for_{wildcards.pid}.sh > {output}"

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
        "data/query/mgdb/python/output_unary_search_descendants_for_{pid}.txt"
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
        "data/query/mgdb/python/output_binary_search_descendants_for_{pid}.txt"
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
        "data/query/mgdb/sql/output_unary_search_descendants_for_{pid}.txt"
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
        "data/query/mgdb/sql/output_binary_search_descendants_for_{pid}.txt"
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
        "data/query/mgdb/cypher/output_unary_search_descendants_for_{pid}.txt"
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
        "data/query/mgdb/cypher/output_binary_search_descendants_for_{pid}.txt"
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
        "data/query/mgdb/sparql/output_unary_search_descendants_for_{pid}.txt"
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
        "data/query/mgdb/sparql/output_binary_search_descendants_for_{pid}.txt"
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
        "data/query/mgdb/clingo/output_unary_search_descendants_for_{pid}.txt"
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
        "data/query/mgdb/clingo/output_binary_search_descendants_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/clingo/benchmark_binary_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_unary_search_descendants_with_bashlog:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv",
        bashlog = "data/query/mgdb/bashlog/unary_search_descendants_for_{pid}.sh"
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/bashlog/output_unary_search_descendants_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/bashlog/benchmark_unary_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    threads: workflow.cores
    shell:
        "bash data/query/mgdb/bashlog/unary_search_descendants_for_{wildcards.pid}.sh > {output}"

rule mgdb_binary_search_descendants_with_bashlog:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv",
        bashlog = "data/query/mgdb/bashlog/binary_search_descendants_for_{pid}.sh",
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/bashlog/output_binary_search_descendants_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/bashlog/benchmark_binary_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    threads: workflow.cores
    shell:
        "bash data/query/mgdb/bashlog/binary_search_descendants_for_{wildcards.pid}.sh > {output}"

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
        "data/query/mgdb/python/output_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt"
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
        "data/query/mgdb/sql/output_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt"
    benchmark:
        repeat("data/query/mgdb/sql/benchmark_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
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
        "data/query/mgdb/clingo/output_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt"
    benchmark:
        repeat("data/query/mgdb/clingo/benchmark_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_entry.py"

rule mgdb_lowest_common_ancestors_with_bashlog:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv",
        bashlog = "data/query/mgdb/bashlog/lowest_common_ancestors_of_{pid1}_and_{pid2}.sh"
    params:
        pid1 = "{pid1}",
        pid2 = "{pid2}"
    output:
        "data/query/mgdb/bashlog/output_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt"
    benchmark:
        repeat("data/query/mgdb/bashlog/benchmark_lowest_common_ancestors_of_{pid1}_and_{pid2}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    threads: workflow.cores
    shell:
        "bash data/query/mgdb/bashlog/lowest_common_ancestors_of_{wildcards.pid1}_and_{wildcards.pid2}.sh > {output}"

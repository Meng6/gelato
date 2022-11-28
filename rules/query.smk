# Preprare bashlog
rule mgdb_prepare_search_ancestors_datalog_for_bashlog:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv"
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/bashlog/search_ancestors_for_{pid}.dlog"
    script:
        "../src/query/mgdb_prepare_search_ancestors_datalog_for_bashlog.py"

rule mgdb_prepare_search_descendants_datalog_for_bashlog:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv"
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/bashlog/search_descendants_for_{pid}.dlog"
    script:
        "../src/query/mgdb_prepare_search_descendants_datalog_for_bashlog.py"

rule mgdb_prepare_bashscript_for_bashlog:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv",
        datalog = "data/query/mgdb/bashlog/{query}_for_{pid}.dlog"
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/bashlog/{query}_for_{pid}.sh"
    script:
        "../src/query/mgdb_prepare_bashscript_for_bashlog.sh"

# Query
# Search ancestors
rule mgdb_search_ancestors_with_python:
    input:
        advised = "data/raw/mgdb/tsv/advised.tsv",
        dissertation = "data/raw/mgdb/tsv/dissertation.tsv",
        person = "data/raw/mgdb/tsv/person.tsv"
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/python/output_search_ancestors_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/python/benchmark_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_search_ancestors_with_python.py"

rule mgdb_search_ancestors_with_sql:
    input:
        database = "data/external/mgdb/mgdb.db"
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/sql/output_search_ancestors_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/sql/benchmark_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_search_ancestors_with_sql.py"

rule mgdb_search_ancestors_with_cypher:
    input:
        optional_mgdb_cypher_input
    params:
        pid = "{pid}",
        database_group = config["MGDB"]["DATA_SOURCE"]["CYPHER"]["DATABASE_GROUP"]
    output:
        "data/query/mgdb/cypher/output_search_ancestors_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/cypher/benchmark_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_search_ancestors_with_cypher.py"

rule mgdb_search_ancestors_with_sparql:
    input:
        optional_mgdb_sparql_input
    params:
        pid = "{pid}",
        namespace="MGDB",
        database_group = config["MGDB"]["DATA_SOURCE"]["SPARQL"]["DATABASE_GROUP"]
    output:
        "data/query/mgdb/sparql/output_search_ancestors_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/sparql/benchmark_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_search_ancestors_with_sparql.py"

rule mgdb_search_ancestors_with_clingo:
    input:
        facts = "data/raw/mgdb/clingo/facts.tsv"
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/clingo/output_search_ancestors_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/clingo/benchmark_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_search_ancestors_with_clingo.py"

rule mgdb_search_ancestors_with_bashlog:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv",
        bashlog = "data/query/mgdb/bashlog/search_ancestors_for_{pid}.sh"
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/bashlog/output_search_ancestors_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/bashlog/benchmark_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    shell:
        "bash data/query/mgdb/bashlog/search_ancestors_for_{wildcards.pid}.sh > {output}"

# Search descendants
rule mgdb_search_descendants_with_python:
    input:
        advised = "data/raw/mgdb/tsv/advised.tsv",
        dissertation = "data/raw/mgdb/tsv/dissertation.tsv",
        person = "data/raw/mgdb/tsv/person.tsv"
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/python/output_search_descendants_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/python/benchmark_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_search_descendants_with_python.py"

rule mgdb_search_descendants_with_sql:
    input:
        database = "data/external/mgdb/mgdb.db"
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/sql/output_search_descendants_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/sql/benchmark_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_search_descendants_with_sql.py"

rule mgdb_search_descendants_with_cypher:
    input:
        optional_mgdb_cypher_input
    params:
        pid = "{pid}",
        database_group = config["MGDB"]["DATA_SOURCE"]["CYPHER"]["DATABASE_GROUP"]
    output:
        "data/query/mgdb/cypher/output_search_descendants_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/cypher/benchmark_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_search_descendants_with_cypher.py"

rule mgdb_search_descendants_with_sparql:
    input:
        optional_mgdb_sparql_input
    params:
        pid = "{pid}",
        namespace="MGDB",
        database_group = config["MGDB"]["DATA_SOURCE"]["SPARQL"]["DATABASE_GROUP"]
    output:
        "data/query/mgdb/sparql/output_search_descendants_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/sparql/benchmark_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_search_descendants_with_sparql.py"

rule mgdb_search_descendants_with_clingo:
    input:
        facts = "data/raw/mgdb/clingo/facts.tsv"
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/clingo/output_search_descendants_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/clingo/benchmark_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_search_descendants_with_clingo.py"

rule mgdb_search_descendants_with_bashlog:
    input:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv",
        bashlog = "data/query/mgdb/bashlog/search_descendants_for_{pid}.sh",
        optional_file = optional_mgdb_search_descendants_with_bashlog_input
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/bashlog/output_search_descendants_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/bashlog/benchmark_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    shell:
        "bash data/query/mgdb/bashlog/search_descendants_for_{wildcards.pid}.sh > {output}"

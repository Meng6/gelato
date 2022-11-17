rule mgdb_search_descendants_with_python:
    input:
        advised = "data/raw/mgdb/python/advised.tsv",
        dissertation = "data/raw/mgdb/python/dissertation.tsv",
        person = "data/raw/mgdb/python/person.tsv"
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/python/output_search_descendants_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/python/benchmark_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_search_descendants_with_python.py"

rule mgdb_search_descendants_with_sql:
    params:
        pid = "{pid}",
        database_group = config["MGDB"]["DATA_SOURCE"]["SQL"]["DATABASE_GROUP"]
    output:
        "data/query/mgdb/sql/output_search_descendants_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/sql/benchmark_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_search_descendants_with_sql.py"

rule mgdb_search_descendants_with_cypher:
    params:
        pid = "{pid}",
        database_group = config["MGDB"]["DATA_SOURCE"]["CYPHER"]["DATABASE_GROUP"]
    output:
        "data/query/mgdb/cypher/output_search_descendants_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/cypher/benchmark_search_descendants_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_search_descendants_with_cypher.py"

rule mgdb_search_ancestors_with_python:
    input:
        advised = "data/raw/mgdb/python/advised.tsv",
        dissertation = "data/raw/mgdb/python/dissertation.tsv",
        person = "data/raw/mgdb/python/person.tsv"
    params:
        pid = "{pid}"
    output:
        "data/query/mgdb/python/output_search_ancestors_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/python/benchmark_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_search_ancestors_with_python.py"

rule mgdb_search_ancestors_with_sql:
    params:
        pid = "{pid}",
        database_group = config["MGDB"]["DATA_SOURCE"]["SQL"]["DATABASE_GROUP"]
    output:
        "data/query/mgdb/sql/output_search_ancestors_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/sql/benchmark_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_search_ancestors_with_sql.py"

rule mgdb_search_ancestors_with_cypher:
    params:
        pid = "{pid}",
        database_group = config["MGDB"]["DATA_SOURCE"]["CYPHER"]["DATABASE_GROUP"]
    output:
        "data/query/mgdb/cypher/output_search_ancestors_for_{pid}.txt"
    benchmark:
        repeat("data/query/mgdb/cypher/benchmark_search_ancestors_for_{pid}.txt", config["BENCHMARK"]["REPEAT_TIMES"])
    script:
        "../src/query/mgdb_search_ancestors_with_cypher.py"

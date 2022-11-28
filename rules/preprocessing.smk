rule mgdb_pull_tsv_files:
    input:
        get_mgdb_external_data
    output:
        "data/raw/mgdb/tsv/{graph_data}.tsv"
    script:
        "../src/preprocess/mgdb_pull_tsv_files.py"

rule mgdb_preprocess_data_for_sparql:
    input:
        advised = "data/raw/mgdb/tsv/advised.tsv",
        dissertation = "data/raw/mgdb/tsv/dissertation.tsv",
        person = "data/raw/mgdb/tsv/person.tsv"
    output:
        "data/raw/mgdb/sparql/mgdb-ttl.txt"
    script:
        "../src/preprocess/mgdb_preprocess_data_for_sparql.py"

rule mgdb_preprocess_data_for_clingo:
    input:
        advised = "data/raw/mgdb/tsv/advised.tsv",
        dissertation = "data/raw/mgdb/tsv/dissertation.tsv",
        person = "data/raw/mgdb/tsv/person.tsv"
    output:
        "data/raw/mgdb/clingo/facts.tsv"
    script:
        "../src/preprocess/mgdb_preprocess_data_for_clingo.py"

rule mgdb_preprocess_data_for_bashlog:
    input:
        advised = "data/raw/mgdb/tsv/advised.tsv",
        dissertation = "data/raw/mgdb/tsv/dissertation.tsv",
        person = "data/raw/mgdb/tsv/person.tsv"
    output:
        advised = "data/raw/mgdb/bashlog/advised.tsv",
        dissertation = "data/raw/mgdb/bashlog/dissertation.tsv",
        person = "data/raw/mgdb/bashlog/person.tsv"
    script:
        "../src/preprocess/mgdb_preprocess_data_for_bashlog.py"

rule load_mgdb_data_to_sqlite:
    input:
        advised = "data/raw/mgdb/tsv/advised.tsv",
        dissertation = "data/raw/mgdb/tsv/dissertation.tsv",
        person = "data/raw/mgdb/tsv/person.tsv"
    output:
        "data/raw/mgdb/sql/mgdb.db"
    script:
        "../src/preprocess/load_mgdb_data_to_sqlite.py"

rule load_mgdb_data_to_neo4j:
    params:
        database_group = config["MGDB"]["DATA_SOURCE"]["CYPHER"]["DATABASE_GROUP"]
    output:
        touch("data/raw/mgdb/cypher/load_mgdb_data_to_neo4j.done")
    script:
        "../src/preprocess/load_mgdb_data_to_neo4j.py"

rule load_mgdb_data_to_blazegraph:
    input:
        "data/raw/mgdb/sparql/mgdb-ttl.txt"
    params:
        namespace="MGDB"
    output:
        touch("data/raw/mgdb/sparql/load_mgdb_data_to_blazegraph.done")
    shell:
        "curl -D- -H 'Content-Type: application/x-turtle-RDR' --upload-file '{input}' -X POST 'http://localhost:9999/blazegraph/namespace/{params.namespace}/sparql'"
    
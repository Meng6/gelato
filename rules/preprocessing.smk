rule mgdb_preprocess_data_for_python:
    input:
        get_mgdb_external_data
    output:
        "data/raw/mgdb/python/{graph_data}.tsv"
    script:
        "../src/preprocess/mgdb_preprocess_data_for_python.py"

# rule preprocess_mgdb_data_for_sql:
#     input:
#         get_external_data
#     output:
#         "data/raw/{graph}/python/{graph_data}.tsv"
#     script:
#         "../src/preprocess/preprocess_mgdb_data_for_sql.py"

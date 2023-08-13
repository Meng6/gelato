import sys
sys.path.append('.')
from tools.template import format_output

def python_entry(*args, **kwargs):

    import pandas as pd
    from fish import python_main

    # Load data
    fish = pd.read_csv(snakemake.input["fish"], sep="\t")

    # Query
    query_function = getattr(python_main, snakemake.params["query"])
    data = query_function(snakemake.params, fish)
    
    return data, None

def sql_entry(*args, **kwargs):

    import sqlite3
    from fish import sql_main

    # Connect to SQLite DB
    conn = sqlite3.connect(snakemake.input["database"])

    # Query
    query_function = getattr(sql_main, snakemake.params["query"])
    data, columns = query_function(snakemake.params, conn)

    # Close the connection
    conn.close ()

    return data, columns

def cypher_entry(credentials, database_group):

    from neo4j import GraphDatabase
    from fish import cypher_main

    # The URI should be in the form protocol://<server location>:<port>. 
    # The supported protocols in URI could either be bolt or neo4j. 
    # - bolt should be used when creating a driver connecting to the Neo4j instance directly. 
    # - neo4j should be used when creating a driver with built-in routing.
    uri = "{protocal}://{host}:{port}".format(protocal=credentials[database_group]["protocol"],
                                              host=credentials[database_group]["host"],
                                              port=credentials[database_group]["port"])
    # Connect to Neo4j DB
    neo4j_driver = GraphDatabase.driver(uri, auth=(credentials[database_group]["user"], credentials[database_group]["password"]))

    # Query
    query_function = getattr(cypher_main, snakemake.params["query"])
    data = query_function(snakemake.params, neo4j_driver)

    # Close the connection
    neo4j_driver.close()

    return data, None

def sparql_entry(credentials, database_group):

    from SPARQLWrapper import SPARQLWrapper
    from fish import sparql_main

    # Connect to Blazegraph DB
    conn = SPARQLWrapper("http://{host}:{port}/blazegraph/namespace/{namespace}/sparql".format(host=credentials[database_group]["host"],
                                                                                port=credentials[database_group]["port"],
                                                                                namespace=snakemake.params["namespace"]))

    # Query
    query_function = getattr(sparql_main, snakemake.params["query"])
    data = query_function(snakemake.params, conn)

    return data, None

def clingo_entry(*args, **kwargs):

    import pandas as pd
    import clingo
    from fish import clingo_main

    # Add facts
    ctl = clingo.Control()
    facts = pd.read_csv(snakemake.input["facts"], sep="\t").to_string(header=False, index=False)
    ctl.add("base", [], facts)

    # Query
    query_function = getattr(clingo_main, snakemake.params["query"])
    data, columns = query_function(snakemake.params, ctl)

    ctl.cleanup()

    return data, columns

def bashlog_entry(*args, **kwargs):

    from fish import bashlog_main
    import subprocess, threading

    # Thread ID
    tid = str(threading.current_thread().ident)

    # File paths: datalog and bashlog
    file_name = snakemake.output[0].replace("output_", "")[:-4]
    datalog_script_path, bashlog_script_path = file_name + ".dlog", file_name + ".sh"

    # Query
    # Create a datalog query
    with open(datalog_script_path, mode="w", encoding="utf-8") as datalog_script:
        query_function = getattr(bashlog_main, snakemake.params["query"])
        (datalog_query, columns) = query_function(snakemake.params, snakemake.input["fish"])
        datalog_script.write(datalog_query)
    # Convert datalog to bashlog
    with open(bashlog_script_path, mode="w", encoding="utf-8") as bashlog_script:
        bashlog_query = subprocess.run(["curl", "--data-binary", "@" + datalog_script_path, "https://www.thomasrebele.org/projects/bashlog/api/datalog?query"], stdout=subprocess.PIPE, text=True).stdout
        bashlog_script.write(bashlog_query.replace(" rm -f tmp/*", "rm -rf tmp").replace("tmp", file_name + "_tmp" + tid))
    
    # Execute bashlog
    data = subprocess.run(["bash", bashlog_script_path], stdout=subprocess.PIPE, text=True).stdout

    return data, columns

# Load credentials
if "database_group" in snakemake.params.keys():
    import yaml
    database_group = snakemake.params["database_group"]
    with open("./credentials.yaml", "r", encoding="utf-8") as f:
        credentials = yaml.safe_load(f)
else:
    credentials, database_group = None, None
data, columns = locals()[snakemake.params["lat"] + "_entry"](credentials, database_group)
formated_data = format_output(data=data, columns=columns, lat=snakemake.params["lat"])
formated_data.to_csv(snakemake.output[0], index=False)

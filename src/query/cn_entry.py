import sys
sys.path.append('.')
from tools.template import format_output
import pandas as pd

def python_entry(credentials, database_group):
    
    from cn import python_main

    # Query
    query_function = getattr(python_main, snakemake.params["query"])
    if database_group:
        data = query_function(snakemake.params, credentials, database_group)
    else:
        data = query_function(snakemake.params, pd.read_csv(snakemake.input["data"]))
    
    return data, None

def sql_entry(*args, **kwargs):

    import sqlite3
    from cn import sql_main

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
    from cn import cypher_main

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
    from cn import sparql_main

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
    from cn import clingo_main

    # Add facts
    ctl = clingo.Control()
    facts = pd.read_csv(snakemake.input["facts"], sep="\t").to_string(header=False, index=False)
    ctl.add("base", [], facts)

    # Query
    query_function = getattr(clingo_main, snakemake.params["query"])
    data, columns = query_function(snakemake.params, ctl)

    ctl.cleanup()

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
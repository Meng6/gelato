def python_entry(*args, **kwargs):

    import pandas as pd
    from fish import python_main

    # Load data
    fish = pd.read_csv(snakemake.input["fish"], sep="\t")

    # Query
    query_function = getattr(python_main, snakemake.params["query"])
    query_function(snakemake.params, fish, snakemake.output[0])
    
    return

def sql_entry(*args, **kwargs):

    import sqlite3
    from fish import sql_main

    # Connect to SQLite DB
    conn = sqlite3.connect(snakemake.input["database"])

    # Query
    query_function = getattr(sql_main, snakemake.params["query"])
    query_function(snakemake.params, conn, snakemake.output[0])

    # Close the connection
    conn.close ()

    return

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
    query_function(snakemake.params["pid"], neo4j_driver, snakemake.output[0])

    # Close the connection
    neo4j_driver.close()

    return

def sparql_entry(credentials, database_group):

    from SPARQLWrapper import SPARQLWrapper
    from fish import sparql_main

    # Connect to Blazegraph DB
    conn = SPARQLWrapper("http://{host}:{port}/blazegraph/namespace/{namespace}/sparql".format(host=credentials[database_group]["host"],
                                                                                port=credentials[database_group]["port"],
                                                                                namespace=snakemake.params["namespace"]))

    # Query
    query_function = getattr(sparql_main, snakemake.params["query"])
    query_function(snakemake.params["pid"], conn, snakemake.output[0])

    return

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
    query_function(snakemake.params, ctl, snakemake.output[0])

    ctl.cleanup()

    return


# Load credentials
if "database_group" in snakemake.params.keys():
    import yaml
    database_group = snakemake.params["database_group"]
    with open("./credentials.yaml", "r", encoding="utf-8") as f:
        credentials = yaml.safe_load(f)
else:
    credentials, database_group = None, None
locals()[snakemake.params["lat"] + "_entry"](credentials, database_group)

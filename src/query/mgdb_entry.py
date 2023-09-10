import sys, os
sys.path.append('.')
from tools.template import format_output, timeit
import pandas as pd
from neo4j import GraphDatabase
from SPARQLWrapper import SPARQLWrapper
import sqlite3, clingo, subprocess

@timeit
def python_load_data(**kwargs):
    advised = pd.read_csv(snakemake.input["advised"], sep="\t")
    dissertation = pd.read_csv(snakemake.input["dissertation"], sep="\t")
    person = pd.read_csv(snakemake.input["person"], sep="\t")
    return advised, dissertation, person

@timeit
def sql_connect_sqlite(**kwargs):
    return sqlite3.connect(snakemake.input["database"])

@timeit
def cypher_connect_neo4j(credentials, database_group, **kwargs):
    # The URI should be in the form protocol://<server location>:<port>. 
    # The supported protocols in URI could either be bolt or neo4j. 
    # - bolt should be used when creating a driver connecting to the Neo4j instance directly. 
    # - neo4j should be used when creating a driver with built-in routing.
    uri = "{protocal}://{host}:{port}".format(protocal=credentials[database_group]["protocol"],
                                              host=credentials[database_group]["host"],
                                              port=credentials[database_group]["port"])
    # Connect to Neo4j DB
    neo4j_driver = GraphDatabase.driver(uri, auth=(credentials[database_group]["user"], credentials[database_group]["password"]))
    return neo4j_driver

@timeit
def sparql_connect_blazegraph(credentials, database_group, **kwargs):
    return SPARQLWrapper("http://{host}:{port}/blazegraph/namespace/{namespace}/sparql".format(host=credentials[database_group]["host"],
                                                                                port=credentials[database_group]["port"],
                                                                                namespace=snakemake.params["namespace"]))

@timeit
def clingo_load_facts(**kwargs):
    ctl = clingo.Control()
    facts = pd.read_csv(snakemake.input["facts"], sep="\t").to_string(header=False, index=False)
    ctl.add("base", [], facts)
    return ctl

@timeit
def datalog_datalog2bashlog(bashlog_script_path, datalog_script_path, file_name, tid, **kwargs):
    # Convert datalog to bashlog
    with open(bashlog_script_path, mode="w", encoding="utf-8") as bashlog_script:
        bashlog_query = subprocess.run(["curl", "--data-binary", "@" + datalog_script_path, "https://www.thomasrebele.org/projects/bashlog/api/datalog?query"], stdout=subprocess.PIPE, text=True).stdout
        bashlog_script.write(bashlog_query.replace(" rm -f tmp/*", "rm -rf tmp").replace("tmp", file_name + "_tmp" + tid))
    return

@timeit
def bashlog_execute_query(bashlog_script_path, **kwargs):
    return subprocess.run(["bash", bashlog_script_path], stdout=subprocess.PIPE, text=True).stdout

@timeit
def logica_execute_query(logica_script_path, output_name, **kwargs):
    logica_output = subprocess.run(["logica", logica_script_path, "run", output_name], capture_output=True, text=True).stdout
    data = [[e.strip() for e in row.split("|")[1:-1]] for row in logica_output.split("\n")[3:-2]]
    return data


def python_entry(*args, **kwargs):

    if os.path.isfile("mgdb/app/python_main.py"):
        from mgdb.app import python_main
    else:
        from mgdb import python_main

    log = []

    # Load data
    (advised, dissertation, person) = python_load_data(log=log)

    # Query
    query_function = getattr(python_main, snakemake.params["query"])
    data = query_function(snakemake.params, advised, dissertation, person, log=log)
    
    return data, None, log

def sql_entry(*args, **kwargs):

    if os.path.isfile("mgdb/app/sql_main.py"):
        from mgdb.app import sql_main
    else:
        from mgdb import sql_main

    log = []

    # Connect to SQLite DB
    conn = sql_connect_sqlite(log=log)

    # Query
    query_function = getattr(sql_main, snakemake.params["query"])
    data, columns = query_function(snakemake.params, conn, log=log)

    # Close the connection
    conn.close ()

    return data, columns, log

def cypher_entry(credentials, database_group):

    if os.path.isfile("mgdb/app/cypher_main.py"):
        from mgdb.app import cypher_main
    else:
        from mgdb import cypher_main

    log = []

    # Connect to Neo4j DB
    neo4j_driver = cypher_connect_neo4j(credentials, database_group, log=log)

    # Query
    query_function = getattr(cypher_main, snakemake.params["query"])
    data = query_function(snakemake.params, neo4j_driver, log=log)

    # Close the connection
    neo4j_driver.close()

    return data, None, log

def sparql_entry(credentials, database_group):

    if os.path.isfile("mgdb/app/sparql_main.py"):
        from mgdb.app import sparql_main
    else:
        from mgdb import sparql_main

    log = []

    # Connect to Blazegraph DB
    conn = sparql_connect_blazegraph(credentials, database_group, log=log)

    # Query
    query_function = getattr(sparql_main, snakemake.params["query"])
    data = query_function(snakemake.params, conn)

    return data, None, log

def clingo_entry(*args, **kwargs):

    if os.path.isfile("mgdb/app/clingo_main.py"):
        from mgdb.app import clingo_main
    else:
        from mgdb import clingo_main

    log = []

    # Add facts
    ctl = clingo_load_facts(log=log)

    # Query
    query_function = getattr(clingo_main, snakemake.params["query"])
    data, columns = query_function(snakemake.params, ctl)

    ctl.cleanup()

    return data, columns, log

def bashlog_entry(*args, **kwargs):

    if os.path.isfile("mgdb/app/bashlog_main.py"):
        from mgdb.app import bashlog_main
    else:
        from mgdb import bashlog_main
    import threading

    log = []

    # Thread ID
    tid = str(threading.current_thread().ident)

    # File paths: datalog and bashlog
    file_name = snakemake.output[0].replace("output_", "")[:-4]
    datalog_script_path, bashlog_script_path = file_name + ".dlog", file_name + ".sh"

    # Query
    # Create a datalog query
    with open(datalog_script_path, mode="w", encoding="utf-8") as datalog_script:
        query_function = getattr(bashlog_main, snakemake.params["query"])
        (datalog_query, columns) = query_function(snakemake.params, snakemake.input["advised"], snakemake.input["dissertation"], snakemake.input["person"])
        datalog_script.write(datalog_query)
    # Convert datalog to bashlog
    datalog_datalog2bashlog(bashlog_script_path, datalog_script_path, file_name, tid, log=log)

    # Execute bashlog
    data = bashlog_execute_query(bashlog_script_path, log=log)

    return data, columns, log

def logica_entry(*args, **kwargs):

    if os.path.isfile("mgdb/app/logica_main.py"):
        from mgdb.app import logica_main
    else:
        from mgdb import logica_main

    log = []

    # Query
    logica_script_path = snakemake.output[0].replace("output_", "").replace("txt", "l")
    with open(logica_script_path, mode="w", encoding="utf-8") as logica_script:
        query_function = getattr(logica_main, snakemake.params["query"])
        (logica_query, output_name, columns) = query_function(snakemake.params, snakemake.input["database"])
        logica_script.write(logica_query)

    # Execute logica
    data = logica_execute_query(logica_script_path, output_name, log=log)

    return data, columns, log

# Load credentials
if "database_group" in snakemake.params.keys():
    import yaml
    database_group = snakemake.params["database_group"]
    with open("./credentials.yaml", "r", encoding="utf-8") as f:
        credentials = yaml.safe_load(f)
else:
    credentials, database_group = None, None
(data, columns, log) = locals()[snakemake.params["lat"] + "_entry"](credentials, database_group)
formated_data = format_output(data=data, columns=columns, lat=snakemake.params["lat"])
# Save the formatted output
formated_data.to_csv(snakemake.output[0], index=False)
# Save execution time per phase
with open(snakemake.log[0], mode="a", encoding="utf-8") as log_file:
    log_file.write("\n".join(x for x in log) + "\n")

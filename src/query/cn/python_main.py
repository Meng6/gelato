import networkx as nx
import pandas as pd
# import numpy as np

def load_graph_from_neo4j_with_query(credentials, database_group, query):
    
    from neo4j import GraphDatabase
    
    # The URI should be in the form protocol://<server location>:<port>. 
    # The supported protocols in URI could either be bolt or neo4j. 
    # - bolt should be used when creating a driver connecting to the Neo4j instance directly. 
    # - neo4j should be used when creating a driver with built-in routing.
    uri = "{protocal}://{host}:{port}".format(protocal=credentials[database_group]["protocol"],
                                              host=credentials[database_group]["host"],
                                              port=credentials[database_group]["port"])
    # Connect to Neo4j DB
    neo4j_driver = GraphDatabase.driver(uri, auth=(credentials[database_group]["user"], credentials[database_group]["password"]))

    # Load data as Graph (NetworkX) from Neo4j DB
    res = neo4j_driver.session().run(query)

    G = nx.DiGraph()
    nodes = list(res.graph()._nodes.values())
    for node in nodes:
        G.add_node(node.element_id, labels=node._labels, properties=node._properties)
    rels = list(res.graph()._relationships.values())
    for rel in rels:
        G.add_edge(rel.start_node.element_id, rel.end_node.element_id, key=rel.element_id, type=rel.type, properties=rel._properties)
    
    # Close the connection
    neo4j_driver.close()

    return G

# Modules to be loaded to cn_entry.py script
def basic_measures(params, credentials, database_group):

    import sys
    sys.path.append('..')
    from tools.network_analysis import basic_measures

    graph = load_graph_from_neo4j_with_query(credentials, database_group, "MATCH (n)-[r]->(c) RETURN *;") # Paper Node only: "MATCH (n:Paper)-[r]->(c:Paper) RETURN *;"
    
    return basic_measures(graph, params["metrics"])

def influential_paper_detection_by_degreecentrality(params, data):

    graph = nx.DiGraph()
    graph = nx.from_pandas_edgelist(data, source="p1.nid", target="p2.nid", create_using=graph)
    node2score = nx.degree_centrality(graph)
    degree_centrality = pd.DataFrame.from_dict({"node": list(node2score.keys()), "score": list(node2score.values())})
    topk = degree_centrality.sort_values(by="score", ascending=False).head(params["k"])[["node"]]

    return topk

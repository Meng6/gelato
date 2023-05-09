import networkx as nx
import pandas as pd
import numpy as np

def compute_average_degree(graph, in_or_out="all"):
    """
    Computes the average degree of a directed or undirected graph.
    If graph is directed, the degree can be either "in", "out", or "all" (default).
    """
    if graph.is_directed() and in_or_out == "in":
        degrees = dict(graph.in_degree())
    elif graph.is_directed() and in_or_out == "out":
        degrees = dict(graph.out_degree())
    else:
        degrees = dict(graph.degree())
    return sum(degrees.values()) / len(degrees)

def compute_diameter(graph):
    # Check if the graph is strongly connected
    if nx.is_strongly_connected(graph):
        diameter = nx.diameter(graph)
    else:
        diameters = []
        for component in nx.strongly_connected_components(graph):
            subgraph = graph.subgraph(component)
            diameters.append(nx.diameter(subgraph))
        diameter = max(diameters)
    return diameter

def compute_average_path_length(graph):
    # Check if the graph is strongly connected
    if nx.is_strongly_connected(graph):
        # Compute the average shortest path length for the entire graph
        avg_path_length = nx.average_shortest_path_length(graph)
    else:
        # If the Graph is not strongly connected, compute average shortest path length for each strongly connected component.
        avg_path_lengths = []
        for component in nx.strongly_connected_components(graph):
            subgraph = graph.subgraph(component)
            avg_path_lengths.append(nx.average_shortest_path_length(subgraph))
        avg_path_length = np.mean(avg_path_lengths)
    return avg_path_length

def basic_measures(graph, metrics):

    basic_measures = pd.DataFrame(columns=["metric", "value"])
    if "numofnodes" in metrics:
        basic_measures.loc[len(basic_measures)] = ["numofnodes", graph.number_of_nodes()]
    if "numofedges" in metrics:
        basic_measures.loc[len(basic_measures)] = ["numofedges", graph.number_of_edges()]
    if "avgdegree" in metrics:
        basic_measures.loc[len(basic_measures)] = ["avgdegree", compute_average_degree(graph, in_or_out="all")]
    if "diameter" in metrics:
        basic_measures.loc[len(basic_measures)] = ["diameter", compute_diameter(graph)]
    if "clusteringcoefficient" in metrics:
        basic_measures.loc[len(basic_measures)] = ["clusteringcoefficient", nx.average_clustering(graph)]
    if "avgpathlength" in metrics:
        basic_measures.loc[len(basic_measures)] = ["avgpathlength", compute_average_path_length(graph)]

    return basic_measures

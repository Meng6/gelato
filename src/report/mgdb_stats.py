import pandas as pd

# Load data
advised = pd.read_csv(snakemake.input["advised"], sep="\t")
dissertation = pd.read_csv(snakemake.input["dissertation"], sep="\t")
person = pd.read_csv(snakemake.input["person"], sep="\t")

# Compute number of nodes and edges
num_of_nodes = len(person["pid"].unique()) + len(dissertation["did"].unique())
num_of_edges = advised.shape[0] + dissertation.shape[0]

# Save the stats of MGDB
mgdb_stats = pd.DataFrame({"graph": ["mgdb"], "num_of_nodes": [num_of_nodes], "num_of_edges": [num_of_edges]})
mgdb_stats.to_csv(snakemake.output[0], index=False)

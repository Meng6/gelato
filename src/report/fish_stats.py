import pandas as pd

rows = []
for graph_path in snakemake.input:
    
    graph_info = graph_path.split("/")[-1][:-4]
    facts = pd.read_csv(graph_path, sep="\t")

    node_set = set(facts["end"])
    node_set.update(set(facts["start"]))
    num_of_nodes = len(node_set)

    num_of_edges = facts.shape[0]

    rows.append((graph_info, num_of_nodes, num_of_edges))

# Save the stats of FISH
fish_stats = pd.DataFrame(columns=["graph", "num_of_nodes", "num_of_edges"], data=rows)
fish_stats.to_csv(snakemake.output[0], index=False)

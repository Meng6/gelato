import pandas as pd

# Load data
fish = pd.read_csv(snakemake.input["fish"], sep="\t").astype(str)
# Facts: fish(end,start,edge).
fish["fact"] = fish.apply(lambda row: 'fish(' + row["end"] + ',' + row["start"] + ',' + row["edge"] + ').', axis=1)

fish[["fact"]].to_csv(snakemake.output[0], sep="\t", index=False)

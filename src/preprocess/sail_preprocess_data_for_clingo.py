import pandas as pd

# Load data
sail = pd.read_csv(snakemake.input["sail"], sep="\t").astype(str)
# Facts: sail(end,start,edge).
sail["fact"] = sail.apply(lambda row: 'sail(' + row["end"] + ',' + row["start"] + ',' + row["edge"] + ').', axis=1)

sail[["fact"]].to_csv(snakemake.output[0], sep="\t", index=False)

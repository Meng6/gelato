import pandas as pd

# Load data
sail = pd.read_csv(snakemake.input["sail"], sep="\t").astype(str)

# Facts
sail["fact"] = sail.apply(lambda row: row["end"] + "\t<points_from>\t" + row["start"] + "\t<edge_color>\t" + row["edge"], axis=1)

fout = open(snakemake.output[0], mode="w", encoding="utf-8")
fout.write("\n".join(sail["fact"].values))
fout.close()

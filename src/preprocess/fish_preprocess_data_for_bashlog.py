import pandas as pd

# Load data
fish = pd.read_csv(snakemake.input["fish"], sep="\t").astype(str)

# Facts
fish["fact"] = fish.apply(lambda row: row["end"] + "\t<points_from>\t" + row["start"] + "\t<edge_color>\t" + row["edge"], axis=1)

fout = open(snakemake.output[0], mode="w", encoding="utf-8")
fout.write("\n".join(fish["fact"].values))
fout.close()

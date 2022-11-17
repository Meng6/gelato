import pandas as pd

FORMAT2SEP = {"csv": ",", "tsv": "\t"}

data = pd.read_csv(snakemake.input[0], sep=FORMAT2SEP[snakemake.input[0].split(".")[1]])
data.to_csv(snakemake.output[0], sep="\t", index=False)


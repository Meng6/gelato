import pandas as pd

# Preprocess string: 
# 1) Replace tab with space
# 2) Replace \ with a space
# 3) Replace consecutive spaces with one space
# 4) Replace " with '
def pstr(text):
    return text.replace('\t', ' ').replace('\\', ' ').replace('  ', ' ').replace('"', "'")

FORMAT2SEP = {"csv": ",", "tsv": "\t"}

# Convert files to TSV format
data = pd.read_csv(snakemake.input[0], sep=FORMAT2SEP[snakemake.input[0].split(".")[1]]).astype(str)
for col in data.columns:
    data[col] = data[col].apply(pstr)

data.to_csv(snakemake.output[0], sep="\t", index=False)

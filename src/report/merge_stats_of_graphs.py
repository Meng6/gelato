import pandas as pd

merged_stats = pd.DataFrame()
for file_path in snakemake.input:
    data = pd.read_csv(file_path)
    merged_stats = pd.concat([merged_stats, data], axis=0)

merged_stats.to_csv(snakemake.output[0], index=False)

import pandas as pd

benchmarks = pd.DataFrame()
for file_path in snakemake.input:
    lat = file_path.split("/")[3]
    query = file_path.split("/")[-1][10:-4]
    benchmark = pd.read_csv(file_path, sep="\t")[["s"]]
    benchmarks.loc[lat, query] = str(round(benchmark.mean().values[0], 2)) + "(" + str(round(benchmark.std().values[0], 2)) + ")"

benchmarks.index.name = "language_or_tool"
benchmarks.to_csv(snakemake.output[0], index=True)
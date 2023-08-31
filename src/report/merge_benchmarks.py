import pandas as pd

benchmarks = pd.DataFrame()
# Split input into two lists, one for total execution time and one for phase execution time
for file_path in snakemake.input:
    lat = file_path.split("/")[3]
    if file_path.endswith(".txt"):
        query = "total#" + file_path.split("/")[-1][10:-4]
        benchmark = pd.read_csv(file_path, sep="\t")[["s"]]
        benchmarks.loc[query, lat] = str(round(benchmark.mean().values[0], 2)) + "(" + str(round(benchmark.std().values[0], 2)) + ")"
    else: # endswith .log
        query = "phase#" + file_path.split("/")[-1][:-4]
        benchmark = pd.read_csv(file_path, header=None, names=["func", "time"])
        avg = benchmark.groupby("func")["time"].mean()
        std = benchmark.groupby("func")["time"].std()
        for phase in avg.index:
            benchmarks.loc[query + "#" + phase, lat] = str(round(avg[phase], 2)) + "(" + str(round(std[phase], 2)) + ")"

benchmarks.index.name = "query"
benchmarks.to_csv(snakemake.output[0], index=True)

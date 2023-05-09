import orjson

print("Start Loading...")
with open(snakemake.input[0]) as fin:
    data = orjson.loads(fin.read())
print("Finished Loading...")

filtered_papers = []
for i, paper in enumerate(data):
    if i % 10000 == 0:
        print(i)
    if "fos" in paper.keys():
        for fos in paper["fos"]:
            if fos["name"].lower() == "natural language processing":
                filtered_papers.append(paper)
                break
print("Finished Filtering...")

with open(snakemake.output[0], "wb") as fout:
    fout.write(orjson.dumps(filtered_papers))

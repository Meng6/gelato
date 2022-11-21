import pandas as pd

# Load data
advised = pd.read_csv(snakemake.input["advised"], sep="\t").astype(str)
dissertation = pd.read_csv(snakemake.input["dissertation"], sep="\t").astype(str)
person = pd.read_csv(snakemake.input["person"], sep="\t").astype(str)

# Person facts
person["fact"] = person.apply(lambda row: row["pid"] + "\t<name>\t" + row["name"] + "\n" + row["pid"] + "\t<country>\t" + row["country"], axis=1)

fout = open(snakemake.output["person"], mode="w", encoding="utf-8")
fout.write("\n".join(person["fact"].values))
fout.close()

# Dissertation facts
dissertation["fact"] = dissertation.apply(lambda row: row["did"] + "\t<author>\t" + row["author"] + "\n" + row["did"] + "\t<title>\t" + row["title"] + "\n" + row["did"] + "\t<university>\t" + row["university"] + "\n" + row["did"] + "\t<year>\t" + row["year"], axis=1)

fout = open(snakemake.output["dissertation"], mode="w", encoding="utf-8")
fout.write("\n".join(dissertation["fact"].values))
fout.close()

# Advised facts
advised["fact"] = advised.apply(lambda row: row["did"] + "\t<advised_by>\t" + row["advisor"] + "\t<ordered>\t" + row["advisororder"], axis=1)

fout = open(snakemake.output["advised"], mode="w", encoding="utf-8")
fout.write("\n".join(advised["fact"].values))
fout.close()
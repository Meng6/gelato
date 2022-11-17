import pandas as pd

def searchAncestorsPython(pid, ancestor_pids):

    dids = dissertation[dissertation["author"] == pid].did.tolist()
    tmp_pids = set(advised[advised["did"].isin(dids)].advisor.tolist()).difference(ancestor_pids)
    ancestor_pids.update(tmp_pids)
    for pid in tmp_pids:
        ancestor_pids = searchAncestorsPython(pid, ancestor_pids)

    return ancestor_pids

# Get params
pid = int(snakemake.params["pid"])
# Load data
advised = pd.read_csv(snakemake.input["advised"], sep="\t")
dissertation = pd.read_csv(snakemake.input["dissertation"], sep="\t")
person = pd.read_csv(snakemake.input["person"], sep="\t")

# Search ancestors
ancestor_pids = searchAncestorsPython(pid=pid, ancestor_pids=set())
ancestors = person[person["pid"].isin(ancestor_pids)]

# Save the results
fout = open(snakemake.output[0], mode="w", encoding="utf-8")
fout.write("The ancestors of " + str(pid) + " are: (" + str(len(ancestors)) + " ancestors)\n")
fout.write(ancestors["name"].to_string())
fout.close()
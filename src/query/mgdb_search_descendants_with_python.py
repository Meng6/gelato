import pandas as pd

def searchDescendantsPython(pid, student_pids):
    advised_dids = advised[advised["advisor"] == pid].did.tolist()
    tmp_pids = set(dissertation[dissertation["did"].isin(advised_dids)].author.tolist()).difference(student_pids)
    student_pids.update(tmp_pids)
    for pid in tmp_pids:
        student_pids = searchDescendantsPython(pid, student_pids)
    return student_pids

# Get params
pid = int(snakemake.params["pid"])
# Load data
advised = pd.read_csv(snakemake.input["advised"], sep="\t")
dissertation = pd.read_csv(snakemake.input["dissertation"], sep="\t")
person = pd.read_csv(snakemake.input["person"], sep="\t")

student_pids = searchDescendantsPython(pid=pid, student_pids=set())
students = person[person["pid"].isin(student_pids)]

# Save the results
fout = open(snakemake.output[0], mode="w", encoding="utf-8")
fout.write("The descendants of " + str(pid) + " are: (" + str(len(students)) + " descendants)\n")
fout.write(students["name"].to_string())
fout.close()
import pandas as pd

def extract_ancestor_pids(pid, ancestor_pids, advised, dissertation):
    dids = dissertation[dissertation["author"] == pid].did.tolist()
    tmp_pids = set(advised[advised["did"].isin(dids)].advisor.tolist()).difference(ancestor_pids)
    ancestor_pids.update(tmp_pids)
    for pid in tmp_pids:
        ancestor_pids = extract_ancestor_pids(pid, ancestor_pids, advised, dissertation)
    return ancestor_pids

def save_ancestors(ancestors, pid, output_file):
    fout = open(output_file, mode="w", encoding="utf-8")
    fout.write("The ancestors of " + str(pid) + " are: (" + str(len(ancestors)) + " ancestors)\n")
    fout.write(ancestors[["pid", "name"]].to_string(index=False, header=False))
    fout.close()
    return

def extract_ancestor_pid_pairs(pid, ancestor_pids, ancestor_pid_pairs, advised, dissertation):
    dids = dissertation[dissertation["author"] == pid].did.tolist()
    current_ancestor_pids = advised[advised["did"].isin(dids)].advisor.tolist()
    tmp_pids = set(current_ancestor_pids).difference(ancestor_pids)
    ancestor_pids.update(tmp_pids)
    tmp_pid_pairs = set([(pid, x) for x in current_ancestor_pids]) #(pid, advisor)
    ancestor_pid_pairs.update(tmp_pid_pairs)
    for pid in tmp_pids:
        ancestor_pid_pairs = extract_ancestor_pid_pairs(pid, ancestor_pids, ancestor_pid_pairs, advised, dissertation)
    return ancestor_pid_pairs

def save_ancestor_pairs(ancestor_pid_pairs, person, pid, output_file):
    # Extract (student, advisor) pairs with pid and name
    student_pids, advisor_pids = zip(*ancestor_pid_pairs)
    ancestor_pairs = pd.DataFrame({"student_pids": student_pids, "advisor_pids": advisor_pids})
    pid2name = person.set_index("pid")["name"].to_dict()
    ancestor_pairs["student_names"], ancestor_pairs["advisor_names"] = ancestor_pairs["student_pids"].map(pid2name), ancestor_pairs["advisor_pids"].map(pid2name)
    # Save pairs
    fout = open(output_file, mode="w", encoding="utf-8")
    fout.write("The ancestor pairs (student, advisor) of " + str(pid) + " are: (" + str(len(ancestor_pairs)) + " ancestor pairs)\n")
    fout.write(ancestor_pairs[["student_pids", "student_names", "advisor_pids", "advisor_names"]].to_string(index=False, header=False))
    fout.close()
    return

def extract_descendant_pids(pid, descendant_pids, advised, dissertation):
    advised_dids = advised[advised["advisor"] == pid].did.tolist()
    tmp_pids = set(dissertation[dissertation["did"].isin(advised_dids)].author.tolist()).difference(descendant_pids)
    descendant_pids.update(tmp_pids)
    for pid in tmp_pids:
        descendant_pids = extract_descendant_pids(pid, descendant_pids, advised, dissertation)
    return descendant_pids

def save_descendants(descendants, pid, output_file):
    fout = open(output_file, mode="w", encoding="utf-8")
    fout.write("The descendants of " + str(pid) + " are: (" + str(len(descendants)) + " descendants)\n")
    fout.write(descendants[["pid", "name"]].to_string(index=False, header=False))
    fout.close()
    return

def extract_descendant_pid_pairs(pid, descendant_pids, descendant_pid_pairs, advised, dissertation):
    advised_dids = advised[advised["advisor"] == pid].did.tolist()
    current_descendant_pids = dissertation[dissertation["did"].isin(advised_dids)].author.tolist()
    tmp_pids = set(current_descendant_pids).difference(descendant_pids)
    # descendant_pids.update(tmp_pids)
    tmp_pid_pairs = set([(x, pid) for x in current_descendant_pids]) #(student, pid)
    descendant_pid_pairs.update(tmp_pid_pairs)
    for pid in tmp_pids:
        descendant_pid_pairs = extract_descendant_pid_pairs(pid, descendant_pids, descendant_pid_pairs, advised, dissertation)
    return descendant_pid_pairs

def save_descendant_pairs(descendant_pid_pairs, person, pid, output_file):
    # Extract (student, advisor) pairs with pid and name
    student_pids, advisor_pids = zip(*descendant_pid_pairs)
    descendant_pairs = pd.DataFrame({"student_pids": student_pids, "advisor_pids": advisor_pids})
    pid2name = person.set_index("pid")["name"].to_dict()
    descendant_pairs["student_names"], descendant_pairs["advisor_names"] = descendant_pairs["student_pids"].map(pid2name), descendant_pairs["advisor_pids"].map(pid2name)
    # Save pairs
    fout = open(output_file, mode="w", encoding="utf-8")
    fout.write("The descendant pairs (student, advisor) of " + str(pid) + " are: (" + str(len(descendant_pairs)) + " descendant pairs)\n")
    fout.write(descendant_pairs[["student_pids", "student_names", "advisor_pids", "advisor_names"]].to_string(index=False, header=False))
    fout.close()
    return

# Modules to be loaded to mgdb_entry.py script
def unary_search_ancestors(pid, advised, dissertation, person, output_file):
    ancestor_pids = extract_ancestor_pids(pid, set(), advised, dissertation)
    ancestors = person[person["pid"].isin(ancestor_pids)]
    save_ancestors(ancestors, pid, output_file)
    return

def binary_search_ancestors(pid, advised, dissertation, person, output_file):
    ancestor_pid_pairs = extract_ancestor_pid_pairs(pid, set(), set(), advised, dissertation)
    save_ancestor_pairs(ancestor_pid_pairs, person, pid, output_file)
    return

def unary_search_descendants(pid, advised, dissertation, person, output_file):
    descendant_pids = extract_descendant_pids(pid, set(), advised, dissertation)
    descendants = person[person["pid"].isin(descendant_pids)]
    save_descendants(descendants, pid, output_file)
    return

def binary_search_descendants(pid, advised, dissertation, person, output_file):
    descendant_pid_pairs = extract_descendant_pid_pairs(pid, set(), set(), advised, dissertation)
    save_descendant_pairs(descendant_pid_pairs, person, pid, output_file)
    return
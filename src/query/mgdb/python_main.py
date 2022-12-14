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
    fout.write(ancestors["name"].to_string())
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
    fout.write(descendants["name"].to_string())
    fout.close()
    return

# Modules to be loaded to mgdb_entry.py script
def unary_search_ancestors(pid, advised, dissertation, person, output_file):
    ancestor_pids = extract_ancestor_pids(pid, set(), advised, dissertation)
    ancestors = person[person["pid"].isin(ancestor_pids)]
    save_ancestors(ancestors, pid, output_file)
    return

def unary_search_descendants(pid, advised, dissertation, person, output_file):
    descendant_pids = extract_descendant_pids(pid, set(), advised, dissertation)
    descendants = person[person["pid"].isin(descendant_pids)]
    save_descendants(descendants, pid, output_file)
    return

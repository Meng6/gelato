import sys
sys.path.append('..')
from tools.template import timeit
import pandas as pd

def extract_ancestor_pids(pid, ancestor_pids, advised, dissertation):
    dids = dissertation[dissertation["author"] == pid].did.tolist()
    tmp_pids = set(advised[advised["did"].isin(dids)].advisor.tolist()).difference(ancestor_pids)
    ancestor_pids.update(tmp_pids)
    for pid in tmp_pids:
        ancestor_pids = extract_ancestor_pids(pid, ancestor_pids, advised, dissertation)
    return ancestor_pids

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

def add_ancestor_names(ancestor_pid_pairs, person):
    # Extract (student, advisor) pairs with pid and name
    student_pids, advisor_pids = zip(*ancestor_pid_pairs)
    ancestor_pairs = pd.DataFrame({"student_pid": student_pids, "advisor_pid": advisor_pids})
    pid2name = person.set_index("pid")["name"].to_dict()
    ancestor_pairs["student_name"], ancestor_pairs["advisor_name"] = ancestor_pairs["student_pid"].map(pid2name), ancestor_pairs["advisor_pid"].map(pid2name)
    return ancestor_pairs[["student_pid", "student_name", "advisor_pid", "advisor_name"]]

def extract_descendant_pids(pid, descendant_pids, advised, dissertation):
    advised_dids = advised[advised["advisor"] == pid].did.tolist()
    tmp_pids = set(dissertation[dissertation["did"].isin(advised_dids)].author.tolist()).difference(descendant_pids)
    descendant_pids.update(tmp_pids)
    for pid in tmp_pids:
        descendant_pids = extract_descendant_pids(pid, descendant_pids, advised, dissertation)
    return descendant_pids

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

def add_descendant_names(descendant_pid_pairs, person):
    # Extract (student, advisor) pairs with pid and name
    student_pids, advisor_pids = zip(*descendant_pid_pairs)
    descendant_pairs = pd.DataFrame({"student_pid": student_pids, "advisor_pid": advisor_pids})
    pid2name = person.set_index("pid")["name"].to_dict()
    descendant_pairs["student_name"], descendant_pairs["advisor_name"] = descendant_pairs["student_pid"].map(pid2name), descendant_pairs["advisor_pid"].map(pid2name)
    return descendant_pairs[["student_pid", "student_name", "advisor_pid", "advisor_name"]]

# Check if pid's child (current descendant) is one of the common ancestors
def is_lowest_common_ancestor(pid, common_ancestors, advised, dissertation):
    advised_dids = advised[advised["advisor"] == pid].did.tolist()
    tmp_pids = dissertation[dissertation["did"].isin(advised_dids)].author.tolist()
    for pid in tmp_pids:
        if pid in common_ancestors:
            return False
    return True

# Modules to be loaded to mgdb_entry.py script
@timeit
def unary_search_ancestors(params, advised, dissertation, person, **kwargs):
    """
    Query unary_search_ancestors: extract ancestors of a researcher with the given PID.
    
    Args:
    - params (dict) : PID could be accessed via params["pid"]
    - advised (Pandas data frame) : did, advisororder, advisor
    - dissertation (Pandas data frame) : did, author, title, university, year
    - person (Pandas data frame) : pid, name, country, onlinedescendants

    Returns:
    - Pandas data frame: pid, name
    """
    ancestor_pids = extract_ancestor_pids(int(params["pid"]), set(), advised, dissertation)
    ancestors = person[person["pid"].isin(ancestor_pids)][["pid", "name"]]
    return ancestors

@timeit
def binary_search_ancestors(params, advised, dissertation, person, **kwargs):
    """
    Query binary_search_ancestors: extract ancestor pairs of a researcher with the given PID.
    
    Args:
    - params (dict) : PID could be accessed via params["pid"]
    - advised (Pandas data frame) : did, advisororder, advisor
    - dissertation (Pandas data frame) : did, author, title, university, year
    - person (Pandas data frame) : pid, name, country, onlinedescendants

    Returns:
    - Pandas data frame: student_pid, student_name, advisor_pid, advisor_name
    """
    ancestor_pid_pairs = extract_ancestor_pid_pairs(int(params["pid"]), set(), set(), advised, dissertation)
    ancestor_pairs = add_ancestor_names(ancestor_pid_pairs, person)
    return ancestor_pairs

@timeit
def unary_search_descendants(params, advised, dissertation, person, **kwargs):
    """
    Query unary_search_descendants: extract descendants of a researcher with the given PID.
    
    Args:
    - params (dict) : PID could be accessed via params["pid"]
    - advised (Pandas data frame) : did, advisororder, advisor
    - dissertation (Pandas data frame) : did, author, title, university, year
    - person (Pandas data frame) : pid, name, country, onlinedescendants

    Returns:
    - Pandas data frame: pid, name
    """
    descendant_pids = extract_descendant_pids(int(params["pid"]), set(), advised, dissertation)
    descendants = person[person["pid"].isin(descendant_pids)][["pid", "name"]]
    return descendants

@timeit
def binary_search_descendants(params, advised, dissertation, person, **kwargs):
    """
    Query binary_search_descendants: extract descendant pairs of a researcher with the given PID.
    
    Args:
    - params (dict) : PID could be accessed via params["pid"]
    - advised (Pandas data frame) : did, advisororder, advisor
    - dissertation (Pandas data frame) : did, author, title, university, year
    - person (Pandas data frame) : pid, name, country, onlinedescendants

    Returns:
    - Pandas data frame: student_pid, student_name, advisor_pid, advisor_name
    """
    descendant_pid_pairs = extract_descendant_pid_pairs(int(params["pid"]), set(), set(), advised, dissertation)
    descendant_pairs = add_descendant_names(descendant_pid_pairs, person)
    return descendant_pairs

@timeit
def lowest_common_ancestors(params, advised, dissertation, person, **kwargs):
    """
    Query lowest_common_ancestors: extract LCAs of a pair of researchers with the given PIDS, i.e., (PID1, PID2).
    
    Args:
    - params (dict) : the (PID1, PID2) pair could be accessed via params["pid1"] and params["pid2"]
    - advised (Pandas data frame) : did, advisororder, advisor
    - dissertation (Pandas data frame) : did, author, title, university, year
    - person (Pandas data frame) : pid, name, country, onlinedescendants

    Returns:
    - Pandas data frame: pid, name
    """
    pid1, pid2 = int(params["pid1"]), int(params["pid2"])
    ancestors_of_pid1 = extract_ancestor_pids(pid1, set(), advised, dissertation)
    ancestors_of_pid2 = extract_ancestor_pids(pid2, set(), advised, dissertation)
    common_ancestors = ancestors_of_pid1.intersection(ancestors_of_pid2)
    
    lowest_common_ancestor_pids = []
    for pid in common_ancestors:
        if is_lowest_common_ancestor(pid, common_ancestors, advised, dissertation):
            lowest_common_ancestor_pids.append(pid)
    lowest_common_ancestors = person[person["pid"].isin(lowest_common_ancestor_pids)][["pid", "name"]]
    return lowest_common_ancestors

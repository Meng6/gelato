def on_model_ancestors(m):
    return m.__str__().split("ancestors_of_pid_with_name")[1:]

def on_model_descendants(m):
    return m.__str__().split("descendants_of_pid_with_name")[1:]

def on_model_lowest_common_ancestors(m):
    return m.__str__().split("lowest_common_ancestors_with_name")[1:]

def on_model_lowest_common_ancestors_path(m):
    return m.__str__().split("path_lca")[1:]

# Modules to be loaded to mgdb_entry.py script
def unary_search_ancestors(params, ctl):
    
    pid = params["pid"]

    rules =  """
    ancestors_of_pid(X, {pid}) :- advise(X, {pid}, _).
    ancestors_of_pid(X, {pid}) :- advise(X, Y, _), ancestors_of_pid(Y, {pid}).
    ancestors_of_pid_with_name(PID, NAME) :- ancestors_of_pid(PID, {pid}), person(PID, NAME, _, _).
    
    #show ancestors_of_pid_with_name/2.
    """.format(pid=str(pid))

    ctl.add("rule", [], rules)
    ctl.ground([("base", []), ("rule", [])])
    ctl.configuration.solve.models="0"

    ancestors = []
    with ctl.solve(yield_=True) as handle:
        for model in handle:
            ancestors = ancestors + on_model_ancestors(model)
    
    columns = ["pid", "name"]

    return ancestors, columns

def binary_search_ancestors(params, ctl):

    pid = params["pid"]
    
    rules =  """
    ancestors_of_pid(X, {pid}, {pid}) :- advise(X, {pid}, _).
    ancestors_of_pid(X, {pid}, Y):- advise(X, Y, _), ancestors_of_pid(Y, {pid}, _).
    ancestors_of_pid_with_name(PID2, NAME2, PID1, NAME1) :- ancestors_of_pid(PID1, {pid}, PID2), person(PID1, NAME1, _, _), person(PID2, NAME2, _, _).
    
    #show ancestors_of_pid_with_name/4.
    """.format(pid=str(pid))

    ctl.add("rule", [], rules)
    ctl.ground([("base", []), ("rule", [])])
    ctl.configuration.solve.models="0"

    ancestors = []
    with ctl.solve(yield_=True) as handle:
        for model in handle:
            ancestors = ancestors + on_model_ancestors(model)
    
    columns = ["student_pid", "student_name", "advisor_pid", "advisor_name"]

    return ancestors, columns

def unary_search_descendants(params, ctl):

    pid = params["pid"]

    rules =  """
    descendants_of_pid({pid}, Y) :- advise({pid}, Y, _).
    descendants_of_pid({pid}, Z) :- descendants_of_pid({pid}, Y), advise(Y, Z, _).
    descendants_of_pid_with_name(PID, NAME) :- descendants_of_pid({pid}, PID), person(PID, NAME, _, _).
    
    #show descendants_of_pid_with_name/2.
    """.format(pid=str(pid))

    ctl.add("rule", [], rules)
    ctl.ground([("base", []), ("rule", [])])
    ctl.configuration.solve.models="0"

    descendants = []
    with ctl.solve(yield_=True) as handle:
        for model in handle:
            descendants = descendants + on_model_descendants(model)

    columns = ["pid", "name"]

    return descendants, columns

def binary_search_descendants(params, ctl):

    pid = params["pid"]

    rules =  """
    descendants_of_pid({pid}, Y, {pid}) :- advise({pid}, Y, _).
    descendants_of_pid({pid}, Z, Y) :- descendants_of_pid({pid}, Y, _), advise(Y, Z, _).
    descendants_of_pid_with_name(PID1, NAME1, PID2, NAME2) :- descendants_of_pid({pid}, PID1, PID2), person(PID1, NAME1, _, _), person(PID2, NAME2, _, _).
    
    #show descendants_of_pid_with_name/4.
    """.format(pid=str(pid))

    ctl.add("rule", [], rules)
    ctl.ground([("base", []), ("rule", [])])
    ctl.configuration.solve.models="0"

    descendants = []
    with ctl.solve(yield_=True) as handle:
        for model in handle:
            descendants = descendants + on_model_descendants(model)
    
    columns = ["student_pid", "student_name", "advisor_pid", "advisor_name"]

    return descendants, columns

def lowest_common_ancestors(params, ctl):

    pid1, pid2 = str(params["pid1"]), str(params["pid2"])

    rules = """
    ancestors_of_pid(X, {pid1}) :- advise(X, {pid1}, _).
    ancestors_of_pid(X, {pid1}) :- advise(X, Y, _), ancestors_of_pid(Y, {pid1}).

    ancestors_of_pid(X, {pid2}) :- advise(X, {pid2}, _).
    ancestors_of_pid(X, {pid2}) :- advise(X, Y, _), ancestors_of_pid(Y, {pid2}).

    common_ancestors({pid1}, {pid2}, X) :- ancestors_of_pid(X, {pid1}), ancestors_of_pid(X, {pid2}).
    not_lowest_common_ancestors({pid1}, {pid2}, X) :- common_ancestors({pid1}, {pid2}, X), common_ancestors({pid1}, {pid2}, Y), advise(X, Y, _).
    lowest_common_ancestors_with_name(X, NAME) :- common_ancestors({pid1}, {pid2}, X), not not_lowest_common_ancestors({pid1}, {pid2}, X), person(X, NAME, _, _).

    #show lowest_common_ancestors_with_name/2.
    """.format(pid1=pid1, pid2=pid2)

    ctl.add("rule", [], rules)
    ctl.ground([("base", []), ("rule", [])])
    ctl.configuration.solve.models="0"

    lowest_common_ancestors = []
    with ctl.solve(yield_=True) as handle:
        for model in handle:
            lowest_common_ancestors = lowest_common_ancestors + on_model_lowest_common_ancestors(model)

    columns = ["pid", "name"]

    return lowest_common_ancestors, columns

def lowest_common_ancestors_path(params, ctl):

    pid1, pid2 = str(params["pid1"]), str(params["pid2"])
    
    # # Allow multiple paths to each LCA
    # rules = """
    # ancestors_of_pid(X, {pid1}, X, {pid1}) :- advise(X, {pid1}, _).
    # ancestors_of_pid(X, {pid1}, X, Y) :- advise(X, Y, _), ancestors_of_pid(Y, {pid1}, M, N).
    # ancestors_of_pid(X, {pid1}, M, N) :- advise(X, Y, _), ancestors_of_pid(Y, {pid1}, M, N).

    # ancestors_of_pid(X, {pid2}, X, {pid2}) :- advise(X, {pid2}, _).
    # ancestors_of_pid(X, {pid2}, X, Y) :- advise(X, Y, _), ancestors_of_pid(Y, {pid2}, M, N).
    # ancestors_of_pid(X, {pid2}, M, N) :- advise(X, Y, _), ancestors_of_pid(Y, {pid2}, M, N).

    # common_ancestors({pid1}, {pid2}, X) :- ancestors_of_pid(X, {pid1}, _, _), ancestors_of_pid(X, {pid2}, _, _).
    # not_lowest_common_ancestors({pid1}, {pid2}, X) :- common_ancestors({pid1}, {pid2}, X), common_ancestors({pid1}, {pid2}, Y), advise(X, Y, _).
    # lowest_common_ancestors(X) :- common_ancestors({pid1}, {pid2}, X), not not_lowest_common_ancestors({pid1}, {pid2}, X).

    # path_lca(APID, ANAME, SPID, SNAME) :- lowest_common_ancestors(X), ancestors_of_pid(X, {pid1}, APID, SPID), person(APID, ANAME, _, _), person(SPID, SNAME, _, _).
    # path_lca(APID, ANAME, SPID, SNAME) :- path_lca(_, _, SPID, _), ancestors_of_pid(SPID, {pid1}, APID, SPID), person(APID, ANAME, _, _), person(SPID, SNAME, _, _).
    
    # path_lca(APID, ANAME, SPID, SNAME) :- lowest_common_ancestors(X), ancestors_of_pid(X, {pid2}, APID, SPID), person(APID, ANAME, _, _), person(SPID, SNAME, _, _).
    # path_lca(APID, ANAME, SPID, SNAME) :- path_lca(_, _, SPID, _), ancestors_of_pid(SPID, {pid2}, APID, SPID), person(APID, ANAME, _, _), person(SPID, SNAME, _, _).

    # #show path_lca/4.
    # """.format(pid1=pid1, pid2=pid2)

    # Only consider the shortest path to each LCA (can have multiple shortest paths as well)
    # ancestors_of_pid(ancestor, pid, _advisor, _student, layer)
    rules = """
    ancestors_of_pid(X, {pid1}, X, {pid1}, 1) :- advise(X, {pid1}, _).
    ancestors_of_pid(X, {pid1}, X, Y, L+1) :- advise(X, Y, _), ancestors_of_pid(Y, {pid1}, M, N, L).
    ancestors_of_pid(X, {pid1}, M, N, L+1) :- advise(X, Y, _), ancestors_of_pid(Y, {pid1}, M, N, L).

    ancestors_of_pid(X, {pid2}, X, {pid2}, 1) :- advise(X, {pid2}, _).
    ancestors_of_pid(X, {pid2}, X, Y, L+1) :- advise(X, Y, _), ancestors_of_pid(Y, {pid2}, M, N, L).
    ancestors_of_pid(X, {pid2}, M, N, L+1) :- advise(X, Y, _), ancestors_of_pid(Y, {pid2}, M, N, L).

    common_ancestors({pid1}, {pid2}, X) :- ancestors_of_pid(X, {pid1}, _, _, _), ancestors_of_pid(X, {pid2}, _, _, _).
    not_lowest_common_ancestors({pid1}, {pid2}, X) :- common_ancestors({pid1}, {pid2}, X), common_ancestors({pid1}, {pid2}, Y), advise(X, Y, _).
    lowest_common_ancestors(X) :- common_ancestors({pid1}, {pid2}, X), not not_lowest_common_ancestors({pid1}, {pid2}, X).

    not_shortest_path(X, {pid1}, L2) :- lowest_common_ancestors(X), ancestors_of_pid(X, {pid1}, APID, SPID, L1), ancestors_of_pid(X, {pid1}, APID, SPID, L2), L1 < L2.
    path_lca_with_layer(APID, ANAME, SPID, SNAME, L) :- lowest_common_ancestors(X), ancestors_of_pid(X, {pid1}, APID, SPID, L), not not_shortest_path(X, {pid1}, L), person(APID, ANAME, _, _), person(SPID, SNAME, _, _).
    path_lca_with_layer(APID, ANAME, SPID, SNAME, L-1) :- path_lca_with_layer(_, _, SPID, _, L), ancestors_of_pid(SPID, {pid1}, APID, SPID, L-1), person(APID, ANAME, _, _), person(SPID, SNAME, _, _).
    
    not_shortest_path(X, {pid2}, L2) :- lowest_common_ancestors(X), ancestors_of_pid(X, {pid2}, APID, SPID, L1), ancestors_of_pid(X, {pid2}, APID, SPID, L2), L1 < L2.
    path_lca_with_layer(APID, ANAME, SPID, SNAME, L) :- lowest_common_ancestors(X), ancestors_of_pid(X, {pid2}, APID, SPID, L), not not_shortest_path(X, {pid2}, L), person(APID, ANAME, _, _), person(SPID, SNAME, _, _).
    path_lca_with_layer(APID, ANAME, SPID, SNAME, L-1) :- path_lca_with_layer(_, _, SPID, _, L), ancestors_of_pid(SPID, {pid2}, APID, SPID, L-1), person(APID, ANAME, _, _), person(SPID, SNAME, _, _).

    path_lca(APID, ANAME, SPID, SNAME) :- path_lca_with_layer(APID, ANAME, SPID, SNAME, _).
    #show path_lca/4.
    """.format(pid1=pid1, pid2=pid2)

    ctl.add("rule", [], rules)
    ctl.ground([("base", []), ("rule", [])])
    ctl.configuration.solve.models="0"

    lowest_common_ancestors_path = []
    with ctl.solve(yield_=True) as handle:
        for model in handle:
            lowest_common_ancestors_path = lowest_common_ancestors_path + on_model_lowest_common_ancestors_path(model)

    columns = ["advisor_pid", "advisor_name", "student_pid", "student_name"]

    return lowest_common_ancestors_path, columns

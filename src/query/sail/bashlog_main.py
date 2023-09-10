# Modules to be loaded to sail_entry.py script
def lowest_common_ancestors(params, sail):
    datalog_query = """

    facts(S, P, O) :~ cat ./{sail_file_path}

    % advise(start_node, end_node).
    advise(START, END) :- facts(END, "<points_from>", START, _, _).

    % ancestor_of_pid(ancestors, PID1).
    ancestor_of_pid(X, {pid1}) :- advise(X, {pid1}).
    ancestor_of_pid(X, {pid1}) :- advise(X, Y), ancestor_of_pid(Y, {pid1}).

    % ancestor_of_pid(ancestors, PID2).
    ancestor_of_pid(X, {pid2}) :- advise(X, {pid2}).
    ancestor_of_pid(X, {pid2}) :- advise(X, Y), ancestor_of_pid(Y, {pid2}).

    % common_ancestors(PID1, PID2, common_ancestor).
    common_ancestors({pid1}, {pid2}, X) :- ancestor_of_pid(X, {pid1}), ancestor_of_pid(X, {pid2}).

    not_lowest_common_ancestors({pid1}, {pid2}, X) :- 
        common_ancestors({pid1}, {pid2}, X), 
        common_ancestors({pid1}, {pid2}, Y), 
        advise(X, Y).

    main(PID) :- 
        common_ancestors({pid1}, {pid2}, PID),
        not not_lowest_common_ancestors({pid1}, {pid2}, PID).

    """.format(
        sail_file_path=sail,
        pid1=params["pid1"],
        pid2=params["pid2"]
    )

    columns = ["number"]

    return datalog_query, columns
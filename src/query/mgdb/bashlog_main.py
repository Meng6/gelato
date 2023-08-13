# Modules to be loaded to mgdb_entry.py script
def unary_search_ancestors(params, advised, dissertation, person):
    datalog_query = """

    facts_p(S, P, O) :~ cat ./{person_file_path}
    facts_d(S, P, O) :~ cat ./{dissertation_file_path}
    facts_a(S, P1, O1, P2, O2) :~ cat ./{advised_file_path}

    % advise(advisor, student).
    advise(A, S) :-
        facts_a(DID, "<advised_by>", A, _, _),
        facts_d(DID, "<author>", S).

    % ancestor_of_pid(ancestors, PID).
    ancestor_of_pid(X, {pid}) :- advise(X, {pid}).
    ancestor_of_pid(X, {pid}) :- advise(X, Y), ancestor_of_pid(Y, {pid}).

    main(PID, NAME) :-
        ancestor_of_pid(PID, {pid}),
        facts_p(PID, "<name>", NAME).

    """.format(
        person_file_path=person,
        dissertation_file_path=dissertation,
        advised_file_path=advised,
        pid=params["pid"]
    )
    
    columns = ["pid", "name"]

    return datalog_query, columns

def binary_search_ancestors(params, advised, dissertation, person):
    datalog_query = """

    facts_p(S, P, O) :~ cat ./{person_file_path}
    facts_d(S, P, O) :~ cat ./{dissertation_file_path}
    facts_a(S, P1, O1, P2, O2) :~ cat ./{advised_file_path}

    % advise(advisor, student).
    advise(A, S) :-
        facts_a(DID, "<advised_by>", A, _, _),
        facts_d(DID, "<author>", S).

    % ancestor_of_pid(ancestor1, PID, ancestor2), where ancestor1 is ancestor2's advisor
    ancestor_of_pid(X, {pid}, {pid}) :- advise(X, {pid}).
    ancestor_of_pid(X, {pid}, Y) :- advise(X, Y), ancestor_of_pid(Y, {pid}, _).

    % Return: ancestor pairs (student, advisor).
    main(PID2, NAME2, PID1, NAME1) :-
        ancestor_of_pid(PID1, {pid}, PID2),
        facts_p(PID1, "<name>", NAME1),
        facts_p(PID2, "<name>", NAME2).

    """.format(
        person_file_path=person,
        dissertation_file_path=dissertation,
        advised_file_path=advised,
        pid=params["pid"]
    )
    
    columns = ["student_pid", "student_name", "advisor_pid", "advisor_name"]

    return datalog_query, columns

def unary_search_descendants(params, advised, dissertation, person):
    datalog_query = """

    facts_p(S, P, O) :~ cat ./{person_file_path}
    facts_d(S, P, O) :~ cat ./{dissertation_file_path}
    facts_a(S, P1, O1, P2, O2) :~ cat ./{advised_file_path}

    % advise(advisor, student).
    advise(A, S) :-
        facts_a(DID, "<advised_by>", A, _, _),
        facts_d(DID, "<author>", S).

    % descendants_for_pid(PID, descendants).
    descendants_for_pid({pid}, Y) :- advise({pid}, Y).
    descendants_for_pid({pid}, Z) :- descendants_for_pid({pid}, Y), advise(Y, Z).

    main(PID, NAME) :-
        descendants_for_pid({pid}, PID),
        facts_p(PID, "<name>", NAME).

    """.format(
        person_file_path=person,
        dissertation_file_path=dissertation,
        advised_file_path=advised,
        pid=params["pid"]
    )

    columns = ["pid", "name"]

    return datalog_query, columns

def binary_search_descendants(params, advised, dissertation, person):
    datalog_query = """

    facts_p(S, P, O) :~ cat ./{person_file_path}
    facts_d(S, P, O) :~ cat ./{dissertation_file_path}
    facts_a(S, P1, O1, P2, O2) :~ cat ./{advised_file_path}

    % advise(advisor, student).
    advise(A, S) :-
        facts_a(DID, "<advised_by>", A, _, _),
        facts_d(DID, "<author>", S).

    % descendants_for_pid(PID, descendant1, descendant2), where descendant2 is descendant1's advisor.
    descendants_for_pid({pid}, Y, {pid}) :- advise({pid}, Y).
    descendants_for_pid({pid}, Z, Y) :- descendants_for_pid({pid}, Y, _), advise(Y, Z).

    % Return: descendant pairs (student, advisor).
    main(PID1, NAME1, PID2, NAME2) :-
        descendants_for_pid({pid}, PID1, PID2),
        facts_p(PID1, "<name>", NAME1),
        facts_p(PID2, "<name>", NAME2).

    """.format(
        person_file_path=person,
        dissertation_file_path=dissertation,
        advised_file_path=advised,
        pid=params["pid"]
    )

    columns = ["student_pid", "student_name", "advisor_pid", "advisor_name"]

    return datalog_query, columns

def lowest_common_ancestors(params, advised, dissertation, person):
    datalog_query = """

    facts_p(S, P, O) :~ cat ./{person_file_path}
    facts_d(S, P, O) :~ cat ./{dissertation_file_path}
    facts_a(S, P1, O1, P2, O2) :~ cat ./{advised_file_path}

    % advise(advisor, student).
    advise(A, S) :-
        facts_a(DID, "<advised_by>", A, _, _),
        facts_d(DID, "<author>", S).

    % ancestor_of_pid(ancestors, PID1).
    ancestor_of_pid(X, {pid1}) :- advise(X, {pid1}).
    ancestor_of_pid(X, {pid1}) :- advise(X, Y), ancestor_of_pid(Y, {pid1}).

    % ancestor_of_pid(ancestors, PID2).
    ancestor_of_pid(X, {pid2}) :- advise(X, {pid2}).
    ancestor_of_pid(X, {pid2}) :- advise(X, Y), ancestor_of_pid(Y, {pid2}).

    % common_ancestors(PID1, PID2, common_ancestor).
    common_ancestors({pid1}, {pid2}, X) :- 
        ancestor_of_pid(X, {pid1}),
        ancestor_of_pid(X, {pid2}).

    nlowest_common_ancestors({pid1}, {pid2}, X) :- 
        common_ancestors({pid1}, {pid2}, X),
        common_ancestors({pid1}, {pid2}, Y),
        advise(X, Y).

    main(PID, NAME) :- 
        common_ancestors({pid1}, {pid2}, PID),
        not nlowest_common_ancestors({pid1}, {pid2}, PID),
        facts_p(PID, "<name>", NAME).

    """.format(
        person_file_path=person,
        dissertation_file_path=dissertation,
        advised_file_path=advised,
        pid1=params["pid1"],
        pid2=params["pid2"]
    )

    columns = ["pid", "name"]

    return datalog_query, columns

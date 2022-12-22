datalog_file = open(snakemake.output[0], encoding="utf-8", mode="w")

datalog_file.write("""
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
common_ancestors({pid1}, {pid2}, X) :- ancestor_of_pid(X, {pid1}), ancestor_of_pid(X, {pid2}).

not_lowest_common_ancestors({pid1}, {pid2}, X) :- 
   common_ancestors({pid1}, {pid2}, X), 
   common_ancestors({pid1}, {pid2}, Y), 
   advise(X, Y).

main(PID, NAME) :- 
   common_ancestors({pid1}, {pid2}, PID),
   not not_lowest_common_ancestors({pid1}, {pid2}, PID),
   facts_p(PID, "<name>", NAME).

""".format(person_file_path=snakemake.input["person"],
           dissertation_file_path=snakemake.input["dissertation"],
           advised_file_path=snakemake.input["advised"],
           pid1=snakemake.params["pid1"],
           pid2=snakemake.params["pid2"]))
datalog_file.close()

datalog_file = open(snakemake.output[0], encoding="utf-8", mode="w")

datalog_file.write("""
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
""".format(person_file_path=snakemake.input["person"],
           dissertation_file_path=snakemake.input["dissertation"],
           advised_file_path=snakemake.input["advised"],
           pid=snakemake.params["pid"]))
datalog_file.close()
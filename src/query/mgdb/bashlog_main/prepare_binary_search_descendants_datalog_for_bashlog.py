datalog_file = open(snakemake.output[0], encoding="utf-8", mode="w")

datalog_file.write("""
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
""".format(person_file_path=snakemake.input["person"],
           dissertation_file_path=snakemake.input["dissertation"],
           advised_file_path=snakemake.input["advised"],
           pid=snakemake.params["pid"]))
datalog_file.close()
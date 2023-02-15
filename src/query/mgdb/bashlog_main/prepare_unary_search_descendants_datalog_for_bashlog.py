datalog_file = open(snakemake.output[0], encoding="utf-8", mode="w")

datalog_file.write("""
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
""".format(person_file_path=snakemake.input["person"],
           dissertation_file_path=snakemake.input["dissertation"],
           advised_file_path=snakemake.input["advised"],
           pid=snakemake.params["pid"]))
datalog_file.close()
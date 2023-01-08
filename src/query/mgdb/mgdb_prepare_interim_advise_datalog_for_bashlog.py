datalog_file = open(snakemake.output[0], encoding="utf-8", mode="w")

datalog_file.write("""
facts_d(S, P, O) :~ cat ./{dissertation_file_path}
facts_a(S, P1, O1, P2, O2) :~ cat ./{advised_file_path}

% advise(advisor, student).
main(A, S) :- 
   facts_a(DID, "<advised_by>", A, _, _),
   facts_d(DID, "<author>", S).

""".format(dissertation_file_path=snakemake.input["dissertation"],
           advised_file_path=snakemake.input["advised"]))
datalog_file.close()

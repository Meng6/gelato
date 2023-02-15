datalog_file = open(snakemake.output[0], encoding="utf-8", mode="w")

datalog_file.write("""
facts_p(S, P, O) :~ cat ./{person_file_path}
% advise(advisor, student).
advise(A, S) :~ cat ./{interim_advise_file_path}
% common_ancestors(PID1, PID2, common_ancestor).
common_ancestors(PID1, PID2, X) :~ cat ./{interim_common_ancestors_file_path}

nlowest_common_ancestors({pid1}, {pid2}, X) :- 
   common_ancestors({pid1}, {pid2}, X), 
   common_ancestors({pid1}, {pid2}, Y),
   advise(X, Y).

main(PID, NAME) :- 
   common_ancestors({pid1}, {pid2}, PID),
   not nlowest_common_ancestors({pid1}, {pid2}, PID), 
   facts_p(PID, "<name>", NAME).

""".format(person_file_path=snakemake.input["person"],
           interim_advise_file_path=snakemake.input["interim_advise"],
           interim_common_ancestors_file_path=snakemake.input["interim_common_ancestors"],
           pid1=snakemake.params["pid1"],
           pid2=snakemake.params["pid2"]))
datalog_file.close()

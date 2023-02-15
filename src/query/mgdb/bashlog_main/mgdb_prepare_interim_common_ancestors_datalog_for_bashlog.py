datalog_file = open(snakemake.output[0], encoding="utf-8", mode="w")

datalog_file.write("""
% advise(advisor, student).
advise(A, S) :~ cat ./{interim_advise_file_path}

% ancestor_of_pid(ancestors, PID1).
ancestor_of_pid(X, {pid1}) :- advise(X, {pid1}).
ancestor_of_pid(X, {pid1}) :- advise(X, Y), ancestor_of_pid(Y, {pid1}).

% ancestor_of_pid(ancestors, PID2).
ancestor_of_pid(X, {pid2}) :- advise(X, {pid2}).
ancestor_of_pid(X, {pid2}) :- advise(X, Y), ancestor_of_pid(Y, {pid2}).

% common_ancestors(PID1, PID2, common_ancestor).
main({pid1}, {pid2}, X) :- 
   ancestor_of_pid(X, {pid1}), ancestor_of_pid(X, {pid2}).

""".format(interim_advise_file_path=snakemake.input["interim_advise"],
           pid1=snakemake.params["pid1"],
           pid2=snakemake.params["pid2"]))
datalog_file.close()

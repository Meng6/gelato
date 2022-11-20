import pandas as pd
import clingo

def on_model(m):
    lt = m.__str__().split("ancestors_of_pid_with_name")
    fout.write(" (" + str(len(lt)-1) + " ancestors)")
    fout.write("\n".join(lt))
    return

def searchAncestorsClingo(pid, ctl):

    rules =  """
    ancestors_of_pid(X, {pid}) :- advise(X, {pid}, _).
    ancestors_of_pid(X, {pid}) :- advise(X, Y, _), ancestors_of_pid(Y, {pid}).
    ancestors_of_pid_with_name(PID, NAME) :- ancestors_of_pid(PID, {pid}), person(PID, NAME, _, _).
    
    #show ancestors_of_pid_with_name/2.
    """.format(pid=str(pid))

    ctl.add("rule", [], rules)
    ctl.ground([("base", []), ("rule", [])])

    ctl.solve(on_model=on_model)

    return

# Add facts
pid = snakemake.params["pid"]
ctl = clingo.Control()
facts = pd.read_csv(snakemake.input["facts"], sep="\t").to_string(header=False, index=False)
ctl.add("base", [], facts)
# Save the results
fout = open(snakemake.output[0], mode="w", encoding="utf-8")
fout.write("The ancestors of " + str(pid) + " are:")
searchAncestorsClingo(pid, ctl)
ctl.cleanup()
fout.close()

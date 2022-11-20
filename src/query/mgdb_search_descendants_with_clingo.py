import pandas as pd
import clingo

def on_model(m):
    lt = m.__str__().split("descendants_of_pid_with_name")
    fout.write(" (" + str(len(lt)-1) + " descendants)")
    fout.write("\n".join(lt))
    return

def searchDescendantsClingo(pid, ctl):

    rules =  """
    descendants_of_pid({pid}, Y) :- advise({pid}, Y, _).
    descendants_of_pid({pid}, Z) :- descendants_of_pid({pid}, Y), advise(Y, Z, _).
    descendants_of_pid_with_name(PID, NAME) :- descendants_of_pid({pid}, PID), person(PID, NAME, _, _).
    
    #show descendants_of_pid_with_name/2.
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
fout.write("The descendants of " + str(pid) + " are:")
searchDescendantsClingo(pid, ctl)
ctl.cleanup()
fout.close()

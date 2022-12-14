def on_model_ancestors(m):
    lt = m.__str__().split("ancestors_of_pid_with_name")
    ancestors = " (" + str(len(lt)-1) + " ancestors)" + "\n".join(lt)
    return ancestors

def on_model_descendants(m):
    lt = m.__str__().split("descendants_of_pid_with_name")
    descendants = " (" + str(len(lt)-1) + " descendants)" + "\n".join(lt)
    return descendants

# Modules to be loaded to mgdb_entry.py script
def unary_search_ancestors(pid, ctl, output_file):
    
    rules =  """
    ancestors_of_pid(X, {pid}) :- advise(X, {pid}, _).
    ancestors_of_pid(X, {pid}) :- advise(X, Y, _), ancestors_of_pid(Y, {pid}).
    ancestors_of_pid_with_name(PID, NAME) :- ancestors_of_pid(PID, {pid}), person(PID, NAME, _, _).
    
    #show ancestors_of_pid_with_name/2.
    """.format(pid=str(pid))

    ctl.add("rule", [], rules)
    ctl.ground([("base", []), ("rule", [])])
    ctl.configuration.solve.models="0"

    ancestors = ""
    with ctl.solve(yield_=True) as handle:
        for model in handle:
            ancestors = ancestors + on_model_ancestors(model)
    
    # Save the results
    fout = open(output_file, mode="w", encoding="utf-8")
    fout.write("The ancestors of " + str(pid) + " are:")
    fout.write(ancestors)
    fout.close()

    return

def unary_search_descendants(pid, ctl, output_file):

    rules =  """
    descendants_of_pid({pid}, Y) :- advise({pid}, Y, _).
    descendants_of_pid({pid}, Z) :- descendants_of_pid({pid}, Y), advise(Y, Z, _).
    descendants_of_pid_with_name(PID, NAME) :- descendants_of_pid({pid}, PID), person(PID, NAME, _, _).
    
    #show descendants_of_pid_with_name/2.
    """.format(pid=str(pid))

    ctl.add("rule", [], rules)
    ctl.ground([("base", []), ("rule", [])])
    ctl.configuration.solve.models="0"

    descendants = ""
    with ctl.solve(yield_=True) as handle:
        for model in handle:
            descendants = descendants + on_model_descendants(model)

    # Save the results
    fout = open(output_file, mode="w", encoding="utf-8")
    fout.write("The descendants of " + str(pid) + " are:")
    fout.write(descendants)
    fout.close()

    return
def on_model_ancestors(m):
    lt = m.__str__().split("ancestors_of_pid_with_name")
    ancestors = " (" + str(len(lt)-1) + " ancestors)" + "\n".join(lt)
    return ancestors

def on_model_descendants(m):
    lt = m.__str__().split("descendants_of_pid_with_name")
    descendants = " (" + str(len(lt)-1) + " descendants)" + "\n".join(lt)
    return descendants

def on_model_lowest_common_ancestors(m):
    lt = m.__str__().split("lowest_common_ancestors_with_name")
    ancestors = "\n".join(lt)
    return ancestors

# Modules to be loaded to mgdb_entry.py script
def unary_search_ancestors(params, ctl, output_file):
    
    pid = params["pid"]

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

def binary_search_ancestors(params, ctl, output_file):

    pid = params["pid"]
    
    rules =  """
    ancestors_of_pid(X, {pid}, {pid}) :- advise(X, {pid}, _).
    ancestors_of_pid(X, {pid}, Y):- advise(X, Y, _), ancestors_of_pid(Y, {pid}, _).
    ancestors_of_pid_with_name(PID2, NAME2, PID1, NAME1) :- ancestors_of_pid(PID1, {pid}, PID2), person(PID1, NAME1, _, _), person(PID2, NAME2, _, _).
    
    #show ancestors_of_pid_with_name/4.
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
    fout.write("The ancestor pairs (student, advisor) of " + str(pid) + " are:")
    fout.write(ancestors.replace("ancestors", "ancestor pairs"))
    fout.close()

    return

def unary_search_descendants(params, ctl, output_file):

    pid = params["pid"]

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

def binary_search_descendants(params, ctl, output_file):

    pid = params["pid"]

    rules =  """
    descendants_of_pid({pid}, Y, {pid}) :- advise({pid}, Y, _).
    descendants_of_pid({pid}, Z, Y) :- descendants_of_pid({pid}, Y, _), advise(Y, Z, _).
    descendants_of_pid_with_name(PID1, NAME1, PID2, NAME2) :- descendants_of_pid({pid}, PID1, PID2), person(PID1, NAME1, _, _), person(PID2, NAME2, _, _).
    
    #show descendants_of_pid_with_name/4.
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
    fout.write("The descendant pairs (student, advisor) of " + str(pid) + " are:")
    fout.write(descendants.replace("descendants", "descendant pairs"))
    fout.close()

    return

def lowest_common_ancestors(params, ctl, output_file):

    pid1, pid2 = str(params["pid1"]), str(params["pid2"])

    rules = """
    ancestors_of_pid(X, {pid1}) :- advise(X, {pid1}, _).
    ancestors_of_pid(X, {pid1}) :- advise(X, Y, _), ancestors_of_pid(Y, {pid1}).

    ancestors_of_pid(X, {pid2}) :- advise(X, {pid2}, _).
    ancestors_of_pid(X, {pid2}) :- advise(X, Y, _), ancestors_of_pid(Y, {pid2}).

    common_ancestors({pid1}, {pid2}, X) :- ancestors_of_pid(X, {pid1}), ancestors_of_pid(X, {pid2}).
    not_lowest_common_ancestors({pid1}, {pid2}, X) :- common_ancestors({pid1}, {pid2}, X), common_ancestors({pid1}, {pid2}, Y), advise(X, Y, _).
    lowest_common_ancestors_with_name(X, NAME) :- common_ancestors({pid1}, {pid2}, X), not not_lowest_common_ancestors({pid1}, {pid2}, X), person(X, NAME, _, _).

    #show lowest_common_ancestors_with_name/2.
    """.format(pid1=pid1, pid2=pid2)

    ctl.add("rule", [], rules)
    ctl.ground([("base", []), ("rule", [])])
    ctl.configuration.solve.models="0"

    lowest_common_ancestors = ""
    with ctl.solve(yield_=True) as handle:
        for model in handle:
            lowest_common_ancestors = lowest_common_ancestors + on_model_lowest_common_ancestors(model)

    # Save the results
    fout = open(output_file, mode="w", encoding="utf-8")
    fout.write("The lowest common ancestor(s) of " + pid1 + " and " + pid2 + ":")
    fout.write(lowest_common_ancestors)
    fout.close()

    return

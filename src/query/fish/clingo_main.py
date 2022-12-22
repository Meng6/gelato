def on_model_lowest_common_ancestors(m):
    lt = m.__str__().split("lowest_common_ancestors")
    ancestors = "\n".join(lt)
    return ancestors

# Modules to be loaded to fish_entry.py script
def lowest_common_ancestors(params, ctl, output_file):

    pid1, pid2 = str(params["pid1"]), str(params["pid2"])

    rules = """

    advise(START, END) :- fish(END, START, _).

    ancestors_of_pid(X, {pid1}) :- advise(X, {pid1}).
    ancestors_of_pid(X, {pid1}) :- advise(X, Y), ancestors_of_pid(Y, {pid1}).

    ancestors_of_pid(X, {pid2}) :- advise(X, {pid2}).
    ancestors_of_pid(X, {pid2}) :- advise(X, Y), ancestors_of_pid(Y, {pid2}).

    common_ancestors({pid1}, {pid2}, X) :- ancestors_of_pid(X, {pid1}), ancestors_of_pid(X, {pid2}).
    not_lowest_common_ancestors({pid1}, {pid2}, X) :- common_ancestors({pid1}, {pid2}, X), common_ancestors({pid1}, {pid2}, Y), advise(X, Y).
    lowest_common_ancestors(X) :- common_ancestors({pid1}, {pid2}, X), not not_lowest_common_ancestors({pid1}, {pid2}, X).

    #show lowest_common_ancestors/1.
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

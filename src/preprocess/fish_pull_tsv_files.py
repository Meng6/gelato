import clingo

def on_model_fish(m):
    lt = m.__str__().replace("(", "").replace(")", "").replace(",", "\t").split("hamming")
    return "end\tstart\tedge" + "\n".join(lt)

ctl = clingo.Control()
facts = """
hamming(1, 1, 1). 
hamming(Y,X, 2) :- hamming(X,_, _), Y = 2*X, Y <= {max_hamming_number}. 
hamming(Y,X, 3) :- hamming(X,_, _), Y = 3*X, Y <= {max_hamming_number}. 
hamming(Y,X, 5) :- hamming(X,_, _), Y = 5*X, Y <= {max_hamming_number}.
""".format(max_hamming_number=snakemake.params["max_hamming_number"])
ctl.add("base", [], facts)
ctl.ground([("base", [])])
ctl.configuration.solve.models="0"

fish = ""
with ctl.solve(yield_=True) as handle:
    for model in handle:
        fish = fish + on_model_fish(model)

fout = open(snakemake.output[0], mode="w", encoding="utf-8")
fout.write(fish)
fout.close()

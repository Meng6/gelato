import clingo, re

def on_model_sail(m):
    lt = re.search(r'hamming(.+,.+,.+)', m.__str__()).group().replace(' ', '\n').replace("hamming(", "").replace(")", "").replace(",", "\t")
    return "end\tstart\tedge\n" + lt

ctl = clingo.Control()
facts = """
in_12(1). in_12(2). in_123(3). in_1235(5). 
in_123(X) :- in_12(X). 
in_1235(X) :- in_123(X).

hamming(1, 1, 1). 
hamming(Y, X, 2) :- hamming(X, _, F), in_12(F), Y = 2 * X , Y <= {max_hamming_number}. 
hamming(Y, X, 3) :- hamming(X, _, F), in_123(F), Y = 3 * X , Y <= {max_hamming_number}. 
hamming(Y, X, 5) :- hamming(X, _, F), in_1235(F), Y = 5 * X , Y <= {max_hamming_number}. 
""".format(max_hamming_number=snakemake.params["max_hamming_number"])
ctl.add("base", [], facts)
ctl.ground([("base", [])])
ctl.configuration.solve.models="0"

sail = ""
with ctl.solve(yield_=True) as handle:
    for model in handle:
        sail = sail + on_model_sail(model)

fout = open(snakemake.output[0], mode="w", encoding="utf-8")
fout.write(sail)
fout.close()

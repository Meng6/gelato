import pandas as pd

# Preprocess string: replace ' with \'
def pstr(text):
    return text.replace("'", "\\'")

# Load data
advised = pd.read_csv(snakemake.input["advised"], sep="\t").astype(str)
dissertation = pd.read_csv(snakemake.input["dissertation"], sep="\t").astype(str)
person = pd.read_csv(snakemake.input["person"], sep="\t").astype(str)

fout = open(snakemake.output[0], encoding="utf-8", mode="w")
fout.write("""
PREFIX : <http://MGDB.com/>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
""")

person["sparql"] = person.apply(lambda row: """
:p{pid} a :Person ;
    :name '{name}' ;
    :country '{country}' ;
    :onlinedescendants {onlinedescendants} .
""".format(pid=row["pid"], name=pstr(row["name"]), country=pstr(row["country"]), onlinedescendants=row["onlinedescendants"]), axis=1)
fout.write("".join(person["sparql"].values))

dissertation["sparql"] = dissertation.apply(lambda row: """
:d{did} a :Dissertation ;
    :title '{title}' ;
    :university '{university}' ;
    :year '{year}' .

:p{author} :writes :d{did} .
""".format(did=row["did"], author=row["author"], title=pstr(row["title"]), university=pstr(row["university"]), year=pstr(row["year"])), axis=1)
fout.write("".join(dissertation["sparql"].values))

advised["sparql"] = advised.apply(lambda row: """
<<:d{did} :advised_by :p{pid}>> :advisororder {order} .
""".format(did=row["did"], pid=row["advisor"], order=row["advisororder"]), axis=1)
fout.write("".join(advised["sparql"].values))
fout.close()
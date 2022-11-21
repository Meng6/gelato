import pandas as pd

# Preprocess string: 
# 1) Replace \ with a space
# 2) Replace consecutive spaces with one space
# 3) Replace " with '
def pstr(text):
    return text.replace('\\', ' ').replace('  ', ' ').replace('"', "'").replace('(', '[').replace(')', ']')

# Load data
advised = pd.read_csv(snakemake.input["advised"], sep="\t").astype(str)
dissertation = pd.read_csv(snakemake.input["dissertation"], sep="\t").astype(str)
person = pd.read_csv(snakemake.input["person"], sep="\t").astype(str)


# Facts: advise(advisor,student,advisororder).
did2author = dissertation.set_index("did")[["author"]].to_dict('index')
advised["student"] = advised["did"].apply(lambda x: did2author[x]["author"])
advised["fact"] = "advise(" + advised["advisor"] + "," + advised["student"] + "," + advised["advisororder"] + ")."

# Facts: dissertation(did,author,title,university,year).
dissertation["fact"] = dissertation.apply(lambda row: 'dissertation(' + row["did"] + ',' + pstr(row["author"]) + ',"' + pstr(row["title"]) + '","' + pstr(row["university"]) + '","' + pstr(row["year"]) + '").', axis=1)

# Facts: person(pid,name,country,onlinedescendants).
person["fact"] = person.apply(lambda row: 'person(' + row["pid"] + ',"' + pstr(row["name"]) + '","' + pstr(row["country"]) + '",' + row["onlinedescendants"] + ').', axis=1)

facts = pd.concat([advised[["fact"]], dissertation[["fact"]], person[["fact"]]], axis=0)
facts.to_csv(snakemake.output[0], sep="\t", index=False)

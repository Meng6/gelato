from SPARQLWrapper import SPARQLWrapper, JSON
import yaml

def searchDescendantsSPARQL(pid, conn):
    
    sparql_query = """
    PREFIX : <http://MGDB.com/>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

    SELECT ?s ?descendant
    WHERE {
        ?s a :Person.
        ?s (:writes|:advised_by)+/:advised_by :p"""+str(pid)+""" .
        ?s :name ?descendant.
    }
    """

    conn.setQuery(sparql_query)
    conn.setReturnFormat(JSON)
    try:
        students = conn.query().convert()["results"]["bindings"]
    except Exception as e:
        print(e)
    return students

# Load credentials
database_group = snakemake.params["database_group"]
with open("./credentials.yaml", "r", encoding="utf-8") as f:
    credentials = yaml.safe_load(f)
# Connect to Blazegraph DB
conn = SPARQLWrapper("http://{host}:{port}/blazegraph/namespace/{namespace}/sparql".format(host=credentials[database_group]["host"],
                                                                            port=credentials[database_group]["port"],
                                                                            namespace=snakemake.params["namespace"]))

# Search descendants
pid = snakemake.params["pid"]
students = searchDescendantsSPARQL(pid, conn)

# Save the results
fout = open(snakemake.output[0], mode="w", encoding="utf-8")
fout.write("The descendants of " + str(pid) + " are: (" + str(len(students)) + " descendants)\n")
for student in students:
    fout.write(student["s"]["value"].split("/p")[1] + ", " + student["descendant"]["value"] + "\n")
fout.close()
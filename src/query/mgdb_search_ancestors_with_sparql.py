from SPARQLWrapper import SPARQLWrapper, JSON
import yaml

def searchAncestorsSPARQL(pid, conn):
    
    sparql_query = """
    PREFIX : <http://MGDB.com/>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

    SELECT ?an ?ancestor
    WHERE {
        ?an a :Person.
        :p63244 :writes/(:writes|:advised_by)+ ?an .
        ?an :name ?ancestor
    }
    """

    conn.setQuery(sparql_query)
    conn.setReturnFormat(JSON)
    try:
        ancestors = conn.query().convert()["results"]["bindings"]
    except Exception as e:
        print(e)
    return ancestors

# Load credentials
database_group = snakemake.params["database_group"]
with open("./credentials.yaml", "r", encoding="utf-8") as f:
    credentials = yaml.safe_load(f)
# Connect to Blazegraph DB
conn = SPARQLWrapper("http://{host}:{port}/blazegraph/namespace/{namespace}/sparql".format(host=credentials[database_group]["host"],
                                                                            port=credentials[database_group]["port"],
                                                                            namespace=snakemake.params["namespace"]))

# Search ancestors
pid = snakemake.params["pid"]
ancestors = searchAncestorsSPARQL(pid, conn)

# Save the results
fout = open(snakemake.output[0], mode="w", encoding="utf-8")
fout.write("The ancestors of " + str(pid) + " are: (" + str(len(ancestors)) + " ancestors)\n")
for ancestor in ancestors:
    fout.write(ancestor["an"]["value"].split("/p")[1] + ", " + ancestor["ancestor"]["value"] + "\n")
fout.close()
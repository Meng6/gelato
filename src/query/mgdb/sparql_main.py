from SPARQLWrapper import JSON

def save_ancestors(ancestors, pid, output_file):
    fout = open(output_file, mode="w", encoding="utf-8")
    fout.write("The ancestors of " + str(pid) + " are: (" + str(len(ancestors)) + " ancestors)\n")
    for ancestor in ancestors:
        fout.write(ancestor["an"]["value"].split("/p")[1] + ", " + ancestor["ancestor"]["value"] + "\n")
    fout.close()
    return

def save_descendants(descendants, pid, output_file):
    fout = open(output_file, mode="w", encoding="utf-8")
    fout.write("The descendants of " + str(pid) + " are: (" + str(len(descendants)) + " descendants)\n")
    for descendant in descendants:
        fout.write(descendant["s"]["value"].split("/p")[1] + ", " + descendant["descendant"]["value"] + "\n")
    fout.close()
    return

# Modules to be loaded to mgdb_entry.py script
def search_ancestors(pid, conn, output_file):

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
    
    save_ancestors(ancestors, pid, output_file)

    return

def search_descendants(pid, conn, output_file):
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
        descendants = conn.query().convert()["results"]["bindings"]
    except Exception as e:
        print(e)
    
    save_descendants(descendants, pid, output_file)

    return
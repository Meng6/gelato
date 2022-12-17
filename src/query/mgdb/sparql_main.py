from SPARQLWrapper import JSON

def save_ancestors(ancestors, pid, output_file):
    fout = open(output_file, mode="w", encoding="utf-8")
    fout.write("The ancestors of " + str(pid) + " are: (" + str(len(ancestors)) + " ancestors)\n")
    for ancestor in ancestors:
        fout.write(ancestor["an"]["value"].split("/p")[1] + ", " + ancestor["ancestor"]["value"] + "\n")
    fout.close()
    return

def save_ancestor_pairs(ancestors, pid, output_file):
    fout = open(output_file, mode="w", encoding="utf-8")
    fout.write("The ancestor pairs (student, advisor) of " + str(pid) + " are: (" + str(len(ancestors)) + " ancestor pairs)\n")
    for ancestor in ancestors:
        fout.write(ancestor["an1"]["value"].split("/p")[1] + ", " + ancestor["an1name"]["value"] + ", " + ancestor["an2"]["value"].split("/p")[1] + ", " + ancestor["an2name"]["value"] + "\n")
    fout.close()
    return

def save_descendants(descendants, pid, output_file):
    fout = open(output_file, mode="w", encoding="utf-8")
    fout.write("The descendants of " + str(pid) + " are: (" + str(len(descendants)) + " descendants)\n")
    for descendant in descendants:
        fout.write(descendant["s"]["value"].split("/p")[1] + ", " + descendant["descendant"]["value"] + "\n")
    fout.close()
    return

def save_descendant_pairs(descendants, pid, output_file):
    fout = open(output_file, mode="w", encoding="utf-8")
    fout.write("The descendant pairs (student, advisor) of " + str(pid) + " are: (" + str(len(descendants)) + " descendant pairs)\n")
    for descendant in descendants:
        fout.write(descendant["de2"]["value"].split("/p")[1] + ", " + descendant["de2name"]["value"] + ", " + descendant["de1"]["value"].split("/p")[1] + ", " + descendant["de1name"]["value"] + "\n")
    fout.close()
    return

# Modules to be loaded to mgdb_entry.py script
def unary_search_ancestors(pid, conn, output_file):

    sparql_query = """
    PREFIX : <http://MGDB.com/>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

    SELECT DISTINCT ?an ?ancestor
    WHERE {
        ?an a :Person.
        :p""" + str(pid) + """ :writes/(:writes|:advised_by)+ ?an .
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

def binary_search_ancestors(pid, conn, output_file):

    sparql_query = """
    PREFIX : <http://MGDB.com/>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

    SELECT DISTINCT ?an1 ?an1name ?an2 ?an2name
    WHERE {
        ?d a :Dissertation.
        ?an1 a :Person.
        ?an2 a :Person.
        :p""" + str(pid) + """ :writes/(:writes|:advised_by)* ?d.
        ?an1 :writes ?d.
        ?d :advised_by ?an2.
        ?an1 :name ?an1name.
        ?an2 :name ?an2name.
    }
    """

    conn.setQuery(sparql_query)
    conn.setReturnFormat(JSON)
    try:
        ancestors = conn.query().convert()["results"]["bindings"]
    except Exception as e:
        print(e)
    
    save_ancestor_pairs(ancestors, pid, output_file)

    return

def unary_search_descendants(pid, conn, output_file):
    sparql_query = """
    PREFIX : <http://MGDB.com/>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

    SELECT DISTINCT ?s ?descendant
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

def binary_search_descendants(pid, conn, output_file):
    sparql_query = """
    PREFIX : <http://MGDB.com/>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

    SELECT DISTINCT ?de2 ?de2name ?de1 ?de1name
    WHERE {
        ?d a :Dissertation.
        ?de a :Person.
        ?de1 a :Person.
        ?de2 a :Person.
        ?de (:writes|:advised_by)+/:advised_by :p"""+str(pid)+""" .
        FILTER (?de1 IN (?de, :p"""+str(pid)+""" ) )
        ?de2 :writes ?d.
        ?d :advised_by ?de1.
        ?de2 :name ?de2name.
        ?de1 :name ?de1name.
    }
    """

    conn.setQuery(sparql_query)
    conn.setReturnFormat(JSON)
    try:
        descendants = conn.query().convert()["results"]["bindings"]
    except Exception as e:
        print(e)
    
    save_descendant_pairs(descendants, pid, output_file)

    return

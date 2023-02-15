from SPARQLWrapper import JSON

# Modules to be loaded to mgdb_entry.py script
def unary_search_ancestors(params, conn):

    pid = params["pid"]

    sparql_query = """
    PREFIX : <http://MGDB.com/>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

    SELECT DISTINCT ?pid ?name
    WHERE {
        ?pid a :Person.
        :p""" + str(pid) + """ :writes/(:writes|:advised_by)+ ?pid .
        ?pid :name ?name
    }
    """

    conn.setQuery(sparql_query)
    conn.setReturnFormat(JSON)
    try:
        ancestors = conn.query().convert()["results"]["bindings"]
    except Exception as e:
        print(e)

    return ancestors

def binary_search_ancestors(params, conn):

    pid = params["pid"]

    sparql_query = """
    PREFIX : <http://MGDB.com/>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

    SELECT DISTINCT ?student_pid ?student_name ?advisor_pid ?advisor_name
    WHERE {
        ?d a :Dissertation.
        ?student_pid a :Person.
        ?advisor_pid a :Person.
        :p""" + str(pid) + """ :writes/(:writes|:advised_by)* ?d.
        ?student_pid :writes ?d.
        ?d :advised_by ?advisor_pid.
        ?student_pid :name ?student_name.
        ?advisor_pid :name ?advisor_name.
    }
    """

    conn.setQuery(sparql_query)
    conn.setReturnFormat(JSON)
    try:
        ancestors = conn.query().convert()["results"]["bindings"]
    except Exception as e:
        print(e)

    return ancestors

def unary_search_descendants(params, conn):

    pid = params["pid"]

    sparql_query = """
    PREFIX : <http://MGDB.com/>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

    SELECT DISTINCT ?pid ?name
    WHERE {
        ?pid a :Person.
        ?pid (:writes|:advised_by)+/:advised_by :p"""+str(pid)+""" .
        ?pid :name ?name.
    }
    """

    conn.setQuery(sparql_query)
    conn.setReturnFormat(JSON)
    try:
        descendants = conn.query().convert()["results"]["bindings"]
    except Exception as e:
        print(e)
    
    return descendants

def binary_search_descendants(params, conn):

    pid = params["pid"]

    sparql_query = """
    PREFIX : <http://MGDB.com/>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

    SELECT DISTINCT ?student_pid ?student_name ?advisor_pid ?advisor_name
    WHERE {
        ?d a :Dissertation.
        ?de a :Person.
        ?advisor_pid a :Person.
        ?student_pid a :Person.
        ?de (:writes|:advised_by)+/:advised_by :p"""+str(pid)+""" .
        FILTER (?advisor_pid IN (?de, :p"""+str(pid)+""" ) )
        ?student_pid :writes ?d.
        ?d :advised_by ?advisor_pid.
        ?student_pid :name ?student_name.
        ?advisor_pid :name ?advisor_name.
    }
    """

    conn.setQuery(sparql_query)
    conn.setReturnFormat(JSON)
    try:
        descendants = conn.query().convert()["results"]["bindings"]
    except Exception as e:
        print(e)

    return descendants

def lowest_common_ancestors(params, conn):

    pid1, pid2 = params["pid1"], params["pid2"]

    # Need to be optimized: too slow
    sparql_query = """
    PREFIX : <http://MGDB.com/>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

    SELECT DISTINCT (?lca AS ?pid) ?name
    WHERE {
    ?lca a :Person .
    ?nlca a :Person .
    ?lca ^(:writes/(:writes|:advised_by)+)  :p"""+str(pid1)+""", :p"""+str(pid2)+""" .
    FILTER NOT EXISTS {
        ?nlca ^(:writes/(:writes|:advised_by)+)  :p"""+str(pid1)+""", :p"""+str(pid2)+""" ;
            :writes/(:writes|:advised_by)+ ?lca .
    }
    ?lca :name ?name .
    }
    """

    conn.setQuery(sparql_query)
    conn.setReturnFormat(JSON)
    try:
        lowest_common_ancestors = conn.query().convert()["results"]["bindings"]
    except Exception as e:
        print(e)

    return lowest_common_ancestors

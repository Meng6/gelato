import sys
sys.path.append('..')
from tools.template import timeit
from SPARQLWrapper import JSON

# Modules to be loaded to mgdb_entry.py script
@timeit
def unary_search_ancestors(params, conn, **kwargs):

    pid = params["pid"]

    sparql_query = """
    PREFIX : <http://MGDB.com/>

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

@timeit
def binary_search_ancestors(params, conn, **kwargs):

    pid = params["pid"]

    sparql_query = """
    PREFIX : <http://MGDB.com/>

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

@timeit
def unary_search_descendants(params, conn, **kwargs):

    pid = params["pid"]

    sparql_query = """
    PREFIX : <http://MGDB.com/>

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

@timeit
def binary_search_descendants(params, conn, **kwargs):

    pid = params["pid"]

    sparql_query = """
    PREFIX : <http://MGDB.com/>

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

@timeit
def lowest_common_ancestors(params, conn, **kwargs):

    pid1, pid2 = params["pid1"], params["pid2"]

    # Below SPARQL query is too slow
    sparql_query = """
    PREFIX : <http://MGDB.com/>

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

    # Below SPARQL query is faster
    sparql_query = """
    PREFIX : <http://MGDB.com/>

    SELECT DISTINCT (?lca1 AS ?pid) ?name
    WHERE {
    ?lca1 a :Person .
    ?lca2 a :Person .
    ?lca1 ^(:writes/(:writes|:advised_by)+) :p"""+str(pid1)+""" .
    ?lca2 ^(:writes/(:writes|:advised_by)+) :p"""+str(pid2)+""" .
    FILTER (?lca1 = ?lca2)
    FILTER NOT EXISTS {
        ?nlca1 a :Person .
        ?nlca2 a :Person .
        ?nlca1 ^(:writes/(:writes|:advised_by)+) :p"""+str(pid1)+""" .
        ?nlca2 ^(:writes/(:writes|:advised_by)+) :p"""+str(pid2)+""" .
        FILTER (?nlca1 = ?nlca2)
        ?nlca1 :writes/(:writes|:advised_by)+ ?lca1 .
    }
    ?lca1 :name ?name .
    }
    """

    conn.setQuery(sparql_query)
    conn.setReturnFormat(JSON)
    try:
        lowest_common_ancestors = conn.query().convert()["results"]["bindings"]
    except Exception as e:
        print(e)

    return lowest_common_ancestors

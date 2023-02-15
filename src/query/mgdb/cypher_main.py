# Modules to be loaded to mgdb_entry.py script
def unary_search_ancestors(params, neo4j_driver):

    pid = params["pid"]

    # Slow: around 7 sec
    # Ref: https://hub.packtpub.com/advanced-cypher-tricks
    # cypher_query = """
    #                 MATCH (p2:Person)-[*]->(d:Dissertation)-[:ADVISED_BY]->(p1:Person)
    #                 WHERE p2.pid='""" + str(pid) + """'
    #                 RETURN DISTINCT p1.pid, p1.name;
    #             """

    cypher_query = """
                    MATCH (p2:Person {pid: '""" + str(pid) + """'})-[*]->(d:Dissertation)
                    WITH DISTINCT d
                    MATCH (d)-[:ADVISED_BY]-(p1:Person)
                    RETURN DISTINCT p1.pid AS pid, p1.name AS name;
                """

    with neo4j_driver.session() as sess:
        ancestors = sess.run(cypher_query).to_df()

    return ancestors

def binary_search_ancestors(params, neo4j_driver):

    # Ref: https://stackoverflow.com/questions/16611723/using-multiple-match-clauses-doesnt-return-any-result-in-neo4j-cypher-query

    pid = params["pid"]
    
    cypher_query = """
                    MATCH (p3:Person {pid: '""" + str(pid) + """'})-[*]->(d:Dissertation)
                    WITH DISTINCT d
                    MATCH (p2:Person)-[:WRITES]->(d)-[:ADVISED_BY]->(p1:Person)
                    RETURN DISTINCT p2.pid AS student_pid, p2.name AS student_name, p1.pid AS advisor_pid, p1.name AS advisor_name
                """
    
    with neo4j_driver.session() as sess:
        ancestors = sess.run(cypher_query).to_df()

    return ancestors

def unary_search_descendants(params, neo4j_driver):

    pid = params["pid"]

    cypher_query = """
                    MATCH (d:Dissertation)-[*]->(p1:Person {pid: '""" + str(pid) + """'})
                    WITH DISTINCT d
                    MATCH (p2:Person)-[:WRITES]-(d)
                    RETURN DISTINCT p2.pid AS pid, p2.name AS name;
                """

    with neo4j_driver.session() as sess:
        descendants = sess.run(cypher_query).to_df()

    return descendants

def binary_search_descendants(params, neo4j_driver):

    pid = params["pid"]

    cypher_query = """
                    MATCH (p1:Person)-[*]->(p3:Person {pid: '""" + str(pid) + """'})
                    WITH DISTINCT p1.pid AS de_pid
                    MATCH (p2:Person)-[:WRITES]->(d:Dissertation)-[:ADVISED_BY]->(p:Person)
                    WHERE p.pid=de_pid OR p.pid='""" + str(pid) + """'
                    RETURN DISTINCT p2.pid AS student_pid, p2.name AS student_name, p.pid AS advisor_pid, p.name AS advisor_name
                """
    
    with neo4j_driver.session() as sess:
        descendants = sess.run(cypher_query).to_df()

    return descendants

def lowest_common_ancestors(params, neo4j_driver):

    pid1, pid2 = str(params["pid1"]), str(params["pid2"])

    cypher_query = """
                    MATCH (p1:Person {pid: '""" + str(pid1) + """'})-[*]->(ca:Person)
                    WITH DISTINCT ca
                    MATCH (p2:Person {pid: '""" + str(pid2) + """'})-[*]->(ca)
                    WITH DISTINCT COLLECT(ca.pid) AS ca_pids
                    MATCH (p3:Person)-[*]->(p4:Person)
                    WHERE p3.pid IN ca_pids AND p4.pid IN ca_pids
                    WITH COLLECT(DISTINCT p3.pid) AS student_pids, COLLECT(DISTINCT p4.pid) AS advisor_pids
                    MATCH (lca:Person)
                    WHERE lca.pid IN [x IN student_pids WHERE NOT(x IN advisor_pids)]
                    RETURN lca.pid AS pid, lca.name AS name
                """
    
    with neo4j_driver.session() as sess:
        lowest_common_ancestors = sess.run(cypher_query).to_df()

    return lowest_common_ancestors

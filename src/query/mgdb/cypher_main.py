def save_ancestors(ancestors, pid, output_file, pairs_flag):
    fout = open(output_file, mode="w", encoding="utf-8")
    if pairs_flag:
        fout.write("The ancestor pairs (student, advisor) of " + str(pid) + " are: (" + str(len(ancestors)) + " ancestor pairs)\n")
    else:
        fout.write("The ancestors of " + str(pid) + " are: (" + str(len(ancestors)) + " ancestors)\n")
    fout.write("\n".join(map(str, ancestors)))
    fout.close()
    return

def save_descendants(descendants, pid, output_file, pairs_flag):
    fout = open(output_file, mode="w", encoding="utf-8")
    if pairs_flag:
        fout.write("The descendant pairs (student, advisor) of " + str(pid) + " are: (" + str(len(descendants)) + " descendant pairs)\n")
    else:
        fout.write("The descendants of " + str(pid) + " are: (" + str(len(descendants)) + " descendants)\n")
    fout.write("\n".join(map(str, descendants)))
    fout.close()
    return

# Modules to be loaded to mgdb_entry.py script
def unary_search_ancestors(pid, neo4j_driver, output_file):

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
                    RETURN DISTINCT p1.pid, p1.name;
                """

    with neo4j_driver.session() as sess:
        nodes = sess.run(cypher_query)
        ancestors = [node.values() for node in nodes]

    save_ancestors(ancestors, pid, output_file, False)

    return

def binary_search_ancestors(pid, neo4j_driver, output_file):
    # Ref: https://stackoverflow.com/questions/16611723/using-multiple-match-clauses-doesnt-return-any-result-in-neo4j-cypher-query
    cypher_query = """
                    MATCH (p3:Person {pid: '""" + str(pid) + """'})-[*]->(d:Dissertation)
                    WITH DISTINCT d
                    MATCH (p2:Person)-[:WRITES]->(d)-[:ADVISED_BY]->(p1:Person)
                    RETURN DISTINCT p2.pid, p2.name, p1.pid, p1.name
                """
    with neo4j_driver.session() as sess:
        nodes = sess.run(cypher_query)
        ancestors = [node.values() for node in nodes]
    save_ancestors(ancestors, pid, output_file, True)

    return

def unary_search_descendants(pid, neo4j_driver, output_file):

    cypher_query = """
                    MATCH (d:Dissertation)-[*]->(p1:Person {pid: '""" + str(pid) + """'})
                    WITH DISTINCT d
                    MATCH (p2:Person)-[:WRITES]-(d)
                    RETURN DISTINCT p2.pid, p2.name;
                """

    with neo4j_driver.session() as sess:
        nodes = sess.run(cypher_query)
        descendants = [node.values() for node in nodes]
    
    save_descendants(descendants, pid, output_file, False)

    return

def binary_search_descendants(pid, neo4j_driver, output_file):

    cypher_query = """
                    MATCH (p1:Person)-[*]->(p3:Person {pid: '""" + str(pid) + """'})
                    WITH DISTINCT p1.pid AS de_pid
                    MATCH (p2:Person)-[:WRITES]->(d:Dissertation)-[:ADVISED_BY]->(p:Person)
                    WHERE p.pid=de_pid OR p.pid='""" + str(pid) + """'
                    RETURN DISTINCT p2.pid, p2.name, p.pid, p.name
                """
    with neo4j_driver.session() as sess:
        nodes = sess.run(cypher_query)
        descendants = [node.values() for node in nodes]
    
    save_descendants(descendants, pid, output_file, True)

    return

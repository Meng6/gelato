def save_ancestors(ancestors, pid, output_file):
    fout = open(output_file, mode="w", encoding="utf-8")
    fout.write("The ancestors of " + str(pid) + " are: (" + str(len(ancestors)) + " ancestors)\n")
    fout.write("\n".join(map(str, ancestors)))
    fout.close()
    return

def save_descendants(descendants, pid, output_file):
    fout = open(output_file, mode="w", encoding="utf-8")
    fout.write("The descendants of " + str(pid) + " are: (" + str(len(descendants)) + " descendants)\n")
    fout.write("\n".join(map(str, descendants)))
    fout.close()
    return

# Modules to be loaded to mgdb_entry.py script
def search_ancestors(pid, neo4j_driver, output_file):

    cypher_query = """
                    MATCH (p2:Person)-[*]->(d:Dissertation)-[:ADVISED_BY]->(p1:Person)
                    WHERE p2.pid='""" + str(pid) + """'
                    RETURN DISTINCT p1.pid, p1.name;
                """
    with neo4j_driver.session() as sess:
        nodes = sess.run(cypher_query)
        ancestors = [node.values() for node in nodes]

    save_ancestors(ancestors, pid, output_file)

    return

def search_descendants(pid, neo4j_driver, output_file):

    cypher_query = """
                    MATCH (p2:Person)-[:WRITES]->(d:Dissertation)-[*]->(p1:Person)
                    WHERE p1.pid='""" + str(pid) + """'
                    RETURN DISTINCT p2.pid, p2.name;
                """
    with neo4j_driver.session() as sess:
        nodes = sess.run(cypher_query)
        descendants = [node.values() for node in nodes]
    
    save_descendants(descendants, pid, output_file)

    return
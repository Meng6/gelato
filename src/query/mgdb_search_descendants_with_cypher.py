from neo4j import GraphDatabase
import yaml

def searchDescendantsCypher(pid, driver):

    cypher_query = """
                MATCH (p2:Person)-[:WRITES]->(d:Dissertation)-[*]->(p1:Person)
                WHERE p1.pid='""" + str(pid) + """'
                RETURN DISTINCT p2.pid, p2.name;
            """
    with driver.session() as sess:
        nodes = sess.run(cypher_query)
        nodes = [node.values() for node in nodes]
    
    return nodes

# Load credentials
database_group = snakemake.params["database_group"]
with open("./credentials.yaml", "r", encoding="utf-8") as f:
    credentials = yaml.safe_load(f)
# The URI should be in the form protocol://<server location>:<port>. 
# The supported protocols in URI could either be bolt or neo4j. 
# - bolt should be used when creating a driver connecting to the Neo4j instance directly. 
# - neo4j should be used when creating a driver with built-in routing.
uri = "{protocal}://{host}:{port}".format(protocal=credentials[database_group]["protocol"],
                                          host=credentials[database_group]["host"],
                                          port=credentials[database_group]["port"])
# Connect to Neo4j DB
neo4j_driver = GraphDatabase.driver(uri, auth=(credentials[database_group]["user"], credentials[database_group]["password"]))

# Search descendants
pid = snakemake.params["pid"]
students = searchDescendantsCypher(pid=pid, driver=neo4j_driver)
# Close the connection
neo4j_driver.close()

# Save the results
fout = open(snakemake.output[0], mode="w", encoding="utf-8")
fout.write("The descendants of " + str(pid) + " are: (" + str(len(students)) + " descendants)\n")
fout.write("\n".join(map(str, students)))
fout.close()

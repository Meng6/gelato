from neo4j import GraphDatabase
import yaml

def loadData2Neo4j(driver):
    
    # Create a constraint so that each Person/Dissertation node has a unique id propery.
    cypher_constraint1_query = """
        CREATE CONSTRAINT personPIdConstraint FOR (person:Person) REQUIRE person.pid IS UNIQUE;
        """
    cypher_constraint2_query = """
        CREATE CONSTRAINT dissertationDIdConstraint FOR (dissertation:Dissertation) REQUIRE dissertation.did IS UNIQUE;
        """

    # Load person.csv
    cypher_person_query = """
        LOAD CSV WITH HEADERS FROM "file:///tsv/person.tsv" AS csvLine FIELDTERMINATOR '\t'
        CREATE (p:Person {pid: csvLine.pid, name: csvLine.name, country: csvLine.country, onlinedescendants: csvLine.onlinedescendants});
        """
    
    # Load dissertation.csv
    cypher_dissertation_query = """
        LOAD CSV WITH HEADERS FROM "file:///tsv/dissertation.tsv" AS csvLine FIELDTERMINATOR '\t'
        CREATE (d:Dissertation {did: csvLine.did, pid: csvLine.author, title: csvLine.title, university: csvLine.university, year: csvLine.year});
        """
    
    # Add WRITES relationship
    cypher_writes_query = """
        MATCH (p:Person)
        MATCH (d:Dissertation)
        WHERE p.pid=d.pid
        MERGE (p)-[:WRITES]->(d);
        """
    
    # Add ADVISED_BY relationship
    cypher_advisedby_query = """
        LOAD CSV WITH HEADERS FROM "file:///tsv/advised.tsv" AS csvLine FIELDTERMINATOR '\t'
        MATCH (p:Person)
        MATCH (d:Dissertation)
        WHERE p.pid=csvLine.advisor and d.did=csvLine.did
        MERGE (d)-[:ADVISED_BY {advisororder: csvLine.advisororder}]->(p);
        """

    with driver.session() as sess:
        sess.run(cypher_constraint1_query)
        sess.run(cypher_constraint2_query)
        sess.run(cypher_person_query)
        sess.run(cypher_dissertation_query)
        sess.run(cypher_writes_query)
        sess.run(cypher_advisedby_query)
    
    return

# Load credentials
database_group = snakemake.params["database_group"]
with open("./credentials.yaml", "r", encoding="utf-8") as f:
    credentials = yaml.safe_load(f)
# The URI should be in the form protocol://<server location>:<port>. 
# The supported protocols in URI could either be bolt or neo4j. 
# - bolt should be used when creating a driver connecting to the Neo4j instance directly. 
# - neo4j should be used when creating a driver with built-in routing.
# Ref: https://neo4j.com/docs/api/python-driver/current/api.html
uri = "{protocal}://{host}:{port}?".format(protocal=credentials[database_group]["protocol"],
                                          host=credentials[database_group]["host"],
                                          port=credentials[database_group]["port"])
# Connect to Neo4j DB
neo4j_driver = GraphDatabase.driver(uri, auth=(credentials[database_group]["user"], credentials[database_group]["password"]))

# Search descendants
loadData2Neo4j(driver=neo4j_driver)
# Close the connection
neo4j_driver.close()

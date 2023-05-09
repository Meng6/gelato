from neo4j import GraphDatabase
import yaml

def loadData2Neo4j(driver):

    # Load dblp_v14_nlp.json file
    cypher_load_cn_data = """
    CALL apoc.load.jsonArray("file:///dblp_v14_nlp.json", "$")
    YIELD value
    WITH value
    MERGE (paper:Paper {nid: value.id})
    MERGE (venue:Venue {nid: value.venue.raw})
    MERGE (paper)-[:PUBLISHED_IN]->(venue)
    SET paper += {title: value.title, abstract: value.abstract, keywords: value.keywords, year: value.year, cntCitation: value.n_citation}

    WITH paper, value
    UNWIND value.authors AS authors
    WITH paper, value, authors.id AS authorIds
    UNWIND authorIds AS authorId
    MERGE (author:Author {nid: authorId})
    MERGE (author)-[:WROTE]->(paper)

    WITH paper, value, value.references AS refs
    WHERE refs IS NOT NULL
    UNWIND refs AS rid
    MERGE (ref:Paper {nid: rid})
    MERGE (paper)-[:CITED]->(ref);
    """

    with driver.session() as sess:
        sess.run(cypher_load_cn_data)
    
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

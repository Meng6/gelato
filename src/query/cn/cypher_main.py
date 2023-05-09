
# Modules to be loaded to cn_entry.py script
def research_question_extraction(params, neo4j_driver):

    keywords = params["keywords"]

    condition = ""
    for keyword in keywords:
        condition = condition + """
                                ((ANY(k in p1.keywords 
                                WHERE toLower(k) = "{keyword}")) 
                                OR (toLower(p1.abstract) CONTAINS "{keyword}")
                                OR (toLower(p1.title) CONTAINS "{keyword}"))
                                AND
                                ((ANY(k in p2.keywords 
                                WHERE toLower(k) = "{keyword}")) 
                                OR (toLower(p2.abstract) CONTAINS "{keyword}")
                                OR (toLower(p2.title) CONTAINS "{keyword}"))        
                            """.format(keyword=keyword)

    cypher_query = """
                    MATCH (p1:Paper)-[r]->(p2:Paper)
                    WHERE
                        {condition}
                    RETURN p1.nid, p2.nid;
                """.format(condition=condition)
    
    with neo4j_driver.session() as sess:
        extracted_research_question = sess.run(cypher_query).to_df()

    return extracted_research_question

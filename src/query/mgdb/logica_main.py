# Modules to be loaded to mgdb_entry.py script
def unary_search_ancestors(params, database):

    logica_query = """

    @Engine("sqlite");
    @AttachDatabase("mgdb", "{database}");
    @Dataset("advised");
    @Dataset("person");
    @Dataset("dissertation");

    #Anc(Advisor,Student)
    Adv_Stu(advisor:, student:author) :- Advised(did:x, advisor:),Dissertation(did:y, author:), x=y;
    Anc(advisor:, student:) :- Adv_Stu(advisor:, student:);
    Anc(advisor:x, author:m) :- Anc(advisor:x, student:y), Adv_Stu(advisor:z, student:m), y=z, m={pid};
    
    """.format(database=database, pid=params["pid"])

    output_name = "Anc"
    columns = ["advisor", "author"]

    return logica_query, output_name, columns
# Modules to be loaded to mgdb_entry.py script
def unary_search_ancestors(params, database):

    logica_query = """
    @Engine("sqlite");
    @AttachDatabase("mgdb", "{database}");
    @Dataset("advised");
    @Dataset("person");
    @Dataset("dissertation");

    Adv_Stu(advisor:, student:author) :- Advised(did:x, advisor:),Dissertation(did:y, author:), x=y;

    @Recursive(Anc,33);
    Anc(ancestor:advisor, student:m) distinct :- Adv_Stu(advisor:, student:m), m={pid};
    Anc(ancestor:x, student:l) distinct:- Adv_Stu(advisor:x, student:y),Anc(ancestor:z, student:l), y=z;

    Anc_with_Name(ancestor_id:x,ancestor_name:n) :-Anc(ancestor:x, student:l), person(pid:y,name:n),x=y;
    """.format(database=database, pid=params["pid"])

    output_name = "Anc_with_Name"
    columns = ["advisor", "author"]

    return logica_query, output_name, columns
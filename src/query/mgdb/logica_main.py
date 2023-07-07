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

    Anc_with_Name(ancestor_id:x,ancestor_name:n) :-Anc(ancestor:x, student:l), Person(pid:y,name:n),x=y;
        
    """.format(
        database=database, pid=params["pid"]
    )

    output_name = "Anc_with_Name"
    columns = ["advisor", "author"]

    return logica_query, output_name, columns


def binary_search_ancestors(params, database):
    logica_query = """

    @Engine("sqlite");
    @AttachDatabase("mgdb", "{database}");
    @Dataset("advised");
    @Dataset("person");
    @Dataset("dissertation");

    Adv_Stu(advisor:, student:author) :- Advised(did:x, advisor:),Dissertation(did:y, author:), x=y;

    @Recursive(Anc,33);
    Anc(ancestor1:advisor, student:m, ancestor2:m) distinct :- Adv_Stu(advisor:, student:m), m={pid};
    Anc(ancestor1:x, student:l, ancestor2:y) distinct:- Adv_Stu(advisor:x, student:y),Anc(ancestor1:z, student:l, ancestor2:), y=z;

    Anc_with_Name(student_id:y,student_name:name2, advisor_id:x,advisor_name:name1) :-Anc(ancestor1:x, student:l, ancestor2:y), Person(pid:pid1,name:name1),Person(pid:pid2,name:name2),x=pid1, y=pid2;        
    """.format(
        database=database, pid=params["pid"]
    )

    output_name = "Anc_with_Name"
    columns = ["student_id", "student_name", "advisor_id", "advisor_name"]

    return logica_query, output_name, columns


def unary_search_descendants(params, database):
    logica_query = """

    @Engine("sqlite");
    @AttachDatabase("mgdb", "{database}");
    @Dataset("advised");
    @Dataset("person");
    @Dataset("dissertation");

    Adv_Stu(advisor:, student:author) :- Advised(did:x, advisor:),Dissertation(did:y, author:), x=y;

    @Recursive(Desc,1);
    Desc(pid:m, descendant:student) distinct :- Adv_Stu(advisor:m, student:), m={pid};
    Desc(pid:, descendant:z) distinct:- Desc(pid:, descendant:l),Adv_Stu(advisor:y, student:z), l=y;

    Desc_with_Name(desc_id:l,desc_name:n) :-Desc(pid:, descendant:l), Person(pid:y,name:n),l=y;
    """.format(
        database=database, pid=params["pid"]
    )

    output_name = "Desc_with_Name"
    columns = ["desc_id", "desc_name"]

    return logica_query, output_name, columns


def binary_search_descendants(params, database):
    logica_query = """

    @Engine("sqlite");

    @AttachDatabase("mgdb","{database}");
    @Dataset("advised");
    @Dataset("person");
    @Dataset("dissertation");

    Adv_Stu(advisor:, student:author) :- Advised(did:x, advisor:),Dissertation(did:y, author:), x=y;

    Desc(pid:m, descendant_1:student, descendant_2:m) distinct :- Adv_Stu(advisor:m, student:), m={pid};
    Desc(pid:, descendant_1:z,descendant_2:y) distinct:- Desc(pid:, descendant_1:l),
    Adv_Stu(advisor:y, student:z), l=y;

    Desc_with_Name(student_id:z,student_name:name2,advisor_id:y,advisor_name:name1) :- 
    Desc(pid:, descendant_1:z,descendant_2:y), Person(pid:y,name:name1),Person(pid:z,name:name2);
    """.format(
        database=database, pid=params["pid"]
    )

    output_name = "Desc_with_Name"
    columns = ["student_id", "student_name", "advisor_id", "advisor_name"]

    return logica_query, output_name, columns
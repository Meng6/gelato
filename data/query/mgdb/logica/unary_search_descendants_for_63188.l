

    @Engine("sqlite");
    @AttachDatabase("mgdb", "data/external/mgdb/mgdb.db");
    @Dataset("advised");
    @Dataset("person");
    @Dataset("dissertation");

    Adv_Stu(advisor:, student:author) :- Advised(did:x, advisor:),Dissertation(did:y, author:), x=y;

    @Recursive(Desc,1);
    Desc(pid:m, descendant:student) distinct :- Adv_Stu(advisor:m, student:), m=63188;
    Desc(pid:, descendant:z) distinct:- Desc(pid:, descendant:l),Adv_Stu(advisor:y, student:z), l=y;

    Desc_with_Name(desc_id:l,desc_name:n) :-Desc(pid:, descendant:l), person(pid:y,name:n),l=y;
    
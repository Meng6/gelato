

    @Engine("sqlite");
    @AttachDatabase("mgdb", "data/external/mgdb/mgdb.db");
    @Dataset("advised");
    @Dataset("person");
    @Dataset("dissertation");

    Adv_Stu(advisor:, student:author) :- Advised(did:x, advisor:),Dissertation(did:y, author:), x=y;
    @Recursive(Anc,33);
    Anc(ancestor:advisor, student:m) distinct :- Adv_Stu(advisor:, student:m), m=63244;
    Anc(ancestor:x, student:l) distinct:- Adv_Stu(advisor:x, student:y),Anc(ancestor:z, student:l), y=z;

    Anc_with_Name(ancestor_id:x,ancestor_name:n) :-Anc(ancestor:x, student:l), person(pid:y,name:n),x=y;
        
    
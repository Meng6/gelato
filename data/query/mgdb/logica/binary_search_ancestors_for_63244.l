

    @Engine("sqlite");
    @AttachDatabase("mgdb", "data/external/mgdb/mgdb.db");
    @Dataset("advised");
    @Dataset("person");
    @Dataset("dissertation");

    Adv_Stu(advisor:, student:author) :- Advised(did:x, advisor:),Dissertation(did:y, author:), x=y;

    @Recursive(Anc,33);
    Anc(ancestor1:advisor, student:m, ancestor2:m) distinct :- Adv_Stu(advisor:, student:m), m=63244;
    Anc(ancestor1:x, student:l, ancestor2:y) distinct:- Adv_Stu(advisor:x, student:y),Anc(ancestor1:z, student:l, ancestor2:), y=z;

    Anc_with_Name(student_id:y,student_name:name2, advisor_id:x,advisor_name:name1) :-Anc(ancestor1:x, student:l, ancestor2:y), person(pid:pid1,name:name1),person(pid:pid2,name:name2),x=pid1, y=pid2;        
    
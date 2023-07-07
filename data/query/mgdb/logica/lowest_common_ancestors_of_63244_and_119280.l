

    @Engine("sqlite");

    @AttachDatabase("mgdb","data/external/mgdb/mgdb.db");
    @Dataset("advised");
    @Dataset("person");
    @Dataset("dissertation");

    Adv_Stu(advisor:, student:author) :- Advised(did:x, advisor:),Dissertation(did:y, author:), x=y;

    @Recursive(Anc_1,33);
    Anc_1(ancestor:advisor, student:m) distinct :- Adv_Stu(advisor:, student:m), m=63244;
    Anc_1(ancestor:x, student:) distinct:- Adv_Stu(advisor:x, student:y),Anc_1(ancestor:y, student:);

    @Recursive(Anc_2,10);
    Anc_2(ancestor:advisor, student:l) distinct :- Adv_Stu(advisor:, student:l), l=119280;
    Anc_2(ancestor:x, student:) distinct:- Adv_Stu(advisor:x, student:y),Anc_2(ancestor:y, student:);

    Common_Ancestors(pid_1:m, pid_2:l, anc_id: x) :- Anc_1(ancestor:x, student:m), 
    Anc_2(ancestor:x, student:l), m=63244,l=119280 ;

    Not_Lowest_Common_Ancestors(pid_1:m, pid_2:l, anc_id: x) :- Common_Ancestors(pid_1:m, pid_2:l, anc_id: x), 
    Common_Ancestors(pid_1:m, pid_2:l, anc_id: y), Adv_Stu(advisor:x, student:y),m=63244,l=119280;

    Lowest_Common_Ancestors_with_Name(pid:x, name:name) :- Common_Ancestors(pid_1:m, pid_2:l, anc_id: x),
    ~Not_Lowest_Common_Ancestors(pid_1:m, pid_2:l, anc_id: x), Person(pid:x,name:name),m=63244,l=119280;
    
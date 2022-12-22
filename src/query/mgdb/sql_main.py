def save_ancestors(ancestors, pid, output_file, pairs_flag):
    fout = open(output_file, mode="w", encoding="utf-8")
    if pairs_flag:
        fout.write("The ancestor pairs (student, advisor) of " + str(pid) + " are: (" + str(len(ancestors)) + " ancestor pairs)\n")
    else:
        fout.write("The ancestors of " + str(pid) + " are: (" + str(len(ancestors)) + " ancestors)\n")
    fout.write("\n".join(map(str, ancestors)))
    fout.close()
    return

def save_descendants(descendants, pid, output_file, pairs_flag):
    fout = open(output_file, mode="w", encoding="utf-8")
    if pairs_flag:
        fout.write("The descendant pairs (student, advisor) of " + str(pid) + " are: (" + str(len(descendants)) + " descendant pairs)\n")
    else:
        fout.write("The descendants of " + str(pid) + " are: (" + str(len(descendants)) + " descendants)\n")
    fout.write("\n".join(map(str, descendants)))
    fout.close()
    return

def save_lowest_common_ancestors(lowest_common_ancestors, pid1, pid2, output_file):
    fout = open(output_file, mode="w", encoding="utf-8")
    fout.write("The lowest common ancestor(s) of " + str(pid1) + " and " + str(pid2) + ": \n")
    fout.write("\n".join(map(str, lowest_common_ancestors)))
    fout.close()
    return

# Modules to be loaded to mgdb_entry.py script
def unary_search_ancestors(params, conn, output_file):
    pid = params["pid"]
    cursor = conn.cursor()
    
    # Recursive CTE
    sql_query = """
    WITH RECURSIVE ancestors AS (
	    SELECT advisor FROM advised a, dissertation d WHERE a.did=d.did AND author = """ + str(pid) + """
	    UNION
	    SELECT a.advisor FROM ancestors an, advised a, dissertation d WHERE a.did=d.did AND d.author=an.advisor
	    )
    SELECT advisor, name FROM ancestors LEFT JOIN person p ON p.pid = advisor
    """

    cursor.execute(sql_query)
    ancestors = cursor.fetchall()
    cursor.close()
    
    save_ancestors(ancestors, pid, output_file, False)
    return

def binary_search_ancestors(params, conn, output_file):
    pid = params["pid"]
    cursor = conn.cursor()
    
    # Recursive CTE
    sql_query = """
    WITH RECURSIVE ancestors AS (
    	SELECT author, advisor FROM advised a, dissertation d WHERE a.did=d.did AND author = """ + str(pid) + """
    	UNION
    	SELECT d.author, a.advisor FROM ancestors an, advised a, dissertation d WHERE a.did=d.did AND d.author=an.advisor
        )
    SELECT author, p1.name, advisor, p2.name
    FROM
    	ancestors 
    	LEFT JOIN person p1 ON p1.pid = author
    	LEFT JOIN person p2 ON p2.pid = advisor
    """

    cursor.execute(sql_query)
    ancestors = cursor.fetchall()
    cursor.close()
    
    save_ancestors(ancestors, pid, output_file, True)
    return

def unary_search_descendants(params, conn, output_file):
    pid = params["pid"]
    cursor = conn.cursor()
    
    # Recursive CTE
    sql_query = """
    WITH RECURSIVE descendants AS (
        SELECT author FROM dissertation d, advised a WHERE a.did=d.did AND advisor = """ + str(pid) + """
        UNION
        SELECT d.author FROM dissertation d, advised a, descendants de WHERE a.did=d.did AND a.advisor=de.author
        )
    SELECT author, name FROM descendants LEFT JOIN person p ON p.pid = author
    """

    cursor.execute(sql_query)
    descendants = cursor.fetchall()
    cursor.close()

    save_descendants(descendants, pid, output_file, False)
    return

def binary_search_descendants(params, conn, output_file):
    pid = params["pid"]
    cursor = conn.cursor()
    
    # Recursive CTE
    sql_query = """
    WITH RECURSIVE descendants AS (
        SELECT author, advisor FROM dissertation d, advised a WHERE a.did=d.did AND advisor = """ + str(pid) + """
        UNION
        SELECT d.author, a.advisor FROM dissertation d, advised a, descendants de WHERE a.did=d.did AND a.advisor=de.author
        )
    SELECT author, p1.name, advisor, p2.name
    FROM
    	descendants 
    	LEFT JOIN person p1 ON p1.pid = author
    	LEFT JOIN person p2 ON p2.pid = advisor
    """

    cursor.execute(sql_query)
    descendants = cursor.fetchall()
    cursor.close()

    save_descendants(descendants, pid, output_file, True)
    return

def lowest_common_ancestors(params, conn, output_file):
    pid1 = params["pid1"]
    pid2 = params["pid2"]
    cursor = conn.cursor()

    # Recursive CTE: ancestors of pid1
    sql_query1 = """
    CREATE VIEW ancestors_of_pid1 AS 
        WITH RECURSIVE ancestors AS (
            SELECT advisor FROM advised a, dissertation d WHERE a.did=d.did AND author = """ + str(pid1) + """
            UNION
            SELECT a.advisor FROM ancestors an, advised a, dissertation d WHERE a.did=d.did AND d.author=an.advisor
            )
        SELECT advisor FROM ancestors;
    """
    cursor.execute(sql_query1)
    
    # Recursive CTE: ancestors of pid2
    sql_query2 = """
    CREATE VIEW ancestors_of_pid2 AS 
        WITH RECURSIVE ancestors AS (
            SELECT advisor FROM advised a, dissertation d WHERE a.did=d.did AND author = """ + str(pid2) + """
            UNION
            SELECT a.advisor FROM ancestors an, advised a, dissertation d WHERE a.did=d.did AND d.author=an.advisor
            )
        SELECT advisor FROM ancestors;
    """
    cursor.execute(sql_query2)

    # Common ancestors of pid1 and pid2
    sql_query3 = """
    CREATE VIEW common_ancestors AS
    SELECT advisor FROM ancestors_of_pid1
    INTERSECT
    SELECT advisor FROM ancestors_of_pid2;
    """
    cursor.execute(sql_query3)


    # Common ancestors with their closest descendants
    sql_query4 = """
    CREATE VIEW common_ancestors_with_their_students AS
    SELECT d.author, ca.advisor
    FROM common_ancestors ca, advised a, dissertation d 
    WHERE ca.advisor=a.advisor AND a.did=d.did;
    """
    cursor.execute(sql_query4)

    # Lowest common ancestors of pid1 and pid2
    sql_query5 = """

    SELECT advisor, name FROM
    
        (SELECT DISTINCT advisor FROM common_ancestors_with_their_students

        EXCEPT

        SELECT DISTINCT advisor FROM common_ancestors_with_their_students
        WHERE author IN (SELECT advisor FROM common_ancestors))

    LEFT JOIN person p ON p.pid = advisor
    """
    cursor.execute(sql_query5)
    lowest_common_ancestors = cursor.fetchall()

    # Drop the views
    cursor.execute("DROP VIEW ancestors_of_pid1;")
    cursor.execute("DROP VIEW ancestors_of_pid2;")
    cursor.execute("DROP VIEW common_ancestors;")
    cursor.execute("DROP VIEW common_ancestors_with_their_students;")
    cursor.close()

    save_lowest_common_ancestors(lowest_common_ancestors, pid1, pid2, output_file)
    return

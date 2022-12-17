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

# Modules to be loaded to mgdb_entry.py script
def unary_search_ancestors(pid, conn, output_file):

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

def binary_search_ancestors(pid, conn, output_file):

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

def unary_search_descendants(pid, conn, output_file):
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

def binary_search_descendants(pid, conn, output_file):
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
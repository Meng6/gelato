def save_ancestors(ancestors, pid, output_file):
    fout = open(output_file, mode="w", encoding="utf-8")
    fout.write("The ancestors of " + str(pid) + " are: (" + str(len(ancestors)) + " ancestors)\n")
    fout.write("\n".join(map(str, ancestors)))
    fout.close()
    return

def save_descendants(descendants, pid, output_file):
    fout = open(output_file, mode="w", encoding="utf-8")
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
	    select advisor from advised a, dissertation d where a.did=d.did and author = """ + str(pid) + """
	    UNION
	    select a.advisor from ancestors an, advised a, dissertation d where a.did=d.did and d.author=an.advisor
	    )
    SELECT advisor, name from ancestors left join person p on p.pid = advisor
    """

    cursor.execute(sql_query)
    ancestors = cursor.fetchall()
    cursor.close()
    
    save_ancestors(ancestors, pid, output_file)
    return

def unary_search_descendants(pid, conn, output_file):
    cursor = conn.cursor()
    
    # Recursive CTE
    sql_query = """
    WITH RECURSIVE descendants AS (
        select author from dissertation d, advised a where a.did=d.did and advisor = """ + str(pid) + """
        UNION
        select di.author from dissertation di, advised a, descendants d where a.did=di.did and a.advisor=d.author
        )
    SELECT author, name from descendants left join person p on p.pid = author
    """

    cursor.execute(sql_query)
    descendants = cursor.fetchall()
    cursor.close()

    save_descendants(descendants, pid, output_file)
    return
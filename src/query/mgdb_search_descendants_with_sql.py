import sqlite3

def searchDescendantsSQL(pid):

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
    students = cursor.fetchall()
    cursor.close()
    return students

# Connect to SQLite DB
conn = sqlite3.connect(snakemake.input["database"])

# Search descendants
pid = snakemake.params["pid"]
students = searchDescendantsSQL(pid=pid)
# Close the connection
conn.close ()

# Save the results
fout = open(snakemake.output[0], mode="w", encoding="utf-8")
fout.write("The descendants of " + str(pid) + " are: (" + str(len(students)) + " descendants)\n")
fout.write("\n".join(map(str, students)))
fout.close()
import sqlite3

def searchAncestorsSQL(pid):

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
    return ancestors

# Connect to SQLite DB
conn = sqlite3.connect(snakemake.input["database"])

# Search ancestors
pid = snakemake.params["pid"]
ancestors = searchAncestorsSQL(pid=pid)
# Close the connection
conn.close ()

# Save the results
fout = open(snakemake.output[0], mode="w", encoding="utf-8")
fout.write("The ancestors of " + str(pid) + " are: (" + str(len(ancestors)) + " ancestors)\n")
fout.write("\n".join(map(str, ancestors)))
fout.close()
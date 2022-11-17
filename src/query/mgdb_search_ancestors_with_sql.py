import psycopg2, yaml

def searchAncestorsSQL(pid):

    cursor = conn.cursor()
    
    # Recursive CTE
    sql_query = """
    WITH RECURSIVE ancestors AS (
	    ( select advisor from advised a, dissertation d where a.did=d.did and author = """ + str(pid) + """ )
	    UNION (
	    select a.advisor from ancestors an, advised a, dissertation d where a.did=d.did and d.author=an.advisor
	    ))
    SELECT advisor, name from ancestors left join person p on p.pid = advisor
    """

    cursor.execute(sql_query)
    ancestors = cursor.fetchall()
    cursor.close()
    return ancestors

# Load credentials
database_group = snakemake.params["database_group"]
with open("./credentials.yaml", "r", encoding="utf-8") as f:
    credentials = yaml.safe_load(f)
# Connect to PostgreSQL DB
conn = psycopg2.connect(database=credentials[database_group]["database"],
                        host=credentials[database_group]["host"],
                        user=credentials[database_group]["user"],
                        password=credentials[database_group]["password"],
                        port=credentials[database_group]["port"])


pid = snakemake.params["pid"]
ancestors = searchAncestorsSQL(pid=pid)

# Save the results
fout = open(snakemake.output[0], mode="w", encoding="utf-8")
fout.write("The ancestors of " + str(pid) + " are: (" + str(len(ancestors)) + " ancestors)\n")
fout.write("\n".join(map(str, ancestors)))
# Close the connection
conn.close ()
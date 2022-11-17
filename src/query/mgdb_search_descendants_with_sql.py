import psycopg2, yaml

def searchDescendantsSQL(pid):

    cursor = conn.cursor()
    
    # Recursive CTE
    sql_query = """
    WITH RECURSIVE descendants AS (
        ( select author from dissertation d, advised a where a.did=d.did and advisor = """ + str(pid) + """)
        UNION (
        select di.author from dissertation di, advised a, descendants d where a.did=di.did and a.advisor=d.author
        ))
    SELECT author, name from descendants left join person p on p.pid = author
    """

    cursor.execute(sql_query)
    students = cursor.fetchall()
    cursor.close()
    return students

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
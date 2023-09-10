import sys
sys.path.append('..')
from tools.template import timeit

# Modules to be loaded to mgdb_entry.py script
@timeit
def unary_search_ancestors(params, conn, **kwargs):
    """
    Query unary_search_ancestors: extract ancestors of a researcher with the given PID.
    
    Args:
    - params (dict) : PID could be accessed via params["pid"]
    - conn (SQLite connection object) : connection to the SQLite database

    Returns:
    - ancestors (list of tuples): result set, i.e., [(pid, name), (pid, name), ...]
    - columns (list of strings): column names of the result set, i.e., pid, name
    """

    pid = params["pid"]
    cursor = conn.cursor()
    
    # Recursive CTE
    sql_query = """
    WITH RECURSIVE ancestors AS (
	    SELECT advisor FROM advised a, dissertation d WHERE a.did=d.did AND author = """ + str(pid) + """
	    UNION
	    SELECT a.advisor FROM ancestors an, advised a, dissertation d WHERE a.did=d.did AND d.author=an.advisor
	    )
    SELECT advisor AS pid, name FROM ancestors LEFT JOIN person p ON p.pid = advisor
    """

    cursor.execute(sql_query)
    ancestors = cursor.fetchall()
    columns = [x[0] for x in cursor.description]
    cursor.close()
    
    return ancestors, columns

@timeit
def binary_search_ancestors(params, conn, **kwargs):
    """
    Query binary_search_ancestors: extract ancestor pairs of a researcher with the given PID.
    
    Args:
    - params (dict) : PID could be accessed via params["pid"]
    - conn (SQLite connection object) : connection to the SQLite database

    Returns:
    - ancestors (list of tuples): result set, i.e., [(student_pid, student_name, advisor_pid, advisor_name), ...]
    - columns (list of strings): column names of the result set, i.e., student_pid, student_name, advisor_pid, advisor_name
    """

    pid = params["pid"]
    cursor = conn.cursor()
    
    # Recursive CTE
    sql_query = """
    WITH RECURSIVE ancestors AS (
    	SELECT author, advisor FROM advised a, dissertation d WHERE a.did=d.did AND author = """ + str(pid) + """
    	UNION
    	SELECT d.author, a.advisor FROM ancestors an, advised a, dissertation d WHERE a.did=d.did AND d.author=an.advisor
        )
    SELECT author AS student_pid, p1.name AS student_name, advisor AS advisor_pid, p2.name AS advisor_name
    FROM
    	ancestors 
    	LEFT JOIN person p1 ON p1.pid = author
    	LEFT JOIN person p2 ON p2.pid = advisor
    """

    cursor.execute(sql_query)
    ancestors = cursor.fetchall()
    columns = [x[0] for x in cursor.description]
    cursor.close()
    
    return ancestors, columns

@timeit
def unary_search_descendants(params, conn, **kwargs):
    """
    Query unary_search_descendants: extract descendants of a researcher with the given PID.
    
    Args:
    - params (dict) : PID could be accessed via params["pid"]
    - conn (SQLite connection object) : connection to the SQLite database

    Returns:
    - descendants (list of tuples): result set, i.e., [(pid, name), (pid, name), ...]
    - columns (list of strings): column names of the result set, i.e., pid, name
    """

    pid = params["pid"]
    cursor = conn.cursor()
    
    # Recursive CTE
    sql_query = """
    WITH RECURSIVE descendants AS (
        SELECT author FROM dissertation d, advised a WHERE a.did=d.did AND advisor = """ + str(pid) + """
        UNION
        SELECT d.author FROM dissertation d, advised a, descendants de WHERE a.did=d.did AND a.advisor=de.author
        )
    SELECT author AS pid, name FROM descendants LEFT JOIN person p ON p.pid = author
    """

    cursor.execute(sql_query)
    descendants = cursor.fetchall()
    columns = [x[0] for x in cursor.description]
    cursor.close()

    return descendants, columns

@timeit
def binary_search_descendants(params, conn, **kwargs):
    """
    Query binary_search_descendants: extract descendant pairs of a researcher with the given PID.
    
    Args:
    - params (dict) : PID could be accessed via params["pid"]
    - conn (SQLite connection object) : connection to the SQLite database

    Returns:
    - descendants (list of tuples): result set, i.e., [(student_pid, student_name, advisor_pid, advisor_name), ...]
    - columns (list of strings): column names of the result set, i.e., student_pid, student_name, advisor_pid, advisor_name
    """

    pid = params["pid"]
    cursor = conn.cursor()
    
    # Recursive CTE
    sql_query = """
    WITH RECURSIVE descendants AS (
        SELECT author, advisor FROM dissertation d, advised a WHERE a.did=d.did AND advisor = """ + str(pid) + """
        UNION
        SELECT d.author, a.advisor FROM dissertation d, advised a, descendants de WHERE a.did=d.did AND a.advisor=de.author
        )
    SELECT author AS student_pid, p1.name AS student_name, advisor AS advisor_pid, p2.name AS advisor_name
    FROM
    	descendants 
    	LEFT JOIN person p1 ON p1.pid = author
    	LEFT JOIN person p2 ON p2.pid = advisor
    """

    cursor.execute(sql_query)
    descendants = cursor.fetchall()
    columns = [x[0] for x in cursor.description]
    cursor.close()

    return descendants, columns

@timeit
def lowest_common_ancestors(params, conn, **kwargs):
    """
    Query lowest_common_ancestors: extract LCAs of a pair of researchers with the given PIDS, i.e., (PID1, PID2).
    
    Args:
    - params (dict) : the (PID1, PID2) pair could be accessed via params["pid1"] and params["pid2"]
    - conn (SQLite connection object) : connection to the SQLite database

    Returns:
    - lowest_common_ancestors (list of tuples): result set, i.e., [(pid, name), (pid, name), ...]
    - columns (list of strings): column names of the result set, i.e., pid, name
    """

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

    SELECT advisor AS pid, name FROM
    
        (SELECT DISTINCT advisor FROM common_ancestors_with_their_students

        EXCEPT

        SELECT DISTINCT advisor FROM common_ancestors_with_their_students
        WHERE author IN (SELECT advisor FROM common_ancestors))

    LEFT JOIN person p ON p.pid = advisor
    """
    cursor.execute(sql_query5)
    lowest_common_ancestors = cursor.fetchall()
    columns = [x[0] for x in cursor.description]

    # Drop the views
    cursor.execute("DROP VIEW ancestors_of_pid1;")
    cursor.execute("DROP VIEW ancestors_of_pid2;")
    cursor.execute("DROP VIEW common_ancestors;")
    cursor.execute("DROP VIEW common_ancestors_with_their_students;")
    cursor.close()

    return lowest_common_ancestors, columns

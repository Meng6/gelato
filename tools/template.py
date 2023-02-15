import pandas as pd
from io import StringIO

def format_output(data, columns, lat):

    formated_data = None

    """
    format_output() function converts the output of any language or tool into the same format: Pandas dataframe.
    :param data: output data
    :param columns: column names
    :param lat: language or tool
    :return: a Pandas dataframe
    """

    if lat == "python" or lat == "cypher":
        formated_data = data
    elif lat == "sql":
        formated_data = pd.DataFrame(data=data, columns=columns)
    elif lat == "sparql":
        formated_data = pd.DataFrame(data).applymap(lambda x: x["value"].replace("http://MGDB.com/p", ""))
    elif lat == "clingo":
        formated_data = pd.read_csv(StringIO("\n".join([x.strip()[1:-1] for x in data])), names=columns)
    else:
        raise ValueError("lat parameter only supports [python, sql, cypher, sparql, clingo] for now.")

    return formated_data
    
    

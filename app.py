import sys, inspect, os, subprocess, yaml, ast
sys.path.append('.')
import streamlit as st
import pandas as pd

code_block = """
import sys
sys.path.append('..')
from tools.template import timeit
import pandas as pd
from SPARQLWrapper import JSON
"""

def ensure_dir_exists(file_path):
    dir_path = os.path.split(file_path)[0]
    if (dir_path != '') and (not os.path.isdir(dir_path)):
        os.makedirs(dir_path)
    return

def interactive_page(query_name, query_function, lat, **kwargs):
    file_path = "src/query/mgdb/app/{lat}_main.py".format(lat=lat)
    with st.expander("See default {lat} code".format(lat=lat)):
        st.code(inspect.getsource(query_function), line_numbers=True)
    with st.expander("Check your {lat} code".format(lat=lat)):
        if os.path.isfile(file_path):
            with open(file_path, "r", encoding="utf-8") as f:
                st.code(f.read(), language="python", line_numbers=True)
    with st.expander("Write your {lat} code".format(lat=lat)):
        st.info(inspect.getdoc(query_function))
        st.text_area("Your {lat} code".format(lat=lat), value="@timeit\ndef {query_name}{variables}:".format(query_name=query_name, variables=str(inspect.signature(query_function))), key="{lat}_code".format(lat=lat), height=300)
        on = st.checkbox("Overwrite {lat} code".format(lat=lat))
        user_code = st.session_state["{lat}_code".format(lat=lat)]
        if on:
            ensure_dir_exists(file_path)
            with open(file_path, "w", encoding="utf-8") as f:
                f.write(code_block+"\n"+user_code+"\n")
        else:
            if os.path.isfile(file_path):
                with open(file_path, "a", encoding="utf-8") as f:
                    f.write("\n"+user_code+"\n")
            else:
                ensure_dir_exists(file_path)
                with open(file_path, "w") as f:
                    f.write(code_block+"\n"+user_code+"\n")
        if st.button("Run"):
            with open("app/config.yaml", "r", encoding="utf-8") as f:
                content = f.read()\
                            .replace("LANGUAGES_AND_TOOLS: {lat}".format(lat=config["MGDB"]["LANGUAGES_AND_TOOLS"]), "LANGUAGES_AND_TOOLS: [{lat}]".format(lat=lat.upper()))\
                            .replace("RUN: {run}", "RUN: [{query_name}]".format(run=config["MGDB"]["QUERIES"]["RUN"], query_name=query_name.upper()))
            if "pids" in kwargs:
                content = content.replace("    {query_name}:\n      PIDS: {pids}".format(query_name=query_name.upper(), pids=config["MGDB"]["QUERIES"][query_name.upper()]["PIDS"]),
                                    "    {query_name}:\n      PIDS: {pids}".format(query_name=query_name.upper(), pids=kwargs["pids"]))
            else: # pid1, pid2
                content = content.replace("    {query_name}:\n      PID1: {pid1}".format(query_name=query_name.upper(), pid1=config["MGDB"]["QUERIES"][query_name.upper()]["PID1"]),
                                    "    {query_name}:\n      PID1: {pid1}".format(query_name=query_name.upper(), pid1=kwargs["pid1"]))\
                                 .replace("    {query_name}:\n      PID2: {pid2}".format(query_name=query_name.upper(), pid2=config["MGDB"]["QUERIES"][query_name.upper()]["PID2"]),
                                    "    {query_name}:\n      PID2: {pid2}".format(query_name=query_name.upper(), pid2=kwargs["pid2"]))
            with open("app/config.yaml", "w", encoding="utf-8") as f:
                f.write(content)

            rule_name = "mgdb_{query_name}_with_{lat}".format(query_name=query_name, lat=lat)
            subprocess.run(["./gelato", "--snakefile", "app/Snakefile", "--configfile", "app/config.yaml", "--allowed-rules", rule_name, "-R", rule_name, "-j1"])
            if "pids" in kwargs:
                for pid in kwargs["pids"]:
                    st.write("**PID: {pid}**".format(pid=pid))
                    st.write(pd.read_csv("data/query/mgdb/{lat}/output_{query_name}_for_{pid}.csv".format(lat=lat, query_name=query_name, pid=pid)))
            else: # pid1, pid2
                st.write("**PID1: {pid1}, PID2: {pid2}**".format(pid1=kwargs["pid1"], pid2=kwargs["pid2"]))
                st.write(pd.read_csv("data/query/mgdb/{lat}/output_{query_name}_of_{pid1}_and_{pid2}.csv".format(lat=lat, query_name=query_name, pid1=kwargs["pid1"], pid2=kwargs["pid2"])))

def select_lat(query_name, **kwargs):
    if lat == "PYTHON":
        from src.query.mgdb import python_main
        query_function = getattr(python_main, query_name)
        interactive_page(query_name, query_function, "python", **kwargs)
    if lat == "SQL":
        from src.query.mgdb import sql_main
        query_function = getattr(sql_main, query_name)
        interactive_page(query_name, query_function, "sql", **kwargs)
    if lat == "CYPHER":
        from src.query.mgdb import cypher_main
        query_function = getattr(cypher_main, query_name)
        interactive_page(query_name,query_function, "cypher", **kwargs)
    if lat == "SPARQL":
        from src.query.mgdb import sparql_main
        query_function = getattr(sparql_main, query_name)
        interactive_page(query_name, query_function, "sparql", **kwargs)
    if lat == "CLINGO":
        from src.query.mgdb import clingo_main
        query_function = getattr(clingo_main, query_name)
        interactive_page(query_name, query_function, "clingo", **kwargs)
    if lat == "BASHLOG":
        from src.query.mgdb import bashlog_main
        query_function = getattr(bashlog_main, query_name)
        interactive_page(query_name, query_function, "bashlog", **kwargs)
    if lat == "LOGICA":
        from src.query.mgdb import logica_main
        query_function = getattr(logica_main, query_name)
        interactive_page(query_name, query_function, "logica", **kwargs)
    return

with st.sidebar:
    st.title("GELATO")
    query_name = st.radio("Query", ("UNARY_SEARCH_ANCESTORS", "BINARY_SEARCH_ANCESTORS", "UNARY_SEARCH_DESCENDANTS", "BINARY_SEARCH_DESCENDANTS", "LOWEST_COMMON_ANCESTORS"))


with open("app/config.yaml", "r", encoding="utf-8") as f:
    config = yaml.safe_load(f)

if query_name in ["UNARY_SEARCH_ANCESTORS", "BINARY_SEARCH_ANCESTORS", "UNARY_SEARCH_DESCENDANTS", "BINARY_SEARCH_DESCENDANTS"]:
    with st.sidebar:
        st.text_input("PIDS", config["MGDB"]["QUERIES"][query_name]["PIDS"], key="pids")
        lat = st.radio("Language or Tool", ("PYTHON", "SQL", "CYPHER", "SPARQL", "CLINGO", "BASHLOG", "LOGICA"))
    select_lat(query_name.lower(), pids=ast.literal_eval(st.session_state.pids))

if query_name == "LOWEST_COMMON_ANCESTORS":
    with st.sidebar:
        st.text_input("PID1", config["MGDB"]["QUERIES"][query_name]["PID1"], key="pid1")
        st.text_input("PID2", config["MGDB"]["QUERIES"][query_name]["PID2"], key="pid2")
        lat = st.radio("Language or Tool", ("PYTHON", "SQL", "CYPHER", "SPARQL", "CLINGO", "BASHLOG", "LOGICA"))
    select_lat(query_name.lower(), pid1=st.session_state.pid1, pid2=st.session_state.pid2)


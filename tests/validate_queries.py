import pandas as pd
import yaml

with open("config.yaml", "r", encoding="utf-8") as f:
    config = yaml.safe_load(f)
graphs = config["GRAPHS"]

for graph in graphs:
    lats = config[graph]["LANGUAGES_AND_TOOLS"]
    if len(lats) <= 1:
        print("WARNING: only have one language or tool for {graph} graph. Cannot validate the results.".format(graph=graph))
    else:
        for query in config[graph]["QUERIES"]["RUN"]:
            if query in ["UNARY_SEARCH_ANCESTORS", "BINARY_SEARCH_ANCESTORS", "UNARY_SEARCH_DESCENDANTS", "BINARY_SEARCH_DESCENDANTS"]:
                for pid in config[graph]["QUERIES"][query]["PIDS"]:
                    base_data, cols, base_lat = None, [], None
                    for lat in lats:
                        curr_data = pd.read_csv("data/query/{graph}/{lat}/output_{query}_for_{pid}.txt".format(graph=graph, lat=lat, query=query.lower(), pid=pid), sep=",")
                        if base_data is None:
                            cols = curr_data.columns.tolist()
                            base_data = curr_data.sort_values(by=cols).reset_index(drop=True).copy()
                        else:
                            curr_data = curr_data.sort_values(by=cols).reset_index()[cols]
                            if not curr_data.equals(base_data):
                                raise ValueError("TEST FAILED! {graph} graph: results of {base_lat} and {lat} for {query} are not consistent.".format(graph=graph, base_lat=base_lat, lat=lat, query=query))
                        base_lat = lat
            if query == "LOWEST_COMMON_ANCESTORS":
                pid1, pid2 = config[graph]["QUERIES"][query]["PID1"], config[graph]["QUERIES"][query]["PID2"]
                base_data, cols, base_lat = None, [], None
                for lat in lats:
                    curr_data = pd.read_csv("data/query/{graph}/{lat}/output_{query}_of_{pid1}_and_{pid2}.txt".format(graph=graph, lat=lat, query=query.lower(), pid1=pid1, pid2=pid2), sep=",")
                    if base_data is None:
                        cols = curr_data.columns.tolist()
                        base_data = curr_data.sort_values(by=cols).reset_index(drop=True).copy()
                    else:
                        curr_data = curr_data.sort_values(by=cols).reset_index()[cols]
                        if not curr_data.equals(base_data):
                            raise ValueError("TEST FAILED! {graph} graph: results of {base_lat} and {lat} for {query} are not consistent.".format(graph=graph, base_lat=base_lat, lat=lat, query=query))
                    base_lat = lat
        print("TEST PASSED! {graph} graph done.".format(graph=graph))
print("ALL TEST PASSED!")
                        

            


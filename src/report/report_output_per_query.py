import pandas as pd
import sys
sys.path.append('.')
from tools.template import report_output_per_query

formated_data = pd.read_csv(snakemake.input["formated_data"], sep=",")

html = report_output_per_query(formated_data=formated_data,
                        title=snakemake.params["query"],
                        edge=snakemake.params["edge"],
                        include_interactive_table=snakemake.params["include_interactive_table"],
                        include_network_graph=snakemake.params["include_network_graph"])

html_file = open(snakemake.output[0], "w", encoding="utf-8")
html_file.write(html)
html_file.close()

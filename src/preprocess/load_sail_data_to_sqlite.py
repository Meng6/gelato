import pandas as pd
import sqlite3

# Create a database
conn = sqlite3.connect(snakemake.output[0])

# Read TSV files
sail = pd.read_csv(snakemake.input["sail"], sep="\t")

# Load data into the database
sail.to_sql("sail_" + snakemake.params["max_hamming_number"], conn, if_exists="replace", index = False)

conn.commit()
conn.close()

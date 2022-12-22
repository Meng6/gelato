import pandas as pd
import sqlite3

# Create a database
conn = sqlite3.connect(snakemake.output[0])

# Read TSV files
fish = pd.read_csv(snakemake.input["fish"], sep="\t")

# Load data into the database
fish.to_sql("fish_" + snakemake.params["max_hamming_number"], conn, if_exists="replace", index = False)

conn.commit()
conn.close()

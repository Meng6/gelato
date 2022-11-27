import pandas as pd
import sqlite3

# Create a database
conn = sqlite3.connect(snakemake.output[0])

# Read TSV files
advised = pd.read_csv(snakemake.input["advised"], sep="\t")
dissertation = pd.read_csv(snakemake.input["dissertation"], sep="\t")
person = pd.read_csv(snakemake.input["person"], sep="\t")

# Load data into the database
advised.to_sql("advised", conn, if_exists="replace", index = False)
dissertation.to_sql("dissertation", conn, if_exists="replace", index = False)
person.to_sql("person", conn, if_exists="replace", index = False)

conn.commit()
conn.close()
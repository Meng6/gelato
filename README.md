# GELATO

**G**raph **E**xploration: **L**anguages **a**nd **To**ols

### Setup 

#### Installation (For Mac)
1. Install [brew](https://brew.sh/)

2. Install MAWK to run Bashlog
```
brew install mawk
```

3. Install miniconda (restart your terminal afterwards)
```
brew install --cask miniconda
conda init zsh # (or conda init bash)
```

4. Create a Python virtual environment
```
cd gelato
conda env create -f environment.yml -n gelato
conda activate gelato
```

5. Make `GELATO` script executable
```
chmod +x gelato
```

6. Check that `GELATO` is working
```
./gelato -j1
```

#### Configuration

Setup `DATABASE_GROUP` and its connection credentials.

  1. Only Neo4j and Blazegraph Databases need the credentials.

  2. If you haven't done so, create an empty file called `credentials.yaml` in your `GELATO` root directory: 

  3. Add the following lines to `credentials.yaml` and replace your database-specific credentials (indentation matters):

  - For Neo4j Database
  ``` yaml
  MGDB_NEO4J:
    protocol: MY_PROTOCOL
    host: MY_HOST
    password: MY_PASSWORD
    port: MY_PORT
    user: MY_USER
  ```
  Usually, we set `MY_PROTOCOL=neo4j`, `MY_HOST=localhost`, `MY_PORT=7687`, and `user=neo4j`.

  - For Blazegraph Database
  ```yaml
  MGDB_BLAZEGRAPH:
    host: MY_HOST
    port: MY_PORT
  ```
  Usually, we set `MY_HOST=localhost` and `MY_PORT=9999`.

#### DATABASES

Neo4j Database

> If you already have data stored in local Neo4j database, start the database before running `GELATO`. Otherwise, please set `{graph}[CYPHER][LOAD_DATA]=True` in `credentials.yaml` file and follow the steps below.

  1. Download and install [Neo4j Desktop](https://neo4j.com/download-neo4j-now/)

  2. Create an empty Database

  - Open Neo4j Desktop and click `Add`, then select `Local DBMS`

  - Fill in `Name` and `Password` fields. Please make sure `Password` field is consistent with `[MY_PASSWORD]` field of `credentials.yaml` file

  - Click `...` at the end of the database you just created and select `Settings`. Then update `dbms.memory.heap.max_size=10G` to get enough memory

  - Start the database you just created

  3. Put the relavant TSV files under the `import` folder

  - In Neo4j Desktop, click `...` at the end of the database you just created, select `Open folder` > `Import`

  - Copy the `data/raw/{graph}/tsv` folder and paste it under the `import` folder

Blazegraph Database

> If you already have data stored in local Blazegraph database, start the database before running `GELATO`. Otherwise, please set `{graph}[SPARQL][LOAD_DATA]=True` in `credentials.yaml` file and follow the steps below.

  1. Download the [blazegraph.jar](https://github.com/blazegraph/database/releases/latest)

  2. Start the Blazegraph server via terminal

  ```
  java -server -Xmx4g -jar blazegraph.jar
  ```

  3. Go to `http://localhost:9999/blazegraph/#namespaces` website and create a namespace named `{graph}` (e.g., MGDB) with `rdr` Mode.

# GELATO

*G*raph *E*xploration: *L*anguages *a*nd *To*ols

!!! tip "Code structure"

  ├── README.md
  ├── Snakefile
  ├── config.yaml
  ├── credentials.yaml
  ├── data
  │   ├── external
  │   │   └── mgdb
  │   │       ├── advised.csv
  │   │       ├── dissertation.tsv
  │   │       └── person.csv
  │   ├── query
  │   │   └── mgdb
  │   │       ├── python
  │   │       │   ├── benchmark_search_descendants_for_63244.txt
  │   │       │   └── output_search_descendants_for_63244.txt
  │   │       └── sql
  │   │           ├── benchmark_search_descendants_for_63244.txt
  │   │           └── output_search_descendants_for_63244.txt
  │   └── raw
  │       └── mgdb
  │           └── python
  │               ├── advised.tsv
  │               ├── dissertation.tsv
  │               └── person.tsv
  ├── environment.yml
  ├── gelato
  ├── rules
  │   ├── common.smk
  │   ├── preprocessing.smk
  │   └── query.smk
  ├── src
  │   ├── preprocess
  │   │   └── mgdb_preprocess_data_for_python.py
  │   └── query
  │       ├── mgdb_search_descendants_with_python.py
  │       └── mgdb_search_descendants_with_sql.py
  └── tools
      └── config.schema.yaml

### Setup 

#### Installation (For Mac)

1. Install miniconda (restart your terminal afterwards)
```
brew install --cask miniconda
conda init zsh # (or conda init bash)
```

2. Create a Python virtual environment
```
cd gelato
conda env create -f environment.yml -n gelato
conda activate gelato
```

3. Make `GELATO` script executable
```
chmod +x gelato
```

4. Check that `GELATO` is working
```
./gelato -j1
```

#### Configuration

Setup `DATABASE_GROUP` and its connection credentials.

  1. If you haven't done so, create an empty file called `#!bash credentials.yaml` in your GELATO root directory: 

  2. Add the following lines to `credentials.yaml` and replace your database-specific credentials (user, password, host, and database):

  ``` yaml
  MY_GROUP:
    database: MY_DATABASE
    host: MY_HOST
    password: MY_PASSWORD
    port: MY_PORT
    user: MY_USER
  ```

  3. Notes
  
    1. The label `[MY_GROUP]` is arbitrary but it has to match the `[DATABASE_GROUP]` attribute of the data stream you choose to use.

    2. Indentation matters

    3. You can have more than one credentials group in `credentials.yaml`

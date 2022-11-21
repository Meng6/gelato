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
  
  - The label `[MY_GROUP]` is arbitrary but it has to match the `[DATABASE_GROUP]` attribute of the data stream you choose to use.

  - Indentation matters

  - You can have more than one credentials group in `credentials.yaml`

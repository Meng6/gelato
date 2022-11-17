import glob

def get_mgdb_external_data(wildcards):
    for file_path in glob.glob(config["MGDB"]["DATA_SOURCE"]["FOLDER_PATH"] + "/*"):
        file_name = file_path.split("/")[-1].split(".")[0]
        if file_name == wildcards.graph_data:
            return file_path
    raise ValueError("MGDB: external file {graph_data} does not exists".format(wildcards.graph_data))


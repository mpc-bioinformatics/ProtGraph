def get_file_from_name(name: list | str):
    if type(name) == list:
        file_txt = name[0]
        file_graphml = file_txt.replace("txt", "graphml")
    else:
        file_graphml = name
    return file_graphml

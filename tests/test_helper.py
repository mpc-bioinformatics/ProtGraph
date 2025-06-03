def get_file_from_name(name_list: list):
    file_txt = name_list[0]
    file_graphml = file_txt.replace("txt", "graphml")
    return file_graphml

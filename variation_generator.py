

def _find_all_paths(graph, start, end):
    # Maybe in C or CPP ? Optimize this, it is not fast enough!
    # Generator for all possible paths
    path  = []
    queue = [(start, end, path)]
    while queue:
        start, end, path = queue.pop()
        path = path + [start]

        for node in set(graph.neighbors(start,mode='OUT')).difference(path):
            queue.append((node, end, path))

        if start == end and len(path) > 0:              
            yield path



def _get_path_count(graph, start, end):
    counter = 0
    for _ in _find_all_paths(graph, start, end):
        counter += 1
    return counter



def get_next_variant(graph_queue, prot_variant_queue):

    while True:
        graph_entry = graph_queue.get()
        [__start_node__] = graph_entry.vs.select(aminoacid="__start__")
        [__stop_node__]  = graph_entry.vs.select(aminoacid="__end__")

        # TODO we currently only count all possible paths
        num = _get_path_count(graph_entry,  __start_node__.index, __stop_node__.index)
        
        prot_variant_queue.put((graph_entry.vs[1]["accession"], num))
import igraph

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

    # TODO remove this!
    # k = []
    # supergraph = igraph.Graph(directed=True) 

    while True:
        try: 
            graph_entry = graph_queue.get(timeout=60)
        except Exception:
            continue

        
            
        ### TODO for generating a graph containing all proteins
        # TODO REMOVE THIS!!
        # supergraph = supergraph + graph_entry
        # k.append(graph_entry)

        # supergraph = igraph.Graph(directed=True) 
        # for entry in k:
        #     supergraph += entry
        # supergraph.write_dot("supergraph.dot")
        ###


        ### TODO Remove_
        # writing a dot file for each graph
        graph_entry.es[:]["qualifiers"] = [",".join([k.type for k in x["qualifiers"]]) if x["qualifiers"] is not None and len(x["qualifiers"]) != 0 else "" for x in list(graph_entry.es[:]) ]
        e = graph_entry.ecount()
        v = graph_entry.vcount()
        acc = graph_entry.vs[1]["accession"]
        graph_entry.write_dot("tempfolder2/" + acc + ".dot")


        [__start_node__] = graph_entry.vs.select(aminoacid="__start__")
        [__stop_node__]  = graph_entry.vs.select(aminoacid="__end__")

        # TODO we currently only count all possible paths
        # TODO exclude currently !!!! NOTE DL this takes too long!
        # num = _get_path_count(graph_entry,  __start_node__.index, __stop_node__.index)


        
        prot_variant_queue.put((graph_entry.vs[1]["accession"], 1))
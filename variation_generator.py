import igraph

def _find_all_paths(graph, start, end, cutoff):
    # Maybe in C or CPP ? Optimize this, it is not fast enough!
    # Generator for all possible paths
    path  = []
    queue = [(start, end, path)]
    while queue:
        start, end, path = queue.pop()
        if len(path) > cutoff: # TODO remove this MODIFICATION of the all path algorithm here!
            continue
        path = path + [start]

        for node in set(graph.neighbors(start,mode='OUT')).difference(path):
            queue.append((node, end, path))

        if start == end and len(path) > 0:              
            yield path

# [__start_node__] = graph_entry.vs.select(aminoacid="__start__")
# [__stop_node__]  = graph_entry.vs.select(aminoacid="__end__")

# TODO we currently only count all possible paths
# TODO exclude currently !!!! NOTE DL this takes too long!
# num = _get_path_count(graph_entry,  __start_node__.index, __stop_node__.index)

def _get_path_count(graph, start, end, cutoff):
    counter = 0
    for _ in _find_all_paths(graph, start, end, cutoff):
        counter += 1
    return counter







def _digest_via_trypsin(graph):

    [__start_node__] = graph.vs.select(aminoacid="__start__")
    [__end_node__] = graph.vs.select(aminoacid="__end__")

    k_s = graph.vs.select(aminoacid="K")
    r_s = graph.vs.select(aminoacid="R")

    k_s_edges = [ graph.es.select(_source=k) for k in k_s ]
    r_s_edges = [ graph.es.select(_source=r) for r in r_s ]


    # Get all edges which would be cleaved
    k_s_edges_new = [y for x in k_s_edges for y in x if graph.vs[y.target]["aminoacid"] != "P"]
    r_s_edges_new = [y for x in r_s_edges for y in x if graph.vs[y.target]["aminoacid"] != "P"]

    # Digest the graph into smaller parts
    # Explicitly we add 

    trypsin_in = [(__start_node__, k.target) for k in k_s_edges_new] # Nodes which should have an edge to start
    trypsin_out = [(k.source, __end_node__) for k in k_s_edges_new] # Node which should have an edge to end


    #  Mark digested edges as cleaved
    cleaved_idcs = [x.index for x in k_s_edges_new + r_s_edges_new]
    graph.es[cleaved_idcs]["cleaved"] = True


    # Add The newly edges for starting and ending
    graph.add_edges(trypsin_in)
    graph.add_edges(trypsin_out)











def get_next_variant(graph_queue, prot_variant_queue, cutoff=60):
    # cutoff, the maximal length of paths




    # TODO remove this!
    # k = []
    # supergraph = igraph.Graph(directed=True) 

    while True:
        try: 
            graph_entry = graph_queue.get(timeout=180)
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


        ## TODO Remove_
        # Use this if gephi is stuck on initializing graph
        # export LIBGL_ALWAYS_SOFTWARE=1

        # writing a dot file for each graph
        # graph_entry.es[:]["qualifiers"] = [",".join([k.type for k in x["qualifiers"]]) if x["qualifiers"] is not None and len(x["qualifiers"]) != 0 else "" for x in list(graph_entry.es[:]) ]
        # e = graph_entry.ecount()
        # v = graph_entry.vcount()
        # acc = graph_entry.vs[1]["accession"]
        # graph_entry.write_dot("protein_graph/" + acc + ".dot")


        # TODO before and after digest!!!?!?

        # Do a digestion via Trypsin (currently only this case)
        _digest_via_trypsin(graph_entry)

        # writing a dot file for each graph
        # graph_entry.es[:]["qualifiers"] = [",".join([k.type for k in x["qualifiers"]]) if x["qualifiers"] is not None and len(x["qualifiers"]) != 0 else "" for x in list(graph_entry.es[:]) ]
        # e = graph_entry.ecount()
        # v = graph_entry.vcount()
        # acc = graph_entry.vs[1]["accession"]
        # graph_entry.write_dot("cleaved_protein_graph/" + acc + ".dot")





        [__start_node__] = graph_entry.vs.select(aminoacid="__start__")
        [__end_node__] = graph_entry.vs.select(aminoacid="__end__")
        # all_paths = graph_entry.get_all_simple_paths(__start_node__, to=__end_node__, cutoff=cutoff)

        all_paths_count = _get_path_count(graph_entry, __start_node__, __end_node__, cutoff)



        # TODO we return the number of possible subgraphs (without the counting of possible variations!)
        prot_variant_queue.put(    (  graph_entry.vs[1]["accession"], all_paths_count )    )
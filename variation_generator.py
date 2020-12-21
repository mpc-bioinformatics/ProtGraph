import igraph

from digestion import digest

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


def _get_path_count(graph, start, end, cutoff):
    counter = 0
    for _ in _find_all_paths(graph, start, end, cutoff):
        counter += 1
    return counter






def _get_mass_dict():
    # So only integers are returned
    # CONVERTED BY FACTOR: 1 000 000 000
    factor = 1000000000
    return dict(
        G=(int(57.021463735*factor),  int(57.05132*factor) ),
        A=(int(71.037113805*factor),  int(71.0779*factor) ),
        S=(int(87.032028435*factor),  int(87.0773*factor) ),
        P=(int(97.052763875*factor),  int(97.11518*factor) ),
        V=(int(99.068413945*factor),  int(99.13106*factor) ),
        T=(int(101.047678505*factor), int(101.10388*factor) ),
        C=(int(103.009184505*factor), int(103.1429*factor) ),
        L=(int(113.084064015*factor), int(113.15764*factor) ),
        I=(int(113.084064015*factor), int(113.15764*factor) ),
        N=(int(114.042927470*factor), int(114.10264*factor) ),
        D=(int(115.026943065*factor), int(115.0874*factor) ),
        Q=(int(128.058577540*factor), int(128.12922*factor) ),
        K=(int(128.094963050*factor), int(128.17228*factor) ),
        E=(int(129.042593135*factor), int(129.11398*factor) ),
        M=(int(131.040484645*factor), int(131.19606*factor) ),
        H=(int(137.058911875*factor), int(137.13928*factor) ),
        F=(int(147.068413945*factor), int(147.17386*factor) ),
        U=(int(150.953633405*factor), int(150.3079*factor) ),
        R=(int(156.101111050*factor), int(156.18568*factor) ),
        Y=(int(163.063328575*factor), int(163.17326*factor) ),
        W=(int(186.079312980*factor), int(186.2099*factor) ),
        O=(int(237.147726925*factor), int(237.29816*factor) ),

        J=(int(113.084064015*factor), int(113.1594*factor) ),
        X=(int(0.0*factor), int(0.0*factor) ), # Unknown Amino Acid
        Z=(int(128.55059*factor), int(128.6231*factor) ),
        B=(int(114.53495*factor), int(114.5962*factor) ),

        __start__=(0, 0),
        __end__=(0, 0)
    )









def get_next_variant(graph_queue, cutoff=60):
    # cutoff, the maximal length of paths

    # TODO we digest here only!

    # get weight dict
    # d["AMINOACID"][0] -> Mono Mass
    # d["AMINOACID"][1] -> Avrg Mass
    mass_dict = _get_mass_dict()

    graph_entry = graph_queue
    # while True:
    #     try: 
    #         graph_entry = graph_queue.get(timeout=180)
    #     except Exception:
    #         continue

    ## TODO Remove_
    # Use this if gephi is stuck on initializing graph
    # export LIBGL_ALWAYS_SOFTWARE=1

    # writing a dot file for each graph
    # graph_entry.es[:]["qualifiers"] = [",".join([k.type for k in x["qualifiers"]]) if x["qualifiers"] is not None and len(x["qualifiers"]) != 0 else "" for x in list(graph_entry.es[:]) ]
    # e = graph_entry.ecount()
    # v = graph_entry.vcount()
    # acc = graph_entry.vs[1]["accession"]
    # graph_entry.write_dot("protein_graph/" + acc + ".dot")


    # Do a digestion via Trypsin (currently only this case)
    # TODO parameter is missing!
    # digest(graph_entry, "trypsin")
    # _digest_via_trypsin(graph_entry)



    #####################
    # Optimize and summarize chains of nodes

    # Get Chains by iterating over DAG
    chain = []
    for e in graph_entry.vs[:]:
        if e.outdegree() == 1 and e.indegree() == 1: # TODO correct?
            prev_node = e.neighbors(mode="IN")[0]
            if prev_node.outdegree() == 1:
                chain.append( (prev_node.index, e.index) )

    # Chain those elements together
    complete_chain = []
    for ele in chain:
        idx = None
        for cc_idx, cc in enumerate(complete_chain):
            if cc[-1] == ele[0]:
                idx = cc_idx
        
        if idx == None:
            complete_chain.append([ele[0], ele[1]])
        else:
            complete_chain[idx].append(ele[1])



    # Generate intermediate Graph Mapping 

    # Special case exclude start and ending:
    [__start_node__] = graph_entry.vs.select(aminoacid="__start__")
    [__end_node__] = graph_entry.vs.select(aminoacid="__end__")
    complete_chain = [[y for y in x if y != __start_node__.index and y != __end_node__.index] for x in complete_chain]

    # Contains [(IDCS, VERTICES, OUT_EDGES)] all attributes in tuple as list
    merge_list = [] 
    for c in complete_chain:
        vs = []
        es = []
        first_node = None
        for x in c:
            vs.append(graph_entry.vs[x].attributes())
            es.append(graph_entry.vs[x].out_edges()[0])
            outs = graph_entry.vs[x].neighbors(mode="IN")
            if outs[0].index not in c:
                first_node = x


        merge_list.append((first_node, c + [], vs, es))



    # Sort them accordingly forom beginning and concat the attributes
    merged_nodes = []
    for first_node, idcs, vertices, edges in merge_list:

        # Create sorted list of edges and nodes (since we may still have some inconsistencies?!? TODO check!? )
        sorted_edges = []
        sorted_nodes = []
        while first_node in idcs:
            fn_idx = idcs.index(first_node)
            sorted_edges.append(edges[fn_idx])
            sorted_nodes.append(vertices[fn_idx])
            first_node = sorted_edges[-1].target
        ### Special case retrieve information of last node also 
        last_node = sorted_edges[-1].source
        if sorted_edges[-1].target_vertex.indegree() == 1:
            if sorted_edges[-1].target_vertex["aminoacid"] != "__end__":
                last_node = sorted_edges[-1].target
                sorted_nodes.append(sorted_edges[-1].target_vertex.attributes())
        
        # TODO check if only one element:
        iso_set = set([x["isoform_accession"] for x in sorted_nodes if x["isoform_accession"] is not None])
        acc_set = set([x["accession"] for x in sorted_nodes  if x["accession"] is not None])

        if len(acc_set) > 1:
            print("DEBUG ME")
        if len(iso_set) > 1:
            print("DEBUG ME") 

        #  Merged Node atributes # TODO
        m_accession = acc_set.pop() if len(acc_set) == 1 else None # Should only Contain 1 Element
        m_isoform_accession = iso_set.pop() if len(iso_set) == 1 else None # Should only Contain 1 Element
        m_position = sorted_nodes[0]["position"]
        m_isoform_position = sorted_nodes[0]["isoform_position"]
        m_aminoacid = "".join([x["aminoacid"] for x in sorted_nodes])

        #  Merges Edges Attributes # TODO

        cle_set = set([x.attributes()["cleaved"] for x in sorted_edges if x.attributes()["cleaved"] is not None])
        if len(acc_set) > 1:
            print("DEBUG ME")

        m_cleaved = cle_set.pop() if len(cle_set) == 1 else None # Should be None # TODO
        m_qualifiers = [y for x in sorted_edges if x.attributes()["qualifiers"] is not None for y in x.attributes()["qualifiers"]]


        # Generate new Node/Edge attrs
        new_node_attrs = dict(
            accession=m_accession, isoform_accession=m_isoform_accession, position=m_position,
            isoform_position=m_isoform_position,  aminoacid=m_aminoacid
        )
        new_edge_attrs = dict(
            cleaved=m_cleaved,
            qualifiers=m_qualifiers
        )

        idcs_to_remove = idcs
        if last_node not in idcs_to_remove:
            idcs_to_remove.append(last_node)
        # save merged information back to list
        merged_nodes.append([
            idcs_to_remove, # Nodes idcs which needs to be removed
            [x.index for x in sorted_edges], # Edges which needs to be removed
            sorted_edges[0].source, # first Node idx where the edges into it needs to be modified
            last_node, # last Node idx where the edge outgoing needs to be modified
            new_node_attrs, # the attributes of the supernode
            new_edge_attrs, # the attributes of the edge ingoing from the supernode (appending at the end!)

        ])


    
    ### Modifiy the graph by adding the new nodes
    ### adding saved edges to the supernode
    ###
    ### removing the edges first in supernode
    ### removing the edges at beginning of supernode (saving them for later)
    ### removing the edges at end of supernode (saving them for later)
    ### removing the nodes in supernode
    
    
    # 1. add all supernodes and its attributes!
    cur_node_count = graph_entry.vcount()
    graph_entry.add_vertices(len(merged_nodes))
    graph_entry.vs[cur_node_count:]["accession"] =  [x[4]["accession"] for x in merged_nodes]
    graph_entry.vs[cur_node_count:]["isoform_accession"] =  [x[4]["isoform_accession"] for x in merged_nodes]
    graph_entry.vs[cur_node_count:]["position"] =  [x[4]["position"] for x in merged_nodes]
    graph_entry.vs[cur_node_count:]["isoform_position"] =  [x[4]["isoform_position"] for x in merged_nodes]
    graph_entry.vs[cur_node_count:]["aminoacid"] =  [x[4]["aminoacid"] for x in merged_nodes]

    supernode_idcs = list(range(cur_node_count, cur_node_count + len(merged_nodes)))


    # 2. Get all edges connected to supernode and them appropiately!
    for supernode_idx, (_, _, first_n, last_n, _, edge_attr) in zip(supernode_idcs, merged_nodes):
        # Get all in edges
        in_e_new_attrs = []
        for in_e in graph_entry.vs[first_n].in_edges():
            attrs = in_e.attributes()
            if attrs["qualifiers"] == None:
                attrs["qualifiers"] = edge_attr["qualifiers"]
            else:
                attrs["qualifiers"] = attrs["qualifiers"] + edge_attr["qualifiers"]
            in_e_new_attrs.append( ((in_e.source, supernode_idx), attrs ) )

        # Get all out edges
        out_e_new_attrs = []
        for out_e in graph_entry.vs[last_n].out_edges():
            attrs = out_e.attributes()
            out_e_new_attrs.append( ((supernode_idx, out_e.target), attrs)  )

        # Add them to the graph
        e_count = graph_entry.ecount()
        graph_entry.add_edges([x[0] for x in in_e_new_attrs + out_e_new_attrs])
        graph_entry.es[e_count:]["qualifiers"] = [x[1]["qualifiers"] for x in in_e_new_attrs + out_e_new_attrs]


    # 3. remove all edges (then nodes)        
    all_nodes_to_remove = [y for x in merged_nodes for y in x[0]]
    all_edges_to_remove = [y for x in merged_nodes for y in x[1]]
    graph_entry.delete_edges(all_edges_to_remove)
    graph_entry.delete_vertices(all_nodes_to_remove)



    
    #####################
    # TODO check DAG in and outdegree specifically!!!

    if not graph_entry.is_dag():
        print("DEBUG ME")
    if graph_entry.indegree().count(0) != 1:
        print("DEBUG ME")
    if graph_entry.outdegree().count(0) != 1:
        print("DEBUG ME")

    #####################

    # mono_masses = [mass_dict[x.target_vertex["aminoacid"]][0] for x in graph_entry.es[:]]
    # avrg_masses = [mass_dict[x.target_vertex["aminoacid"]][1] for x in graph_entry.es[:]]
    mono_masses = [sum([mass_dict[y][0] for y in x.target_vertex["aminoacid"].replace("__start__", "").replace("__end__", "") ]) for x in graph_entry.es[:]]
    avrg_masses = [sum([mass_dict[y][1] for y in x.target_vertex["aminoacid"].replace("__start__", "").replace("__end__", "") ]) for x in graph_entry.es[:]]


    


    graph_entry.es[:]["mono_weight"] = mono_masses
    graph_entry.es[:]["avrg_weight"] = avrg_masses


    graph_entry.vs[:]["mono_end_weight"] = float("inf")
    graph_entry.vs[:]["avrg_end_weight"] = float("inf")
    graph_entry.vs.select(aminoacid="__end__")["mono_end_weight"] = 0
    graph_entry.vs.select(aminoacid="__end__")["avrg_end_weight"] = 0

    # Get distance to end ("all pair shrotest path manner")
    sorted_nodes = graph_entry.topological_sorting(mode="IN")

    # Doing it for mono mass
    for node in sorted_nodes:
        for edge in graph_entry.es.select(_target=node):
            temp_weight = edge.target_vertex["mono_end_weight"] + edge["mono_weight"]
            if graph_entry.vs[edge.source]["mono_end_weight"] > temp_weight:
                graph_entry.vs[edge.source]["mono_end_weight"] = temp_weight

    # Doing it for avrg mass
    for node in sorted_nodes:
        for edge in graph_entry.es.select(_target=node):
            temp_weight = edge.target_vertex["avrg_end_weight"] + edge["avrg_weight"]
            if graph_entry.vs[edge.source]["avrg_end_weight"] > temp_weight:
                graph_entry.vs[edge.source]["avrg_end_weight"] = temp_weight




    #### drag the end_weights to the edges
    mono_end_masses = [x.target_vertex["mono_end_weight"] for x in graph_entry.es[:]]
    avrg_end_masses = [x.target_vertex["avrg_end_weight"] for x in graph_entry.es[:]]
    graph_entry.es[:]["mono_end_weight"] = mono_end_masses
    graph_entry.es[:]["avrg_end_weight"] = avrg_end_masses




    #########################


    ### Get the Number of all possible simple Paths!

    # First get top sort
    sorted_nodes = graph_entry.topological_sorting()

    # Get Number of Variable Paths (Init to 0)
    var_paths = [0]*graph_entry.vcount()
    first = sorted_nodes[0]
    var_paths[first] = 1 # Path to itself is zero! For convenience we set it to one 
    for v in graph_entry.neighbors(first, mode="OUT"):
        var_paths[v] = 1
    
    # Iterative approach look how many paths are possible from preve to itself
    for v in sorted_nodes[1:]:
        summed = 0
        for v_prev in graph_entry.neighbors(v, mode="IN"):
            summed += var_paths[v_prev]

        var_paths[v] = summed

    var_paths[first] = 0 # Path to itself is zero!
    num_of_paths = var_paths[sorted_nodes[-1]] # This contains the number of Paths







    ######################################
    # writing a dot file for each graph
    # graph_entry.es[:]["qualifiers"] = [",".join([k.type for k in x["qualifiers"]]) if x["qualifiers"] is not None and len(x["qualifiers"]) != 0 else "" for x in list(graph_entry.es[:]) ]
    # e = graph_entry.ecount()
    # v = graph_entry.vcount()
    # acc = graph_entry.vs[1]["accession"]
    # graph_entry.write_dot("cleaved_protein_graph/" + acc + ".dot")





    # [__start_node__] = graph_entry.vs.select(aminoacid="__start__")
    # [__end_node__] = graph_entry.vs.select(aminoacid="__end__")
    # all_paths = graph_entry.get_all_simple_paths(__start_node__, to=__end_node__, cutoff=cutoff)

    # all_paths_count = _get_path_count(graph_entry, __start_node__, __end_node__, cutoff)



    # TODO we return the number of possible subgraphs (without the counting of possible variations!)
    # prot_variant_queue.put(    (  graph_entry.vs[1]["accession"], None )    )
    # prot_variant_queue.put(  (graph_entry, num_of_paths)   )
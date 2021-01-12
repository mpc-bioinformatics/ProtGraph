
def merge_aminoacids(graph_entry):
    """ This merging reduces the amount of nodes and edges drastically.

        Each node in the graph represents a single aminoacid when built. It is very
        likely that each graph has many chains of aminoacids without interfering edges 
        inbetween them. 

        In essence this method does the following:
        E.G: 
             M -v                   ,-> K
                A -> M -> I -> N -> O 
             K -^                   `-> P
            (9 Nodes, 8 Edges)

        can be summarized to:
             M -v   ,-> K
                AMINO 
             K -^   `-> P
            (5 Nodes, 4 Edges)


        NOTE: If there are Proteins which can be completely summarized into a single node,
        then those are split into 3 Nodes and 2 Edges, so that the __start__ and __end__ node
        are present in all graphs. E.G:   __start__ -> PROTEIN -> __end__
    """
    #####################
    # Optimize and summarize chains of nodes

    # Get Chains by iterating over DAG
    chain = []
    for e in graph_entry.vs[:]:
        if e.outdegree() == 1 and e.indegree() == 1: # TODO correct?
            prev_node = e.neighbors(mode="IN")[0]
            if prev_node.outdegree() == 1:
                chain.append((prev_node.index, e.index))

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
    complete_chain = [
        [y for y in x if y != __start_node__.index and y != __end_node__.index]
        for x in complete_chain
    ]

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
            in_e_new_attrs.append(((in_e.source, supernode_idx), attrs))

        # Get all out edges
        out_e_new_attrs = []
        for out_e in graph_entry.vs[last_n].out_edges():
            attrs = out_e.attributes()
            out_e_new_attrs.append(((supernode_idx, out_e.target), attrs))

        # Add them to the graph
        e_count = graph_entry.ecount()
        graph_entry.add_edges([x[0] for x in in_e_new_attrs + out_e_new_attrs])
        graph_entry.es[e_count:]["qualifiers"] = [x[1]["qualifiers"] for x in in_e_new_attrs + out_e_new_attrs]


    # 3. remove all edges (then nodes)        
    all_nodes_to_remove = [y for x in merged_nodes for y in x[0]]
    all_edges_to_remove = [y for x in merged_nodes for y in x[1]]
    graph_entry.delete_edges(all_edges_to_remove)
    graph_entry.delete_vertices(all_nodes_to_remove)



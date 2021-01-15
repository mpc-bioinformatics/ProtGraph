def merge_aminoacids(graph_entry):
    """
    This merging reduces the amount of nodes and edges drastically.

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

    TODO DL this is very slow (takes ~23 minutes for complete human (without FTs, with trypsin))
    (~8 minutes with digestion skip)
    """
    # TODO Optimize and summarize chains of nodes, this is  slow (somewhere)!

    # Get single chain elements by iterating over the DAG
    # Here we retrieve those with in and out degree 1 and its prev node
    chain_elements = []
    for node in graph_entry.vs[:]:
        if node.outdegree() == 1 and node.indegree() == 1:  # TODO correct?
            prev_node = node.neighbors(mode="IN")[0]
            if prev_node.outdegree() == 1:
                chain_elements.append((prev_node.index, node.index))

    # Chain those elements together (get chains)
    complete_chain = []
    for ele in chain_elements:
        idx = None
        # For each tuple check if a node end
        # matches to a node beginning
        for cc_idx, cc in enumerate(complete_chain):
            if cc[-1] == ele[0]:
                idx = cc_idx

        if idx is None:
            # No matching, therefore we add both nodes
            complete_chain.append([ele[0], ele[1]])
        else:
            # Matching, we append the corresponding end node
            complete_chain[idx].append(ele[1])

    # Special case exclude start and stop Nodes:
    [__start_node__] = graph_entry.vs.select(aminoacid="__start__")
    [__end_node__] = graph_entry.vs.select(aminoacid="__end__")
    complete_chain = [
        [
            y  # Exclude start and end nodes here
            for y in x
            if y != __start_node__.index and y != __end_node__.index
        ]
        for x in complete_chain  # For each complete chain
    ]

    # Now we generate the intermediate graph mapping
    # This contains [(FIRST_NODE, IDCS, VERTICES, OUT_EDGES)]. All attributes in this tuple are lists
    merge_list = []
    for c in complete_chain:
        # For each chain
        vs = []
        es = []
        first_node = None
        for x in c:
            # For each element in chain
            # We add the vertex, its out going edeges
            vs.append(graph_entry.vs[x].attributes())
            es.append(graph_entry.vs[x].out_edges()[0])
            # And the first node as well as the indices (available through c)
            outs = graph_entry.vs[x].neighbors(mode="IN")
            if outs[0].index not in c:
                first_node = x
        merge_list.append((first_node, c + [], vs, es))

    # Sort them accordingly from beginning and concat the attributes
    merged_nodes = []
    for first_node, idcs, vertices, edges in merge_list:
        # For each complete chain -> merged information list entry

        # We create sorted lists of edges and nodes (since we may still have some inconsistencies?!? TODO check!?)
        # We reiterate the chains here again (for consistency!)
        sorted_edges = []
        sorted_nodes = []
        while first_node in idcs:
            fn_idx = idcs.index(first_node)
            sorted_edges.append(edges[fn_idx])
            sorted_nodes.append(vertices[fn_idx])
            first_node = sorted_edges[-1].target
        # Special case: Retrieve information of last node, sometimes it may be the stop node
        last_node = sorted_edges[-1].source
        if sorted_edges[-1].target_vertex.indegree() == 1:
            if sorted_edges[-1].target_vertex["aminoacid"] != "__end__":
                last_node = sorted_edges[-1].target
                sorted_nodes.append(sorted_edges[-1].target_vertex.attributes())

        # Then we merge node attributes
        # First, the standard attributes, which are always available
        m_aminoacid = "".join([x["aminoacid"] for x in sorted_nodes])  # Concat aminoacids
        m_position = sorted_nodes[0]["position"]  # Retrieve only the first position
        m_accession = _get_single_set_element(sorted_nodes, "accession")  # Get the ONLY accession
        # Now the attributes in nodes, which may be present
        m_isoform_accession = _get_single_set_element(sorted_nodes, "isoform_accession")  # Get the ONLY iso_accession
        m_isoform_position = [sorted_nodes[0]["isoform_position"] if "isoform_position" in sorted_nodes[0] else None]

        # Merge edges attributes, which also may be present!
        sorted_nodes_attrs = [x.attributes() for x in sorted_edges]
        m_cleaved = _get_single_set_element(sorted_nodes_attrs, "cleaved")  # cleaved consistency check. It is unused
        if "qualifiers" in sorted_nodes_attrs[0]:
            m_qualifiers = [
                y  # Special merging for qualifiers, we concat the (non-None) lists here
                for x in sorted_edges
                if x.attributes()["qualifiers"] is not None
                for y in x.attributes()["qualifiers"]
            ]
        else:
            m_qualifiers = None

        # Generate new node/edge dict attrs
        # Here we set all available! (So this may need to be appended for new attrs)
        new_node_attrs = dict(
            accession=m_accession, isoform_accession=m_isoform_accession, position=m_position,
            isoform_position=m_isoform_position, aminoacid=m_aminoacid,
        )
        new_edge_attrs = dict(cleaved=m_cleaved, qualifiers=m_qualifiers)

        # Retrieve the indices (from original graph), whihc need to be removed
        idcs_to_remove = idcs
        if last_node not in idcs_to_remove:
            idcs_to_remove.append(last_node)
        # Save merged information back to list
        merged_nodes.append(
            [
                idcs_to_remove,  # Nodes idcs which needs to be removed
                [x.index for x in sorted_edges],  # Edges which needs to be removed
                sorted_edges[0].source,  # First Node idx where the edges into it needs to be modified
                last_node,  # Last Node idx where the edge outgoing needs to be modified
                new_node_attrs,  # The attributes of the supernode
                new_edge_attrs,  # The attributes of the edges ingoing from the supernode (appending at the end!)
            ]
        )

    # Modifiy the graph by adding the new nodes in the following order
    # 1. adding saved supernodes and their attributes
    # 2. Get all edges and add them with attributes
    # 3. removing all edges first then nodes

    # 1. add all supernodes and its attributes!
    cur_node_count = graph_entry.vcount()
    graph_entry.add_vertices(len(merged_nodes))
    # First add the always available attributes
    _add_node_attributes(graph_entry, merged_nodes, cur_node_count, "aminoacid")
    _add_node_attributes(graph_entry, merged_nodes, cur_node_count, "position")
    _add_node_attributes(graph_entry, merged_nodes, cur_node_count, "accession")
    # Then add the possibly available attributes
    _add_node_attributes(graph_entry, merged_nodes, cur_node_count, "isoform_accession")
    _add_node_attributes(graph_entry, merged_nodes, cur_node_count, "isoform_position")

    # 2. Get all edges connected to supernode and add them appropiately!
    supernode_idcs = list(range(cur_node_count, cur_node_count + len(merged_nodes)))
    for supernode_idx, (_, _, first_n, last_n, _, edge_attr) in zip(supernode_idcs, merged_nodes):
        # Get all in edges
        in_e_new_attrs = []
        for in_e in graph_entry.vs[first_n].in_edges():
            attrs = in_e.attributes()
            new_attrs = _get_edge_attrs(attrs, edge_attr["qualifiers"])
            in_e_new_attrs.append(((in_e.source, supernode_idx), new_attrs))

        # Get all out edges
        out_e_new_attrs = []
        for out_e in graph_entry.vs[last_n].out_edges():
            attrs = out_e.attributes()
            out_e_new_attrs.append(((supernode_idx, out_e.target), attrs))

        # Add edges to the graph and add all possibly available attributes
        e_count = graph_entry.ecount()
        graph_entry.add_edges([x[0] for x in in_e_new_attrs + out_e_new_attrs])
        _add_edge_attributes(graph_entry, in_e_new_attrs + out_e_new_attrs, e_count, "qualifiers")
        _add_edge_attributes(graph_entry, in_e_new_attrs + out_e_new_attrs, e_count, "cleaved")

    # 3. remove removable edges (then removable nodes) at once
    all_nodes_to_remove = [y for x in merged_nodes for y in x[0]]
    all_edges_to_remove = [y for x in merged_nodes for y in x[1]]
    graph_entry.delete_edges(all_edges_to_remove)
    graph_entry.delete_vertices(all_nodes_to_remove)


def _get_edge_attrs(edge_attrs, concat_qualifiers):
    """ get edge attrs, returns for qualifiers always a list  """
    attrs = dict()
    if "qualifiers" not in edge_attrs:
        attrs["qualifiers"] = []
    elif edge_attrs["qualifiers"] is None:
        attrs["qualifiers"] = edge_attrs["qualifiers"]
        if attrs["qualifiers"] is None:
            attrs["qualifiers"] = []
    else:
        attrs["qualifiers"] = edge_attrs["qualifiers"]

    # add here the qualifiers afterwards from merged supernodes
    if concat_qualifiers is not None:
        attrs["qualifiers"] += concat_qualifiers

    for key in edge_attrs.keys():
        if key != "qualifiers":
            attrs[key] = edge_attrs[key]

    return attrs


def _add_edge_attributes(graph_entry, merged_edge_list, edge_start_count, key):
    """ We check here, if from the original graph a key is already set and set it then here appropiately"""
    if key not in graph_entry.es[0].attributes():
        # We skip adding the "new" attribute in the graph
        return

    # Add attributes
    graph_entry.es[edge_start_count:][key] = [x[1][key] for x in merged_edge_list]


def _add_node_attributes(graph_entry, merged_nodes_list, node_start_count, key):
    """ We check here, if from the original graph a key is already set and set it then here appropiately"""
    if key not in graph_entry.vs[0].attributes():
        # We skip adding the "new" attribute in the graph
        return

    # Add attributes
    graph_entry.vs[node_start_count:][key] = [x[4][key] for x in merged_nodes_list]


def _get_single_set_element(node_edges, key):
    """ Only return element if key exists and is unique! TODO """

    if key not in node_edges[0]:
        # If key not in nodes than we return simply None
        return None

    # Get list of keys, skipping None values
    attrs_list = [x[key] for x in node_edges]  # TODO do we skip None Values here?
    # attrs_list = [x[key] for x in node_edges if x is not None]
    attrs_set = set(attrs_list)

    if len(attrs_set) > 1:
        # TODO custom exception, this should not happen!
        # track information like key, which node/edge especially
        # retrieve the accession explicitly here too for debugging!
        raise Exception("Multiple entries in set!")

    return attrs_set.pop()

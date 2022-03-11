def traverse_to_end(graph_entry, complete_chain, single_nodes, single_in, next_node, c):
    # iterate as long as possible.
    while next_node in single_nodes:
        single_nodes.remove(next_node)
        c.append(next_node)
        next_node = graph_entry.vs[next_node].neighbors(mode="OUT")[0].index

    # Special case for single in
    if next_node in single_in:
        single_in.remove(next_node)
        c.append(next_node)
        next_node = graph_entry.vs[next_node].neighbors(mode="OUT")[0].index

    # Skip chains containing only 1 element
    if len(c) != 1:
        complete_chain.append(c)


def find_chains(graph_entry):
    """
    Retrive chain of nodes.
    Here we essentially extract all chains and return their vertices indices
    E.G Graph:
         M -v                   ,-> K
            A -> M -> I -> N -> O
         K -^                   `-> P
        (9 Nodes, 8 Edges)

    would give us the node indices of [A M I N O] in a list (and all other available
    chains in the graph)
    """
    [__start_node__] = graph_entry.vs.select(aminoacid="__start__")

    # Sort all nodes into 3 possible bins
    single_nodes = set()
    single_out = set()
    single_in = set()
    for idx, (a, b) in enumerate(zip(graph_entry.vs.indegree(), graph_entry.vs.outdegree())):
        if a == 1 and b == 1:
            single_nodes.add(idx)  # case single
        if a == 1 and b > 1:
            single_in.add(idx)  # case one in but multiple out
        if a > 1 and b == 1:
            single_out.add(idx)  # case one out but multiple in

    # save all chains in this list
    complete_chain = []

    # CASE 1: Chain is starting with a single out node
    for so in single_out:
        c = [so]
        next_node = graph_entry.vs[so].neighbors(mode="OUT")[0].index

        traverse_to_end(graph_entry, complete_chain, single_nodes, single_in, next_node, c)

    # CASE 2: it may happen that the start node is at a beginning of chain.
    # Here we do a intersection of remaining nodes in single_nodes with the
    # nodes having a direct connection to start
    start_set = {x.index for x in __start_node__.neighbors(mode="OUT")}
    single_start_points = single_nodes.intersection(start_set)
    for sn in single_start_points:
        c = [sn]
        next_node = graph_entry.vs[sn].neighbors(mode="OUT")[0].index
        single_nodes.remove(sn)

        traverse_to_end(graph_entry, complete_chain, single_nodes, single_in, next_node, c)

    # return complete chain
    return complete_chain


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

    TODO DL This has been reimplemented (since 0.2.1) to be a bit faster (~20 times faster).
    but can probably be further optimized.
    """
    # Retrive all chains of nodes (sorted by chain order in list of lists)
    # TODO this is still the slowest part. It may be improved
    complete_chain = find_chains(graph_entry)

    # Generate merged nodes information (iow supernodes attributes)
    merged_nodes = []
    for idcs in complete_chain:
        # Retrieve nodes and edges information
        sorted_nodes = [graph_entry.vs[x].attributes() for x in idcs]
        sorted_edges = [graph_entry.vs[x].out_edges()[0] for x in idcs[:-1]]

        # Then we merge node attributes
        # First, the standard attributes, which are always available
        m_aminoacid = "".join([x["aminoacid"] for x in sorted_nodes])  # Concat aminoacids
        m_position = sorted_nodes[0]["position"]  # Retrieve only the first position
        m_accession = _get_single_set_element(sorted_nodes, "accession")  # Get the ONLY accession
        # Now the attributes in nodes, which may be present
        m_isoform_accession = _get_single_set_element(sorted_nodes, "isoform_accession")  # Get the ONLY iso_accession
        m_isoform_position = sorted_nodes[0]["isoform_position"] if "isoform_position" in sorted_nodes[0] else None
        m_delta_mass = sum(x["delta_mass"] if x["delta_mass"] else 0 for x in sorted_nodes) \
            if "delta_mass" in sorted_nodes[0] else None

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
            isoform_position=m_isoform_position, aminoacid=m_aminoacid, delta_mass=m_delta_mass
        )
        new_edge_attrs = dict(cleaved=m_cleaved, qualifiers=m_qualifiers)

        # Save merged information back to list
        merged_nodes.append(
            [
                idcs,  # Nodes idcs which needs to be removed
                [x.index for x in sorted_edges],  # Edges which need to be removed
                idcs[0],  # First Node idx where the edges into it needs to be modified
                idcs[-1],  # Last Node idx where the edges outgoing needs to be modified
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
    _add_node_attributes(graph_entry, merged_nodes, cur_node_count, "delta_mass")

    # 2. Get all edges connected to supernode and add them appropiately!
    supernode_idcs = list(range(cur_node_count, cur_node_count + len(merged_nodes)))
    new_edges_with_attrs = []  # List of all edges to be added
    # creating mapping for bulk edge inserts
    first_ns_dict = {
        first_n: supernode_idx
        for supernode_idx, (_, _, first_n, _, _, edge_attr) in zip(supernode_idcs, merged_nodes)
    }
    last_ns_dict = {
        last_n: supernode_idx
        for supernode_idx, (_, _, _, last_n, _, edge_attr) in zip(supernode_idcs, merged_nodes)
    }
    registered_edges = set()  # Set saving already added edges
    for supernode_idx, (_, _, first_n, last_n, _, edge_attr) in zip(supernode_idcs, merged_nodes):
        # For each supernode
        # Get all in edges
        for in_e in graph_entry.vs[first_n].in_edges():
            attrs = in_e.attributes()
            new_attrs = _get_edge_attrs(attrs, edge_attr["qualifiers"])
            # And add edge iff not already (with mapping to supernode)
            if in_e.source in last_ns_dict:
                if (last_ns_dict[in_e.source], supernode_idx) not in registered_edges:
                    new_edges_with_attrs.append(((last_ns_dict[in_e.source], supernode_idx), new_attrs))
                    registered_edges.add((last_ns_dict[in_e.source], supernode_idx))
            else:
                new_edges_with_attrs.append(((in_e.source, supernode_idx), new_attrs))

        # Get all out edges
        for out_e in graph_entry.vs[last_n].out_edges():
            attrs = out_e.attributes()
            # And add edge iff not already (with mapping to supernode)
            if out_e.target in first_ns_dict:
                if (supernode_idx, first_ns_dict[out_e.target]) not in registered_edges:
                    new_edges_with_attrs.append(((supernode_idx, first_ns_dict[out_e.target]), attrs))
                    registered_edges.add((supernode_idx, first_ns_dict[out_e.target]))
            else:
                new_edges_with_attrs.append(((supernode_idx, out_e.target), attrs))

    # Bulk add edges to the graph and add all possibly available attributes
    e_count = graph_entry.ecount()
    graph_entry.add_edges([x[0] for x in new_edges_with_attrs])
    _add_edge_attributes(graph_entry, new_edges_with_attrs, e_count, "qualifiers")
    _add_edge_attributes(graph_entry, new_edges_with_attrs, e_count, "cleaved")

    #####################################
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

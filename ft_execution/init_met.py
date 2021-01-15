def execute_init_met(graph, init_met_feature):
    """
    This function adds ONLY edges to skip the initiator metheonine.

    NOTE: This transforms the graph without returning it!

    Following Keys are set here:
    Nodes: <None>
    Edges: "qualifiers" ( -> adds INIT_MET)
    """
    # Get start node
    [start] = graph.vs.select(aminoacid="__start__")

    # Get all neighbors of start (next nodes)
    pos_first_aas = graph.neighbors(start, mode="out")

    # Filtering is important, since we may skip an aminoacid of an already cleaved protein
    met_aas = [
        x
        for x in pos_first_aas
        if graph.vs[x]["aminoacid"] == "M" and graph.vs[x]["position"] == 1
    ]

    # Get possible features of remaining nodes
    features = [_get_qualifiers(graph, start, x) for x in met_aas]

    # Get the next nodes which should be after them
    targets = [graph.neighbors(x, mode="out") for x in met_aas]

    # Get current number of edges
    cur_count = graph.ecount()

    # Generate the edges list
    edge_list = [(start, y) for x in targets for y in x]
    # Generate the corresponding qualifiers
    qualifiers_list = [
        [init_met_feature, *features[x_idx][y_idx]]
        for x_idx, x in enumerate(targets)
        for y_idx, _ in enumerate(x)
    ]

    # Add new edges (bulk)
    graph.add_edges(edge_list)
    graph.es[cur_count:]["qualifiers"] = qualifiers_list


def _get_qualifiers(graph, source_node: int, target_node: int):
    """ Retrieve the qualifiers list (empty list if nonexistent) """
    if "qualifiers" in graph.es[0].attributes():
        qualifiers_list = graph.es.select(_source=source_node, _target=target_node)["qualifiers"]
        qualifiers_list = [x if x is not None else [] for x in qualifiers_list]
        return qualifiers_list
    else:
        return [[]]

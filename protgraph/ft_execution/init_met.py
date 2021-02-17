from protgraph.unexpected_exception import UnexpectedException


def execute_init_met(graph, init_met_feature):
    """
    This function adds ONLY edges to skip the initiator methionine.

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
    # Depending on the set reference we need to distinguish between multiple INIT_METs for one protein
    met_aas = _get_methionines(graph, pos_first_aas, init_met_feature)

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


def _get_methionines(graph, pos_first_aas, init_met_feature):
    """ Wrapper function to retrieve the methionines, depending on reference and isoforms """
    # Check if graph has isoforms
    has_isoforms = True if "isoform_accession" in graph.vs[0].attributes() and \
        "isoform_position" in graph.vs[0].attributes() else False

    _aas_set = True  # Bool to check if we were able to retrieve aminoacids from the graph
    if init_met_feature.ref is None:
        # Canonical M
        if not has_isoforms:
            # for no isoforms
            met_aas = [
                x
                for x in pos_first_aas
                if graph.vs[x]["aminoacid"] == "M" and graph.vs[x]["position"] == 1
            ]
        else:
            # for isoforms
            met_aas = [
                x
                for x in pos_first_aas
                if graph.vs[x]["aminoacid"] == "M" and graph.vs[x]["position"] == 1 and
                graph.vs[x]["isoform_accession"] is None
            ]
    elif has_isoforms:
        # Referenced isoform M depending on ref
        met_aas = [
            x
            for x in pos_first_aas
            if graph.vs[x]["aminoacid"] == "M" and graph.vs[x]["isoform_position"] == 1 and
            graph.vs[x]["isoform_accession"] == init_met_feature.ref
        ]
    else:
        if init_met_feature.ref is not None and not has_isoforms:
            print(
                "Warning, INIT_MET could not applied on isoform {} (isoform is missing in graph)"
                .format(init_met_feature.ref)
            )
            met_aas = []
            _aas_set = False

        else:
            # Unknown case
            raise UnexpectedException(
                accession=graph.vs[0]["accession"],
                position=1,
                message="Unknown Case for the feature INIT_MET. No Aminoacid found to be skipped",
                additional_info=str(init_met_feature)
            )

    # Check if we found a M to skip, if not print a Warning, since this may be the case for
    # unreviewed entries...
    if _aas_set and len(met_aas) == 0:
        print(
            "WARNING: Protein '{}' does not have Methionine at the beginning, while "
            "trying to skip it via the feature 'INIT_MET'. Skipping ...".format(graph.vs[0]["accession"])
        )

    return met_aas

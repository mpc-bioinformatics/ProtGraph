from Bio.SeqFeature import UnknownPosition

from protgraph.ft_execution import _get_qualifiers
from protgraph.unexpected_exception import UnexpectedException


def execute_signal(graph, signal_feature):
    """
    This function adds ONLY edges to skip the the signal peptide.

    NOTE: This transforms the graph without returning it!

    Following Keys are set here:
    Nodes: <None>
    Edges: "qualifiers" ( -> adds SIGNAL)
    """
    if isinstance(signal_feature.location.end, UnknownPosition):
        # The Position of the end is not known. Therefore we skip
        # this entry simply. It does not contain any useful information
        return

    # Get end node
    [__stop_node__] = graph.vs.select(aminoacid="__end__")

    # Get start and end position of signal peptide
    # NOTE: + 1, since the start node occupies the position 0
    start_position, end_position = (
        signal_feature.location.start + 1,
        signal_feature.location.end + 0,
    )

    # Get all nodes with position start_position and their corresponding end_position(s)
    all_start_points = list(graph.vs.select(position=start_position))
    all_end_points = [
        _get_nodes_from_position(graph, x, end_position) for x in all_start_points
    ]

    # Check if only one end point exists for each start
    for x in all_end_points:
        if len(x) > 1:
            # This is unexpected and should not happen!
            raise UnexpectedException(
                accession=graph.vs[0]["accession"],
                position=signal_feature.location.end,
                message="Multiple End-Points for Signal-Feature defined. All end points are listed below",
                additional_info=str(all_end_points)
            )

    # TODO should the signal peptide be exactly the same as in canonical? Or can we leave it as is for isoforms?
    # Should we check this? If so: do this here!!
    # TODO can we associate the end points with the start points, according to its index?
    # TODO this can be probably summarized in generic (with INIT_MET!)

    # Create edge list
    all_edges = []
    for start_point, end_point_list in zip(all_start_points, all_end_points):
        for end_point in end_point_list:
            # For each start and end point

            # Get the corresponding edges
            edges_in = list(graph.es.select(_target=start_point))
            edges_out = list(graph.es.select(_source=end_point))

            # And add a new edge to skip the signal
            for ei in edges_in:
                for eo in edges_out:
                    all_edges.append(
                        (
                            (ei.source_vertex, eo.target_vertex),
                            [*_get_qualifiers(ei), signal_feature],
                        )
                    )

            # Special case the signal end can go directly to the stop node
            all_edges.append(((end_point, __stop_node__), [signal_feature]))

    # Bulk adding of edges into the graph
    cur_edges = graph.ecount()
    graph.add_edges([x[0] for x in all_edges])
    graph.es[cur_edges:]["qualifiers"] = [x[1] for x in all_edges]


def _get_nodes_from_position(graph, node, pos):
    """ Retrieve the nodes depending on isoform, whether it is set or not """
    if "isoform_accession" in node.attributes():
        return list(
            graph.vs.select(isoform_accession=node["isoform_accession"], position=pos)
        )
    else:
        return list(graph.vs.select(position=pos))

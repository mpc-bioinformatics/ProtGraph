from protgraph.unexpected_exception import UnexpectedException

from protgraph.ft_execution import _combine_vertices, _get_qualifiers, _get_all_vertices_before_after


def execute_variant(graph, variant_feature):
    """
    This function adds vertices and edges depending on the feature variant.
    Effectively, the variant information get executed. This method can skip
    feature variants which may reference ONLY isoforms. Skipped variants are printed.

    NOTE: This transforms the graph without returning it!

    Following Keys are set here:
    Nodes: <None>
    Edges: "qualifiers" ( -> adds VARIANT)
    """
    [__start_node__] = graph.vs.select(aminoacid="__start__")
    [__stop_node__] = graph.vs.select(aminoacid="__end__")
    # Get start and end position first
    # NOTE: Shifted by 1 due to the __start__ node beeing at 0
    aa_before = variant_feature.location.start + 1
    aa_after = variant_feature.location.end + 0

    # Get all vertices which are at beginning (before)
    # and all vertices which are at the end (after)
    vertices_before, vertices_after = _get_all_vertices_before_after(
        graph, aa_before, aa_after, variant_feature.ref
    )
    if len(vertices_before) == 0 or len(vertices_after) == 0:
        # Check if we have vertices, if not simply skip
        print("No Vertices retrieved for protein {}, using VARIANT: {} (referencing: {}). Skipping...".format(
            graph.vs[0]["accession"], variant_feature.id, variant_feature.ref))
        return

    # Now we check if we skip or add nodes
    text = variant_feature.qualifiers["note"]
    edge_list = []  # Here we append all edges, which should be added at the end
    if text.lower().startswith("missing"):
        # A sequence is missing! Just append an edge and its information
        # Get the aminoacid position before and after it
        # NOTE: Shifted by 1 due to the __start__ node at 0
        aa_before = variant_feature.location.start + 1
        aa_after = variant_feature.location.end + 0

        # TODO is such a combination enough?
        # Here we iterate over all possiblites over two pairs of nodes and its edges
        for aa_in_list, aa_out_list in zip(vertices_before, vertices_after):
            for aa_in in aa_in_list:
                for aa_edge_in in list(graph.es.select(_target=aa_in)):  # Get all incoming edges
                    for aa_out in aa_out_list:
                        for aa_edge_out in list(graph.es.select(_source=aa_out)):  # Get all outgoing edges
                            # Add corresponding edges and the qualifiers information
                            # for that edge (At least the VARIANT feature)
                            # Only if they not point to start_end directly! (#special case e.g. in P49782)
                            if aa_edge_in.source != __start_node__.index or aa_edge_out.target != __stop_node__.index:
                                edge_list.append(
                                    (
                                        (aa_edge_in.source, aa_edge_out.target),
                                        [*_get_qualifiers(aa_edge_in), variant_feature],
                                    )
                                )

    else:
        # A Sequence (or AA) is added
        # Get   X -> Y   Information
        # TODO duplicated in VAR_SEQ?
        idx = text.find("(")
        if idx != -1:
            text = text[:idx]
        xy = text.split("->")
        assert len(xy) == 2
        y = xy[1].strip().replace(" ", "")

        # TODO is such a combination enough?
        for aa_in_list, aa_out_list in zip(vertices_before, vertices_after):
            if len(aa_out_list) == 0 or len(aa_in_list) == 0:
                # Skip this entry, since we do not have complete information
                # -> Not possible to link either start or end or both
                continue

            # Add each individual amino acid as a node
            y_idcs = []
            for entry in y:
                vertex = graph.add_vertex()
                graph.vs[vertex.index]["aminoacid"] = entry
                graph.vs[vertex.index]["accession"] = graph.vs[1]["accession"]
                if "isoform_accession" in aa_in_list[0].attributes():
                    graph.vs[vertex.index]["isoform_accession"] = aa_in_list[0]["isoform_accession"]
                y_idcs.append(vertex.index)

            # Add edges between them (if needed)
            for idx, n in enumerate(y_idcs[:-1]):
                graph.add_edges([(n, y_idcs[idx + 1])])

            # Get the first and last node index
            first_node, last_node = y_idcs[0], y_idcs[-1]

            # And add then to the edges list (to connect them to the rest of the graph)
            for aa_in in aa_in_list:
                for aa_edge_in in list(graph.es.select(_target=aa_in)):  # Get all incoming edges
                    edge_list.append(
                        ((aa_edge_in.source, first_node), [*_get_qualifiers(aa_edge_in), variant_feature],)
                    )
            for aa_out in aa_out_list:
                for aa_edge_out in list(graph.es.select(_source=aa_out)):  # Get all outgoing edges
                    edge_list.append(
                        ((last_node, aa_edge_out.target), [])
                    )

    # Finally bulk add of the remaining edges
    cur_edges = graph.ecount()
    graph.add_edges([x[0] for x in edge_list])
    graph.es[cur_edges:]["qualifiers"] = [x[1] for x in edge_list]

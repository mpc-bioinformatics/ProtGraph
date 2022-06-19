from protgraph.ft_execution import _get_qualifiers, get_content
from protgraph.unexpected_exception import UnexpectedException


def execute_variant(graph, variant_feature):
    """ Wrapper for VARIANT """
    _execute_generic_feature(graph, variant_feature, "(")


def execute_mutagen(graph, mutagen_feature):
    """ Wrapper for MUTAGEN """
    _execute_generic_feature(graph, mutagen_feature, ":")


def execute_conflict(graph, conflict_feature):
    """ Wrapper for CONFLICT """
    _execute_generic_feature(graph, conflict_feature, "(")


def _execute_generic_feature(graph, generic_feature, beginning):
    """
    This function adds vertices and edges depending on the feature variant.
    Effectively, the variant information get executed. This method can skip
    features which may reference ONLY isoforms. Skipped features are printed (no thrown error).

    NOTE: This transforms the graph without returning it!

    Following Keys are set here:
    Nodes: <None>
    Edges: "qualifiers" ( -> adds VARIANT|MUTAGEN|CONFLICT)
    """
    # Get vertices before and after the chain (including "null"-chain)
    vertices_before, vertices_after = _get_vertices_before_after(graph, generic_feature)
    if vertices_after is None:
        return

    # Now we check if we skip or add nodes
    text = generic_feature.qualifiers["note"]
    edge_list = []  # Here we append all edges, which should be added at the end
    if text.lower().startswith("missing"):
        _append_edge_list_missing(
            graph, generic_feature, edge_list, vertices_before, vertices_after
        )
    else:
        _append_edge_list_chain(
            graph, text, generic_feature, edge_list, vertices_before, vertices_after, beginning=beginning
        )

    # Finally bulk add of the remaining edges
    cur_edges = graph.ecount()
    graph.add_edges([x[0] for x in edge_list])
    graph.es[cur_edges:]["qualifiers"] = [x[1] for x in edge_list]


def _get_vertices_before_after(graph, generic_feature):
    """ TODO DESC """
    # Get start and end position first
    # NOTE: Shifted by 1 due to the __start__ node beeing at 0
    aa_before = generic_feature.location.start + 1
    aa_after = generic_feature.location.end + 0

    # Get all vertices which are at beginning (before)
    # and all vertices which are at the end (after)
    vertices_before, vertices_after = _get_all_vertices_before_after(
        graph, aa_before, aa_after, generic_feature.ref
    )
    if len(vertices_before) == 0 or len(vertices_after) == 0:
        # Check if we have vertices, if not simply skip
        print("No Vertices retrieved for protein {}, using {}: {} (referencing: {}). Skipping...".format(
            graph.vs[0]["accession"], generic_feature.type, generic_feature.id, generic_feature.ref))
        return None, None

    return vertices_before, vertices_after


def _append_edge_list_chain(
    graph, text, generic_feature, edge_list, vertices_before, vertices_after, beginning="(", delimiter="->"
):
    # Get to be replaced amino_acids
    y_s = get_content(text, beginning, delimiter)

    # Generalizing: It can reference multiple substitutions
    for y in y_s.split(","):
        y = y.strip().replace(" ", "")
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
                        ((aa_edge_in.source, first_node), [*_get_qualifiers(aa_edge_in), generic_feature],)
                    )
            for aa_out in aa_out_list:
                for aa_edge_out in list(graph.es.select(_source=aa_out)):  # Get all outgoing edges
                    edge_list.append(
                        ((last_node, aa_edge_out.target), _get_qualifiers(aa_edge_out))
                    )


def _append_edge_list_missing(graph, generic_feature, edge_list, v_before, v_after):
    """ TODO DESC! """
    [__start_node__] = graph.vs.select(aminoacid="__start__")
    [__stop_node__] = graph.vs.select(aminoacid="__end__")
    # A sequence is missing! Just append an edge and its information

    # TODO is such a combination enough?
    # Here we iterate over all possiblites over two pairs of nodes and its edges
    for aa_in_list, aa_out_list in zip(v_before, v_after):
        for aa_in in aa_in_list:
            for aa_edge_in in list(graph.es.select(_target=aa_in)):  # Get all incoming edges
                for aa_out in aa_out_list:
                    for aa_edge_out in list(graph.es.select(_source=aa_out)):  # Get all outgoing edges
                        # Add corresponding edges and the qualifiers information
                        # for that edge (At least the generic feature)
                        # Only if they not point to start_end directly! (#special case e.g. in P49782)
                        if aa_edge_in.source != __start_node__.index or aa_edge_out.target != __stop_node__.index:
                            edge_list.append(
                                (
                                    (aa_edge_in.source, aa_edge_out.target),
                                    [*_get_qualifiers(aa_edge_in), generic_feature],
                                )
                            )


def _get_all_vertices_before_after(graph, aa_before: int, aa_after: int, reference: str):
    """
    Get the vertices which are at the beginning and end of the referencing feature.
    We explicitly check here if we need to take the isoform position (and accession)
    via the reference attribute, or if we simply query the graph for its position
    attribute.

    Returns two lists:
        before -> Aminoacid at the beginning
        after  -> Aminoacid at the end
    """
    if reference is None:
        # Get list of all aa before and after the specified position (no isoform)
        vertices_before_raw = list(graph.vs.select(position=aa_before))
        vertices_after_raw = list(graph.vs.select(position=aa_after))
    elif "isoform_accession" in graph.vs[0].attributes() and \
         "isoform_position" in graph.vs[0].attributes():
        # Get list of all aa before and after the specified position from the isoforms,
        # since this mutagen references explicitly a isoform!
        vertices_before_raw = list(graph.vs.select(isoform_position=aa_before, isoform_accession=reference))
        vertices_after_raw = list(graph.vs.select(isoform_position=aa_after, isoform_accession=reference))
    else:
        # No vertex found return empty lists.
        # this case happens e.g. if the user wants mutagen but no isoforms,
        # but since there are mutagen specifically referencing isoforms, we return
        # nothing
        vertices_before_raw, vertices_after_raw = [], []

    # TODO is such a combination enough and correct?
    # Combine the retrieved vertices to their corresponding isoforms (if available)
    vertices_before, vertices_after = _combine_vertices(vertices_before_raw, vertices_after_raw)

    return vertices_before, vertices_after


def _combine_vertices(list_a, list_b):
    """ Sorts vertices by their isoforms """
    # Sort the first list (incoming, a)
    out_d = {}
    for a in list_a:
        key = a["isoform_accession"] if "isoform_accession" in a.attributes() else None
        if key not in out_d:
            out_d[key] = dict(inn=[a])
        else:
            # We cannot simply associate via accession?
            # Does this case ever happen !?!
            raise UnexpectedException(
                accession=None,
                position=None,
                message="Exisiting Key would be overwritten. Additional Info contains dict, lista and listb",
                additional_info=str([out_d, list_a, list_b])
            )

    # Sort the second list (outgoing, b)
    for b in list_b:
        key = b["isoform_accession"] if "isoform_accession" in b.attributes() else None
        if key not in out_d:
            out_d[key] = dict(out=[b])
        else:
            if "out" not in out_d[key]:
                out_d[key]["out"] = [b]
            else:
                out_d[key]["out"].append(b)
            if len(out_d[key]) > 2:
                # TODO we cannot simply associate via accession!?!?
                # TODO custom exception
                raise Exception("multiple Entries Found")
    # TODO Does these two cases ever happen?

    # Generate the corresponding output lists
    # E.G.
    # a_s,          b_s
    # [PXXXX-1],    [PXXXX-1]
    # [PXXXX-2],    []
    # [],           [PXXXX-5]
    # [None],       [None]
    k = list(out_d.items())
    a_s = [x[1]["inn"] if "inn" in x[1] else [] for x in k]
    b_s = [x[1]["out"] if "out" in x[1] else [] for x in k]
    return a_s, b_s

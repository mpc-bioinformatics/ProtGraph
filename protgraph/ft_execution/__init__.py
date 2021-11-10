def _get_all_vertices_before_after(graph, aa_before: int, aa_after: int, reference: str):
    """
    Get the vertices which are at the beginning and end of the referencing mutagen.
    We explicitly check here if we need to take the isoform position (and accesion)
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


def _get_qualifiers(edge):
    """ A simple method to retrieve qualifiers. It always returns a list """
    qualifiers = edge["qualifiers"] if "qualifiers" in edge.attributes() else []
    qualifiers = [] if qualifiers is None else qualifiers
    return qualifiers

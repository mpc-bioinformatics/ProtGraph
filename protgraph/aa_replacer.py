def replace_aa(graph_entry, substitutions: list):
    """
    Replaces amino acids to other amino acids.

    E.G: Set Substitutions to [(X,  [Y])] (or from CLI via X -> Y)
    will replace all aminoacid with 'X' in the graph with 'Y'

    Multiple Substitutions are allowed: [(X,  [Y, Z])] (or from CLI via X -> Y,Z)
    which then replaces 'X' with 'Y' and additionally 'X' with 'Z'.

    Replacements can also be chained. E.G.: J -> I   and I -> L (as list: [(J, [I]), (I, [L])]
    would first replace all amino acids with 'J' to 'I' afterwared replacing all 'I' with 'L'.
    Effectively, this example replaces J -> L.
    """
    # Check if we have something to replace
    if substitutions is None:
        return  # No replacement given, we return simply

    # For each replacement
    for source, targets in substitutions:
        # Set the list for the removable edges
        edges_to_be_removed = []

        # Get all amino acids which should be replaced
        vertices_to_replace = graph_entry.vs.select(aminoacid=source)
        for v in vertices_to_replace:
            # Get all in and out edges
            i_edges = [x.index for x in v.in_edges()]
            o_edges = [x.index for x in v.out_edges()]
            # And also add them explicitly here to be removed later
            edges_to_be_removed.extend(
                [(graph_entry.es[e].source, graph_entry.es[e].target) for e in i_edges] +
                [(graph_entry.es[e].source, graph_entry.es[e].target) for e in o_edges]
            )

            # Add new nodes depending on the number of replacement candidates
            vc = graph_entry.vcount()
            graph_entry.add_vertices(len(targets))

            # Set the attributes of the new node(s) (bulk)
            graph_entry.vs[vc:]["aminoacid"] = targets
            graph_entry.vs[vc:]["position"] = [v["position"]]*len(targets)
            graph_entry.vs[vc:]["accession"] = [v["accession"]]*len(targets)
            # If the graph (or vertex) has isoform-information add it here too!
            if "isoform_accession" in v.attributes():
                graph_entry.vs[vc:]["isoform_accession"] = [v["isoform_accession"]]*len(targets)
            if "isoform_position" in v.attributes():
                graph_entry.vs[vc:]["isoform_position"] = [v["isoform_position"]]*len(targets)

            # Generate list of edges to be added (we use all possible combinations)
            edges_to_add = [(graph_entry.es[e].source, x.index) for x in graph_entry.vs[vc:] for e in i_edges] +\
                [(x.index, graph_entry.es[e].target) for x in graph_entry.vs[vc:] for e in o_edges]

            # Add edges and (bulk) add their corresponding qualifier information if needed
            ec = graph_entry.ecount()
            graph_entry.add_edges(edges_to_add)
            if "qualifiers" in graph_entry.es[0].attributes():
                # Generate an additional edges_attribute list for qualifiers
                edges_attrs_to_add = [graph_entry.es[e]["qualifiers"] for x in graph_entry.vs[vc:] for e in i_edges] +\
                    [graph_entry.es[e]["qualifiers"] for x in graph_entry.vs[vc:] for e in o_edges]
                graph_entry.es[ec:]["qualifiers"] = edges_attrs_to_add

        # Finally delete all the edges and nodes which we have replaced by the new nodes!
        graph_entry.delete_edges(edges_to_be_removed)
        graph_entry.delete_vertices(vertices_to_replace)  # We also remove the vertices which now got replaced

from collections import defaultdict

from protgraph.unexpected_exception import UnexpectedException


class Or(list):
    """
    This is a wrapper class over the type list.

    If we want to collapse edges, then we need to have construct,
    which shows that both options are valid to get to the other node.

    This Or class represents this. If something like
    Or[[Feat1], [Feat2, Feat3]] is encountered on a edge, it indicates that either by
    using the Feat1 OR (Feat1 and then Feat3) can go to the same node.
    """
    def __repr__(self):
        return "Or" + super().__repr__()


def collapse_parallel_edges(graph):
    """
    Collapses edges. This has a direct connection to the issue:
    https://github.com/mpc-bioinformatics/ProtGraph/issues/41

    In some rare cases (currently ~10 Proteins in UniProt) parallel edges are drawn due to their
    annotated features. These parallel edges are NOT wrongly added but simply show that it is
    possible to use different features while "connecting" to the same target node.

    TODO we may need to sort Variants, Mutagens and Conflicts by the starting point in ascending order,
    since it may happen that if they are applied unordered, that an parallel edge is never drawn.
    (tested on Q9R2E6, 1 out of 6 order possibilities did not yield the parallel edge! (132 for reference))
    """
    # Go through each node
    for x in range(0, graph.vcount()):
        k_edges = graph.vs[x].out_edges()
        k_targets = [x.target for x in k_edges]
        k_set = set(k_targets)
        # And check if any of the outgoing edges go to the same node
        if len(k_targets) != len(k_set):

            # If that is the case, we resort them into bins
            # so we do something more sophisticated only iff there are
            # parallel edges
            dups = defaultdict(list)
            for e, t in zip(k_edges, k_targets):
                dups[t].append(e)

            # And iterate over the bins. If the size is larger than 1 we replace an edge
            remove_dup_edges = []
            for edge_list in dups.values():
                if len(edge_list) > 1:

                    if "cleaved" in graph.es[0].attributes():
                        # Set cleavage information of collapsed edge
                        s_cleaved = set(x["cleaved"] for x in edge_list)
                        if len(s_cleaved) != 1:
                            # This should not happen!
                            raise UnexpectedException(
                                accession=graph.vs[0]["accession"],
                                position=graph.vs[edge_list[0].source]["position"],
                                message="Different cleavage information in edges, which would be collapsed found.",
                                additional_info=str(edge_list)
                            )
                        cleaved = s_cleaved.pop()

                    if "qualifiers" in graph.es[0].attributes():
                        # Set qualfiers information of collapsed edge
                        collapsed_qualifiers = [x["qualifiers"] for x in edge_list]
                        if collapsed_qualifiers:
                            # Sometimes we may need to compare and summarize
                            # qualifiers. It can happen that qualfiers
                            # like: "Or[VAR_012345|VAR_012345]" are generated,
                            # which we would want to simplify to: "VAR_012345"
                            qualifiers = _collapse_qualifier_information(collapsed_qualifiers)
                        else:
                            qualifiers = None

                    # Set collapsed edge and remove other edges from graph
                    e_to_keep = edge_list[0].index
                    if "cleaved" in graph.es[0].attributes():
                        graph.es[e_to_keep]["cleaved"] = cleaved
                    if "qualifiers" in graph.es[0].attributes():
                        graph.es[e_to_keep]["qualifiers"] = qualifiers

                    # Add remaining edges to remove list
                    remove_dup_edges.extend([x.index for x in edge_list[1:]])

            # We finally remove all edges at once for the current node we looked into
            graph.delete_edges(remove_dup_edges)


def _collapse_qualifier_information(qualifiers_list):
    """
    Collapses qualifier information by iterating twice over
    the list and saving only one of the duplicated qualifiers.
    Returns either the qualifier itself (if only one remains)
    or an Or-Object.

    NOTE: the qualfier_list, MUST be at least of lenght 2!
    """
    include_qualifiers = []
    already_added = [False]*len(qualifiers_list)

    for i in range(len(qualifiers_list)-1):
        # Skip entry if we already added it
        if already_added[i]:
            continue

        for j in range(i+1, len(qualifiers_list)):
            # Skip entry if we already added it
            if already_added[j]:
                continue
            if qualifiers_list[i] is None or qualifiers_list[j] is None:
                if qualifiers_list[i] == qualifiers_list[j]:
                    already_added[j] = True
                    if not already_added[i]:
                        include_qualifiers.append(qualifiers_list[i])
                        already_added[i] = True
                continue

            if len(qualifiers_list[i]) == len(qualifiers_list[j]) and \
                all([
                    True if x == y else None
                    for x, y in zip(qualifiers_list[i], qualifiers_list[j])
                    ]):

                already_added[j] = True
                if not already_added[i]:
                    include_qualifiers.append(qualifiers_list[i])
                    already_added[i] = True

    for x, y in zip(qualifiers_list, already_added):
        if not y:
            include_qualifiers.append(x)

    if len(include_qualifiers) == 1:
        return include_qualifiers[0]
    else:
        return [Or(include_qualifiers)]

from collections import defaultdict


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
                    s_cleaved = set(x["cleaved"] for x in edge_list)
                    if len(s_cleaved) != 1:
                        # TODO custom exception, this should not happen!
                        # track information like key, which node/edge especially
                        # retrieve the accession explicitly here too for debugging!
                        raise Exception("Multiple entries in set!")

                    # Get and set all attributes in edges on the first edge
                    cleaved = s_cleaved.pop()
                    qualifiers = [Or([x["qualifiers"] for x in edge_list])]

                    e_to_keep = edge_list[0].index
                    graph.es[e_to_keep]["cleaved"] = cleaved
                    graph.es[e_to_keep]["qualifiers"] = qualifiers

                    # Add remaining edges to remove list
                    remove_dup_edges.extend([x.index for x in edge_list[1:]])

            # We finally remove all edges at once for the current node we looked into
            graph.delete_edges(remove_dup_edges)

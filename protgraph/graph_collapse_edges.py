from collections import defaultdict

class Or(list):
    """
    TODO DODODODODOD
    """
    def __repr__(self):
        return "Or" + super().__repr__()




def collapse_parallel_edges(graph):

    # Go through each node
    for x in range(0, graph.vcount()):
        k_edges = graph.vs[x].out_edges()
        k_targets = [x.target for x in k_edges]
        k_set = set(k_targets)
        # and check if any of the outgoing edges go to the same node
        if len(k_targets) != len(k_set):
            
            # if that is the case, we resort them into bins
            dups = defaultdict(list)
            for e, t in zip(k_edges, k_targets):
                dups[t].append(e)

            # and iterate over the bins. If the size larger than 1 we replace an edge
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

            graph.delete_edges(remove_dup_edges)


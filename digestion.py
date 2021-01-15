def digest(graph, enzyme: str):
    """
    Select digestion method depending on enzyme

    Each digestion returns the number of cleaved edges (marked)

    Following Keys are set here (except for skip):
    Nodes: <None>
    Edges: "cleaved" ( -> either True or False)
    """
    return dict(
        skip=_digest_via_skip,
        trypsin=_digest_via_trypsin
        # Add more Enzymes if needed here!
    )[enzyme](graph)


def _digest_via_skip(graph):
    """ Skipping digestion """
    return 0


def _digest_via_trypsin(graph):
    """
    Digestion via Trypsin.

    Each edge from a node with the aminoacid K or R as source gets cleaved (marked as TRUE)
    except when a P is followed (-> set as target).

    Additionally two new edges are added:
    1: One edge from __start__ to the next nodes target (Beginning of a peptide)
    2: One edge from K or R to __end__ and (Ending of a peptide)

    NOTE: Multiple edges can go in or out from K and R and are therefore also considered.
    Also, 'infinitly' many miscleavages are represented!
    """

    # Get start and end node
    [__start_node__] = graph.vs.select(aminoacid="__start__")
    [__end_node__] = graph.vs.select(aminoacid="__end__")

    # Get all aminoacid edges where source is K or R
    k_s = graph.vs.select(aminoacid="K")
    r_s = graph.vs.select(aminoacid="R")
    k_s_edges = [graph.es.select(_source=k) for k in k_s]
    r_s_edges = [graph.es.select(_source=r) for r in r_s]

    # Filter out edges, which have P as target
    k_s_edges_remaining = [
        y for x in k_s_edges for y in x if graph.vs[y.target]["aminoacid"] != "P"
    ]
    r_s_edges_remaining = [
        y for x in r_s_edges for y in x if graph.vs[y.target]["aminoacid"] != "P"
    ]

    # Mark edges as cleaved
    cleaved_idcs = [x.index for x in k_s_edges_remaining + r_s_edges_remaining]
    graph.es[cleaved_idcs]["cleaved"] = True

    # Digest the graph into smaller parts
    # Explicitly we add ->
    trypsin_in = [
        (__start_node__, k.target) for k in k_s_edges_remaining + r_s_edges_remaining
    ]  # edges for Nodes which should have an edge to __start__
    trypsin_out = [
        (k.source, __end_node__) for k in k_s_edges_remaining + r_s_edges_remaining
    ]  # edges for Nodes which should have an edge to __end__

    # Add the newly created edges to the graph
    graph.add_edges(trypsin_in + trypsin_out)

    # Return the number of cleaved edges
    return len(cleaved_idcs)

from typing import List


def digest(graph, enzymes: List[str]):
    """
    Select digestion method depending on enzyme

    Each digestion returns the number of cleaved edges (marked)

    Following Keys are set here (except for skip):
    Nodes: <None>
    Edges: "cleaved" ( -> either True or False/None)
    """
    if enzymes is None:  # Default is Trypsin, therefor we add it
        enzymes = ["trypsin"]

    sum_of_cleavages = 0
    for i in enzymes:
        sum_of_cleavages += dict(
            skip=_digest_via_skip,
            trypsin=_digest_via_trypsin,
            full=_digest_via_full
            # Add more Enzymes if needed here!
        )[i](graph)

    return sum_of_cleavages


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
        y for x in k_s_edges for y in x
        if graph.vs[y.target]["aminoacid"] != "P" and graph.vs[y.target]["aminoacid"] != "__end__"
    ]
    r_s_edges_remaining = [
        y for x in r_s_edges for y in x
        if graph.vs[y.target]["aminoacid"] != "P" and graph.vs[y.target]["aminoacid"] != "__end__"
    ]

    # Mark edges as cleaved
    cleaved_idcs = [x.index for x in k_s_edges_remaining + r_s_edges_remaining]
    graph.es[cleaved_idcs]["cleaved"] = True

    # Digest the graph into smaller parts
    # Explicitly we add ->
    trypsin_in = [
        (__start_node__.index, k.target) for k in k_s_edges_remaining + r_s_edges_remaining
    ]  # edges for Nodes which should have an edge to __start__
    trypsin_out = [
        (k.source, __end_node__.index) for k in k_s_edges_remaining + r_s_edges_remaining
    ]  # edges for Nodes which should have an edge to __end__

    # Check if edge already exists, if so skip it:
    remaining_edges = []
    for x in set(trypsin_out).union(set(trypsin_in)):
        try:
            # Skip if found
            graph.es.find(_between=((x[0],), (x[1],)))
        except ValueError:
            remaining_edges.append(x)

    # Add the newly created edges to the graph
    graph.add_edges(remaining_edges)

    # Return the number of cleaved edges
    return len(cleaved_idcs)


def _digest_via_full(graph):
    """
    Digestion via Full. Here we digest at every possible and available edge in the graph.

    First, retrieve all possible Nodes which DO NOT have an direct edge to either
    the start or end node. Then add edges inbetween them, so that they have an direct
    edge from start AND to end.
    """
    # Get start and end node
    [__start_node__] = graph.vs.select(aminoacid="__start__")
    [__end_node__] = graph.vs.select(aminoacid="__end__")

    # Get all Nodes which should have an IN edge from start
    start_in = graph.vs.indices
    # Remove all nodes which already have an edge (and are start and end itself)
    start_in.remove(__start_node__.index)
    start_in.remove(__end_node__.index)
    for i in graph.neighbors(__start_node__, mode="OUT"):
        try:
            start_in.remove(i)
        except Exception:
            print(
                "WARNING: Graph has parallel edges, which may be due to the input file. "
                "Please check that features are not repeated! (Entry: {})".format(graph.vs[0]["accession"])
            )

    # Do the same for end to get all nodes which need OUT edges
    end_out = graph.vs.indices
    # As above, remove all nodes to end
    end_out.remove(__start_node__.index)
    end_out.remove(__end_node__.index)
    for i in graph.neighbors(__end_node__, mode="IN"):
        try:
            end_out.remove(i)
        except Exception:
            print(
                "WARNING: Graph has parallel edges, which may be due to the input file. "
                "Please check that features are not repeated! (Entry: {})".format(graph.vs[0]["accession"])
            )

    # Get all edges, which should be cleaved and mark them
    cleaved_edges = [
        x.index
        for x in graph.es[:]
        if x.source != __start_node__.index and
        x.target != __end_node__.index
    ]
    graph.es[cleaved_edges]["cleaved"] = True

    # Generate the edges, which should be added:
    start_in_edges = [(__start_node__.index, x) for x in start_in]
    end_out_edges = [(x, __end_node__.index) for x in end_out]

    # Add the edges to the graph
    graph.add_edges(start_in_edges + end_out_edges)

    # Return the number of cleaved edges
    return len(cleaved_edges)

from typing import List


def digest(graph, enzymes: List[str], entry_dict):
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
        sum_of_cleavages += DIGESTION_MAP[i](graph)

    entry_dict["num_cleavages"] = sum_of_cleavages


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
    if "qualifiers" in graph.es[0].attributes():
        qualifiers_info = [
            k["qualifiers"] for k in k_s_edges_remaining + r_s_edges_remaining
        ] + [None]*len(cleaved_idcs)

    # Add the newly created edges to the graph
    e_count = graph.ecount()
    graph.add_edges(trypsin_in + trypsin_out)
    if "qualifiers" in graph.es[0].attributes():
        graph.es[e_count:]["qualifiers"] = qualifiers_info

    # Return the number of cleaved edges
    return len(cleaved_idcs)


def _digest_via_glu_c(graph):
    """
    Digestion via Glu-C.

    Each edge from a node with the aminoacid D or E as source gets cleaved (marked as TRUE)
    except when a P is followed (-> set as target).

    Additionally two new edges are added:
    1: One edge from __start__ to the next nodes target (Beginning of a peptide)
    2: One edge from D or E to __end__ and (Ending of a peptide)

    NOTE: Multiple edges can go in or out from K and R and are therefore also considered.
    Also, 'infinitly' many miscleavages are represented!
    """

    # Get start and end node
    [__start_node__] = graph.vs.select(aminoacid="__start__")
    [__end_node__] = graph.vs.select(aminoacid="__end__")

    # Get all aminoacid edges where source is K or R
    d_s = graph.vs.select(aminoacid="D")
    e_s = graph.vs.select(aminoacid="E")
    d_s_edges = [graph.es.select(_source=d) for d in d_s]
    e_s_edges = [graph.es.select(_source=e) for e in e_s]

    # Filter out edges, which have P as target
    d_s_edges_remaining = [
        y for x in d_s_edges for y in x
        if graph.vs[y.target]["aminoacid"] != "P" and graph.vs[y.target]["aminoacid"] != "__end__"
    ]
    e_s_edges_remaining = [
        y for x in e_s_edges for y in x
        if graph.vs[y.target]["aminoacid"] != "P" and graph.vs[y.target]["aminoacid"] != "__end__"
    ]

    # Mark edges as cleaved
    cleaved_idcs = [x.index for x in d_s_edges_remaining + e_s_edges_remaining]
    graph.es[cleaved_idcs]["cleaved"] = True

    # Digest the graph into smaller parts
    # Explicitly we add ->
    gluc_in = [
        (__start_node__.index, e.target) for e in d_s_edges_remaining + e_s_edges_remaining
    ]  # edges for Nodes which should have an edge to __start__
    gluc_out = [
        (e.source, __end_node__.index) for e in d_s_edges_remaining + e_s_edges_remaining
    ]  # edges for Nodes which should have an edge to __end__
    if "qualifiers" in graph.es[0].attributes():
        qualifiers_info = [
            k["qualifiers"] for k in d_s_edges_remaining + e_s_edges_remaining
        ] + [None]*len(cleaved_idcs)

    # Add the newly created edges to the graph
    e_count = graph.ecount()
    graph.add_edges(gluc_in + gluc_out)
    if "qualifiers" in graph.es[0].attributes():
        graph.es[e_count:]["qualifiers"] = qualifiers_info

    # Return the number of cleaved edges
    return len(cleaved_idcs)


def _digest_via_full(graph):
    # TODO DL add qualfiers!!
    """
    Digestion via Full. Here we digest at every possible and available edge in the graph.

    First, retrieve all possible Nodes which DO NOT have an direct edge to either
    the start or end node. Then add edges inbetween them, so that they have an direct
    edge from start AND to end.
    """
    # Get start and end node
    [__start_node__] = graph.vs.select(aminoacid="__start__")
    [__end_node__] = graph.vs.select(aminoacid="__end__")

    # Prune start and end
    all_edges = list(graph.es[:])
    all_edges_wo_start = [x for x in all_edges if x.source != __start_node__.index]
    all_edges_wo_start_end = [x for x in all_edges_wo_start if x.target != __end_node__.index]

    # Mark cleaved edges
    graph.es[[x.index for x in all_edges_wo_start_end]]["cleaved"] = True

    # Generate the edges, which should be added:
    start_in_edges = [(__start_node__.index, x.target) for x in all_edges_wo_start_end]
    end_out_edges = [(x.source, __end_node__.index) for x in all_edges_wo_start_end]
    if "qualifiers" in graph.es[0].attributes():
        qualifiers = [x["qualifiers"] for x in all_edges_wo_start_end] + [None] * len(all_edges_wo_start_end)

    # Add the edges to the graph
    # Add the newly created edges to the graph
    e_count = graph.ecount()
    graph.add_edges(start_in_edges + end_out_edges)
    if "qualifiers" in graph.es[0].attributes():
        graph.es[e_count:]["qualifiers"] = qualifiers

    # Return the number of cleaved edges
    return len(all_edges_wo_start_end)


DIGESTION_MAP = dict(
    skip=_digest_via_skip,
    trypsin=_digest_via_trypsin,
    gluc=_digest_via_glu_c,
    full=_digest_via_full
    # Add more Enzymes if needed here!
)

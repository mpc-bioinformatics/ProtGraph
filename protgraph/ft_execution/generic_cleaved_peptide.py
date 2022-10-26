from Bio.SeqFeature import UnknownPosition

from protgraph.ft_execution import _get_qualifiers
from protgraph.ft_execution.generic import _get_all_vertices_before_after


def execute_propeptide(graph, propeptide_feature):
    """ Wrapper function to execute propeptide_features"""
    execute_generic_cleaved_peptide(graph, propeptide_feature)


def execute_peptide(graph, peptide_feature):
    """ Wrapper function to execute peptide_features"""
    # NOTE: a few annotations in "proteins" reference the start and end
    # of it completely, making parallel edges at start or ending?
    # Are these "Proteins" in UniProt actually correct to be annotated as such?
    # These are actually only pepitdes...
    # E.G.: P85078 or P37046 (there are in total a few thousands probably)
    execute_generic_cleaved_peptide(graph, peptide_feature)


def execute_chain(graph, chain_feature):
    """ Wrapper function to execute chain_features"""
    # TODO
    execute_generic_cleaved_peptide(graph, chain_feature)


def execute_generic_cleaved_peptide(graph, generic_cleaved_feature):
    """
    This function adds ONLY edges to cleave the propeptide.

    NOTE: This transforms the graph without returning it!

    Following Keys are set here:
    Nodes: <None>
    Edges: "qualifiers" ( -> adds PROPEP|PEPTIDE)
    """
    if isinstance(generic_cleaved_feature.location.end, UnknownPosition) or \
       isinstance(generic_cleaved_feature.location.start, UnknownPosition):
        # The Position of the end or start is not known. Therefore we skip
        # this entry simply. It does not contain any useful information
        return

    # Get start and end node
    [__start_node__] = graph.vs.select(aminoacid="__start__")
    [__stop_node__] = graph.vs.select(aminoacid="__end__")

    # Get start and end position of signal peptide
    # NOTE: + 1, since the start node occupies the position 0
    start_position, end_position = (
        generic_cleaved_feature.location.start + 1,
        generic_cleaved_feature.location.end + 0,
    )

    # Get the corrsponding start and end nodes of the referenced peptide
    # (including isoforms, if not specified or only isoforms)
    start_nodes, end_nodes = _get_all_vertices_before_after(
        graph, start_position, end_position, generic_cleaved_feature.ref
    )

    # Create the edge-list (cleaving the referenced peptide)
    edge_list, edge_fts = _create_edges_list_and_feature(
        start_nodes, end_nodes, __start_node__.index, __stop_node__.index, generic_cleaved_feature
    )

    # Bulk adding of edges into the graph
    cur_edges = graph.ecount()
    graph.add_edges(edge_list)
    graph.es[cur_edges:]["qualifiers"] = edge_fts


def _create_edges_list_and_feature(start_nodes, end_nodes, start_idx, stop_idx, feature):
    """
    Cleaves the ingoing edges of start_nodes and the outgoing edges of end_nodes.

    This function returns two lists containing the edges and the features,
    which can be later applied on the graph itself
    """
    # Create the edge-list (cleaving the referenced peptide)
    edge_list = []
    edge_feature = []

    # Cleave the beginning
    for sn_iso in start_nodes:
        for sn in sn_iso:
            edge_list.append((start_idx, sn.index))  # Add Ingoing edge from start
            edge_feature.append([feature])
            for ie in sn.in_edges():
                if ie.source != start_idx:  # Check if ingoing node is start
                    edge_list.append((ie.source, stop_idx))  # Add outgoing edge from all in-going nodes
                    edge_feature.append([*_get_qualifiers(ie), feature])

    # Cleave the end
    for en_iso in end_nodes:
        for en in en_iso:
            edge_list.append((en.index, stop_idx))  # Add Outgoing edge to end
            edge_feature.append([feature])
            for oe in en.out_edges():
                if oe.target != stop_idx:  # Check if ingoing node is start
                    edge_list.append((start_idx, oe.target))  # Add outgoing edge from all in-going nodes
                    edge_feature.append([*_get_qualifiers(oe), feature])

    return edge_list, edge_feature

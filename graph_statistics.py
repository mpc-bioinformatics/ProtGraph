def get_statistics(graph, **kwargs):
    """
        TODO

        returns #Node, #Edges, #Num_Of_Paths
    """

    # Get the number of nodes and edges (an be done instantly)
    num_edges = _get_edge_count(graph)
    num_nodes = _get_node_count(graph)

    # Get the number of possible paths if set
    num_possible_paths = _num_of_possible_paths(graph) if kwargs["calc_num_possibilities"] else None

    # TODO can we calculate more statistics?

    # Return all information
    return num_nodes, num_edges, num_possible_paths


def _get_edge_count(graph_entry):
    """ Get number of edges"""
    return graph_entry.ecount()


def _get_node_count(graph_entry):
    """ Get number of nodes"""
    return graph_entry.vcount()


def _num_of_possible_paths(graph_entry):
    """ Get the Number of all possible simple Paths for a Protein or Peptide.
        A dynamic programming approach is taken here. We can minimize this problem
        into subproblems. The goal is to find the number of paths from the start of
        a protein to the end of a protein. 

        This can be divided into the number of paths from the start to a node in the protein,
        which sums the number of possible paths from its previous nodes. This can be continued up to
        the end node, yielding the number of possible (non-repetative) paths to the end of a protein.

        This algorithm therefore needs to iterate over the graph, which can be done with the help of 
        the topological sort. The runtime of the dynamic programming part should be O(n^2) TODO is this 100% correct?
    """

    # First get topological sorting of the graph
    sorted_nodes = graph_entry.topological_sorting()

    # Create list with num of paths
    var_paths = [0] * graph_entry.vcount()

    # Initialize Path from the very first node! For convenience we set it to one (actually 0!)
    # The very first node should always be the __start__ of a protein
    first = sorted_nodes[0]
    var_paths[first] = 1

    # Iterative approach look how many paths are possible from previous to itself (O(n^2)) # TODO is this 100% correct?
    # Get next node in topological sorted nodes
    for v in sorted_nodes[1:]:
        # Get the sum of possible paths from it, looking at each node which points to it
        summed = 0
        for v_prev in graph_entry.neighbors(v, mode="IN"):
            summed += var_paths[v_prev]

        # Set its number of paths as the calculated sum
        var_paths[v] = summed

    # We changed the first nodes value to 1, which is not correct, we set it back here!
    var_paths[first] = 0  # Path to itself is zero!

    # Returning the last node (which should always be the __end__ of a protein)
    return var_paths[sorted_nodes[-1]]  # This contains the number of paths

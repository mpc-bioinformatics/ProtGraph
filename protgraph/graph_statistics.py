from operator import add

from protgraph.graph_collapse_edges import Or

def get_statistics(graph, **kwargs):
    """
    TODO can we retrieve even more information!?
    returns #Node, #Edges, #Num_Of_Paths
    """

    # Get the number of nodes and edges (an be done instantly)
    num_edges = _get_edge_count(graph)
    num_nodes = _get_node_count(graph)

    # Get the number of possible paths if set
    num_possible_paths = (_num_of_possible_paths(graph) if kwargs["calc_num_possibilities"] else None)

    # Get the number of possible paths (with all miscleavages) if set
    num_possible_paths_all_mis = (
        _num_of_possible_paths_all_miscleavages(graph)
        if kwargs["calc_num_possibilities_miscleavages"] else None
    )

    # Get the number of possible paths (with number of hops) if set
    num_possible_paths_all_hops = (
        _num_of_possible_paths_all_hops(graph)
        if kwargs["calc_num_possibilities_hops"] else None
    )

    # Get the number of possible paths by specific features
    num_possible_paths_variant  = (
        _num_of_possible_paths_feature_type(
            graph, feature_type="VARIANT", 
            or_count=kwargs["calc_num_possibilites_or_count"]
        )
        if kwargs["calc_num_possibilities_variant"] else None
    )
    num_possible_paths_mutagen  = (
        _num_of_possible_paths_feature_type(
            graph, feature_type="MUTAGEN", 
            or_count=kwargs["calc_num_possibilites_or_count"]
        )
        if kwargs["calc_num_possibilities_mutagen"] else None
    )
    num_possible_paths_conflict = (
        _num_of_possible_paths_feature_type(
            graph, feature_type="CONFLICT", 
            or_count=kwargs["calc_num_possibilites_or_count"]
        )
        if kwargs["calc_num_possibilities_conflict"] else None
    )

    # TODO can we calculate more statistics?

    # Return all information
    return num_nodes, num_edges, num_possible_paths, \
        num_possible_paths_all_mis, num_possible_paths_all_hops, \
        num_possible_paths_variant, num_possible_paths_mutagen, \
        num_possible_paths_conflict


def _get_edge_count(graph_entry):
    """ Get number of edges"""
    return graph_entry.ecount()


def _get_node_count(graph_entry):
    """ Get number of nodes"""
    return graph_entry.vcount()


def _num_of_possible_paths(graph_entry):
    """
    Get the Number of all possible simple Paths for a Protein or Peptide.
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

    # Iterative approach look how many paths are possible from previous to itself (O(n^2))
    # TODO is the runtime 100% correct?
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


def _num_of_possible_paths_all_miscleavages(graph_entry):
    """
    Get the Number of all possible simple Paths with cleavages for a Protein or Peptide.
    A dynamic programming approach is taken here. We can minimize this problem
    into subproblems. The goal is to find the number of paths from the start of
    a protein to the end of a protein.

    This can be divided into the number of paths from the start to a node in the protein,
    which sums the number of possible paths from its previous nodes. This can be continued up to
    the end node, yielding the number of possible (non-repetative) paths to the end of a protein.

    For counting the miscleavages, we simply use a list instead of an value and shift, if an edge is
    an cleaved one. NOTE: Due to this approach, the calculation can be memory heavy!

    This algorithm therefore needs to iterate over the graph, which can be done with the help of
    the topological sort. The runtime of the dynamic programming part should be O(n^2) TODO is this 100% correct?
    """
    if "cleaved" not in graph_entry.es[0].attributes():
        # Case: we do not have a property to execute this statistic
        return None

    # First get topological sorting of the graph
    sorted_nodes = graph_entry.topological_sorting()

    # Create list with num of paths as LISTS
    var_paths = [[]] * graph_entry.vcount()

    # Initialize Path from the very first node! For convenience we set it to one (actually 0!)
    # The very first node should always be the __start__ of a protein
    first = sorted_nodes[0]
    var_paths[first] = [1]  # as LIST entry!

    # Iterative approach look how many paths are possible from previous to itself (O(n^2))
    # TODO is the runtime 100% correct?
    # Get next node in topological sorted nodes
    for v in sorted_nodes[1:]:
        # Get the sum of possible paths from it, looking at each node which points to it
        summed = []
        for e_in in graph_entry.vs[v].in_edges():
            # This is the only change in the algorithm
            if e_in["cleaved"]:
                # if cleaved then shift by one!
                summed = _add_lists(summed, [0] + var_paths[e_in.source])
            else:
                # simply add lists
                summed = _add_lists(summed, var_paths[e_in.source])

        # Set its number of paths as the calculated sum
        var_paths[v] = summed

    # We changed the first nodes value to 1, which is not correct, we set it back here!
    var_paths[first] = [0]  # Path to itself is zero! (lastly as LIST entry)

    # Returning the last node (which should always be the __end__ of a protein)
    # Here the index of the list gives us the number of how many cleavages we have missed.
    # Sum each element up to retrieve the number of all possible paths ("infinite" many miscleavages)
    return var_paths[sorted_nodes[-1]]


def _num_of_possible_paths_all_hops(graph_entry):
    """
    Get the Number of all possible simple Paths with the number of hops for a Protein or Peptide.
    A dynamic programming approach is taken here. We can minimize this problem
    into subproblems. The goal is to find the number of paths from the start of
    a protein to the end of a protein.

    This can be divided into the number of paths from the start to a node in the protein,
    which sums the number of possible paths from its previous nodes. This can be continued up to
    the end node, yielding the number of possible (non-repetative) paths to the end of a protein.

    For counting the hops, we use a similar approach as in counting the miscleavages.
    NOTE: Due to this approach, the calculation can be memory heavy!

    This can be very usefull to determine the depth for each graph and at which depth
    the number of possible peptides explode.

    This algorithm therefore needs to iterate over the graph, which can be done with the help of
    the topological sort. The runtime of the dynamic programming part should be O(n^2) TODO is this 100% correct?
    """
    # First get topological sorting of the graph
    sorted_nodes = graph_entry.topological_sorting()

    # Create list with num of paths as LISTS
    var_paths = [[]] * graph_entry.vcount()

    # Initialize Path from the very first node! For convenience we set it to one (actually 0!)
    # The very first node should always be the __start__ of a protein
    first = sorted_nodes[0]
    var_paths[first] = [1]  # as LIST entry!

    # Iterative approach look how many paths are possible from previous to itself (O(n^2))
    # TODO is the runtime 100% correct?
    # Get next node in topological sorted nodes
    for v in sorted_nodes[1:-1]:
        # Get the sum of possible paths from it, looking at each node which points to it
        summed = []
        for e_in in graph_entry.vs[v].in_edges():
            # Always shift by one, except for the last element in list
            summed = _add_lists(summed, [0] + var_paths[e_in.source])

        # Set its number of paths as the calculated sum
        var_paths[v] = summed
    # CASE last element in list
    summed = []
    for e_in in graph_entry.vs[sorted_nodes[-1]].in_edges():
        # Always shift by one, except for the last element in list
        summed = _add_lists(summed, var_paths[e_in.source])
    var_paths[sorted_nodes[-1]] = summed

    # We changed the first nodes value to 1, which is not correct, we set it back here!
    var_paths[first] = [0]  # Path to itself is zero! (lastly as LIST entry)

    # Returning the last node (which should always be the __end__ of a protein)
    # Here the index of the list gives us the number of how many cleavages we have missed.
    # Sum each element up to retrieve the number of all possible paths ("infinite" many miscleavages)
    return var_paths[sorted_nodes[-1]]



def _num_of_possible_paths_feature_type(graph_entry, feature_type="Variant", or_count=min):
    """
    Get the Number of all possible simple Paths with cleavages for a Protein or Peptide.
    A dynamic programming approach is taken here. We can minimize this problem
    into subproblems. The goal is to find the number of paths from the start of
    a protein to the end of a protein.

    This can be divided into the number of paths from the start to a node in the protein,
    which sums the number of possible paths from its previous nodes. This can be continued up to
    the end node, yielding the number of possible (non-repetative) paths to the end of a protein.

    For counting the miscleavages, we simply use a list instead of an value and shift, if an edge is
    an cleaved one. NOTE: Due to this approach, the calculation can be memory heavy!

    This algorithm therefore needs to iterate over the graph, which can be done with the help of
    the topological sort. The runtime of the dynamic programming part should be O(n^2) TODO is this 100% correct?
    """
    if "qualifiers" not in graph_entry.es[0].attributes():
        # Case: we do not have a property to execute this statistic
        return None

    # First get topological sorting of the graph
    sorted_nodes = graph_entry.topological_sorting()

    # Create list with num of paths as LISTS
    var_paths = [[]] * graph_entry.vcount()

    # Initialize Path from the very first node! For convenience we set it to one (actually 0!)
    # The very first node should always be the __start__ of a protein
    first = sorted_nodes[0]
    var_paths[first] = [1]  # as LIST entry!

    # Iterative approach look how many paths are possible from previous to itself (O(n^2))
    # TODO is the runtime 100% correct?
    # Get next node in topological sorted nodes
    for v in sorted_nodes[1:]:
        # Get the sum of possible paths from it, looking at each node which points to it
        summed = []
        for e_in in graph_entry.vs[v].in_edges():
            # This is the only change in the algorithm
            fts = e_in["qualifiers"]
            if fts is None:
                # simply add lists
                summed = _add_lists(summed, var_paths[e_in.source])
            else:
                offset = _resolve_or(fts, feature_type, or_count)
                # if cleaved then shift by one!
                summed = _add_lists(summed, [0]*offset + var_paths[e_in.source])# TODO

        # Set its number of paths as the calculated sum
        var_paths[v] = summed

    # We changed the first nodes value to 1, which is not correct, we set it back here!
    var_paths[first] = [0]  # Path to itself is zero! (lastly as LIST entry)

    # Returning the last node (which should always be the __end__ of a protein)
    # Here the index of the list gives us the number of how many cleavages we have missed.
    # Sum each element up to retrieve the number of all possible paths ("infinite" many miscleavages)
    return var_paths[sorted_nodes[-1]]



def _resolve_or(fts, feature_type, or_count):
    count = 0
    for ft in fts:
        if isinstance(ft, Or):
            t = [_resolve_or(or_ft, feature_type, or_count) for or_ft in ft]
            count += or_count(t)
        else:
            if ft.type == feature_type:
                count += 1

    return count


def _get_weight(seq, dict, mono=True):
    idx = 0 if mono else 1
    return sum([dict[x][idx] for x in seq.replace("__start__", "").replace("__end__", "")])


def _add_lists(list_a, list_b):
    """ Appends 0 if list A or B is too short. Does the operation '+' to two lists (vectors, etc.) """
    if len(list_a) > len(list_b):
        # Swap elements if the other is larger
        t = list_a
        list_a = list_b
        list_b = t
    return list(map(
        add,
        list_a + [0]*(len(list_b) - len(list_a)),
        list_b
    ))

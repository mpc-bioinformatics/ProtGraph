from operator import add

from protgraph.graph_collapse_edges import Or

# A list of statistic methods with csv keys, which should be calculated.
STATISTICS_METHOD_LIST = [
    ("calc_num_possibilities", lambda graph, or_count: _count_pos_paths(graph), "num_paths"),  # noqa: E501
    ("calc_num_possibilities_miscleavages", lambda graph, or_count: _count_miscleavages_list(graph), "list_paths_miscleavages"),  # noqa: E501
    ("calc_num_possibilities_hops", lambda graph, or_count: _count_hops_list(graph), "list_paths_hops"),  # noqa: E501
    ("calc_num_possibilities_variant", lambda graph, or_count: _count_feature_list(graph, feature_type="VARIANT", or_count=or_count), "list_paths_variant"),  # noqa: E501
    ("calc_num_possibilities_mutagen", lambda graph, or_count: _count_feature_list(graph, feature_type="MUTAGEN", or_count=or_count), "list_paths_mutagen"),  # noqa: E501
    ("calc_num_possibilities_conflict", lambda graph, or_count: _count_feature_list(graph, feature_type="CONFLICT", or_count=or_count), "list_paths_conflict"),  # noqa: E501
]


def get_statistics(graph, entry_dict, **kwargs):
    """
    TODO can we retrieve even more information!?
    returns #Node, #Edges, #Num_Of_Paths
    """

    # Get the number of nodes and edges (an be done instantly)
    entry_dict["num_edges"] = _get_edge_count(graph)
    entry_dict["num_nodes"] = _get_node_count(graph)

    for calculate_bool, method, entry_dict_key in STATISTICS_METHOD_LIST:
        if kwargs[calculate_bool]:
            entry_dict[entry_dict_key] = method(graph, kwargs["calc_num_possibilites_or_count"])

    # TODO can we calculate more statistics?


def _get_edge_count(graph_entry):
    """ Get number of edges"""
    return graph_entry.ecount()


def _get_node_count(graph_entry):
    """ Get number of nodes"""
    return graph_entry.vcount()


def _dynamic_programming(graph_entry, kernel_func, init_first_val=[1], _summed=[]):
    """
    Get the Number of all possible simple Paths depending on kernel function for a Protein or Peptide.
    A dynamic programming approach is taken here. We can minimize this problem
    into subproblems. The goal is to find the number of paths from the start of
    a protein to the end of a protein.

    E.G.: This can be divided into the number of paths from the start to a node in the protein,
    which sums the number of possible paths from its previous nodes. This can be continued up to
    the end node, yielding the number of possible (non-repetative) paths to the end of a protein.

    For counting the miscleavages, we simply use a list instead of an value and shift, if an edge is
    an cleaved one. NOTE: Due to such a approach, the calculation can be memory heavy!

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
    var_paths[first] = init_first_val  # as LIST entry!

    # Iterative approach look how many paths are possible from previous to itself (O(n^2))
    # TODO is the runtime 100% correct?
    # Get next node in topological sorted nodes
    for v in sorted_nodes[1:]:
        # Get the summed number from it, looking at each node which points to it
        summed = _summed
        for e_in in graph_entry.vs[v].in_edges():
            summed = kernel_func(e_in, summed, var_paths[e_in.source])

        # Set its number of paths as the calculated sum
        var_paths[v] = summed

    # We changed the first nodes value to 1, which is not correct, we set it back here!
    var_paths[first] = [0]  # Path to itself is zero! (lastly as LIST entry)

    # Returning the last node (which should always be the __end__ of a protein)
    return var_paths[sorted_nodes[-1]]


def _count_feature_list(graph_entry, feature_type="VARIANT", or_count=min):
    """ Wrapper for counting features like VARIANT|MUTAGEN|CONFLICT """
    if "qualifiers" not in graph_entry.es[0].attributes():
        # Case: we do not have a property to execute this statistic
        return None

    def kernel(edge, s_list, edge_list):
        offset = _count_feature(edge["qualifiers"], feature_type, or_count)
        # if feature then shift by offset!
        return _add_lists(s_list, [0]*offset + edge_list)

    return _dynamic_programming(graph_entry, kernel)


def _count_miscleavages_list(graph_entry):
    """ Wrapper for counting Miscleavages """
    if "cleaved" not in graph_entry.es[0].attributes():
        # Case: we do not have a property to execute this statistic
        return None

    def kernel(edge, s_list, edge_list):
        if edge["cleaved"]:  # if cleaved then shift by one!
            return _add_lists(s_list, [0] + edge_list)
        else:  # simply add lists
            return _add_lists(s_list, edge_list)

    return _dynamic_programming(graph_entry, kernel)


def _count_hops_list(graph_entry):
    """ Wrapper for counting hops (nodes) """
    # Special Case at the end, we need to remove first element
    return _dynamic_programming(graph_entry, lambda _, a, b: _add_lists(a, [0] + b))[1:]


def _count_pos_paths(graph_entry):
    """ Wrapper for counting number of possible paths """
    # Here we return directly the value, since list will always have len==1
    return _dynamic_programming(graph_entry, lambda _, a, b: _add_lists(a, b))[0]


def _extend_lists(graph_entry):
    """ DL EXREMELY BIG TODO, for a possible graph splitting algorithm! """
    def kernel(edge, a, b):
        return (
            a[0] + b[0],
            [*a[1], b[0]]
        )
    return _dynamic_programming(graph_entry, kernel, init_first_val=(1, [1]), _summed=(0, []))[0]


def _count_feature(fts, feature_type, or_count):
    count = 0
    if fts is None:
        return 0
    for ft in fts:
        if isinstance(ft, Or):
            t = [_count_feature(or_ft, feature_type, or_count) for or_ft in ft]
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

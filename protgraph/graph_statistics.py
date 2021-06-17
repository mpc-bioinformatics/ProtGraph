from operator import add

from collections import defaultdict
from protgraph.aa_masses_annotation import _get_mass_dict

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

    # TODO DL this is unfinished and may not work for some specific proteins!
    if kwargs["calc_possible_weigths"]:
        mass_dict = _get_mass_dict(factor=kwargs["mass_dict_factor"], type=kwargs["mass_dict_type"])
        set_of_weights = _possible_weights(graph, mass_dict)
    else: 
        set_of_weights = None

    # TODO can we calculate more statistics?

    # Return all information
    return num_nodes, num_edges, num_possible_paths, \
        num_possible_paths_all_mis, num_possible_paths_all_hops, set_of_weights


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



def _possible_weights(graph_entry, mass_dict):
    """
    Get the set of all possible weights inside a graph
    """

    from tqdm import tqdm
    # First get the reverse topological sorting of the graph
    sorted_nodes = graph_entry.topological_sorting(mode="IN")

    # Debug create your own top sort!
    # sorted_min = []
    # s = set([x for x, y in zip(range(graph_entry.vcount()), graph_entry.vs.outdegree()) if y == 0])
    # marked_edges = [False]*graph_entry.ecount()

    # while len(s) != 0:
    #     t = [(graph_entry.vs[x].outdegree(), x) for x in s]
    #     n = min(t, key=lambda x: x[0])[1]
    #     s.remove(n)

    #     sorted_min.append(n)
    #     for e_in in graph_entry.vs[n].in_edges():
    #         marked_edges[e_in.index] = True
    #         in_out_edges = [x.index for x in graph_entry.vs[e_in.source].out_edges()]
    #         if all([marked_edges[x] for x in in_out_edges]):
    #             s.add(e_in.source)

    # # Debug create your own top sort!
    # sorted_max = []
    # s = set([x for x, y in zip(range(graph_entry.vcount()), graph_entry.vs.outdegree()) if y == 0])
    # marked_edges = [False]*graph_entry.ecount()

    # while len(s) != 0:
    #     t = [(graph_entry.vs[x].outdegree(), x) for x in s]
    #     n = max(t, key=lambda x: x[0])[1]
    #     s.remove(n)

    #     sorted_max.append(n)
    #     for e_in in graph_entry.vs[n].in_edges():
    #         marked_edges[e_in.index] = True
    #         in_out_edges = [x.index for x in graph_entry.vs[e_in.source].out_edges()]
    #         if all([marked_edges[x] for x in in_out_edges]):
    #             s.add(e_in.source)


    # # Debug create your own top sort!
    # sorted_high = []
    # s = set([x for x, y in zip(range(graph_entry.vcount()), graph_entry.vs.outdegree()) if y == 0])
    # marked_edges = [False]*graph_entry.ecount()

    # while len(s) != 0:
    #     t = [(graph_entry.vs[x].outdegree(), x) for x in s]
    #     n = max(t, key=lambda x: x[1])[1]
    #     s.remove(n)

    #     sorted_high.append(n)
    #     for e_in in graph_entry.vs[n].in_edges():
    #         marked_edges[e_in.index] = True
    #         in_out_edges = [x.index for x in graph_entry.vs[e_in.source].out_edges()]
    #         if all([marked_edges[x] for x in in_out_edges]):
    #             s.add(e_in.source)


    # # Debug create your own top sort!
    # sorted_low = []
    # s = set([x for x, y in zip(range(graph_entry.vcount()), graph_entry.vs.outdegree()) if y == 0])
    # marked_edges = [False]*graph_entry.ecount()

    # while len(s) != 0:
    #     t = [(graph_entry.vs[x].outdegree(), x) for x in s]
    #     n = min(t, key=lambda x: x[1])[1]
    #     s.remove(n)

    #     sorted_low.append(n)
    #     for e_in in graph_entry.vs[n].in_edges():
    #         marked_edges[e_in.index] = True
    #         in_out_edges = [x.index for x in graph_entry.vs[e_in.source].out_edges()]
    #         if all([marked_edges[x] for x in in_out_edges]):
    #             s.add(e_in.source)

    # # Debug create your own top sort!
    # sorted_most_marked = []
    # s = set([x for x, y in zip(range(graph_entry.vcount()), graph_entry.vs.outdegree()) if y == 0])
    # marked_edges = [False]*graph_entry.ecount()

    # while len(s) != 0:
    #     t_ids = [x for x in s]
    #     t = []
    #     for t_temp in t_ids:
    #         count_t = 0
    #         for e_in in graph_entry.vs[t_temp].in_edges():
    #             in_out_edges = [x.index for x in graph_entry.vs[e_in.source].out_edges()]
    #             count_t += sum([marked_edges[x] for x in in_out_edges])
    #         t.append((count_t, t_temp))

    #     n = min(t, key=lambda x: x[0])[1]
    #     s.remove(n)

    #     sorted_most_marked.append(n)
    #     for e_in in graph_entry.vs[n].in_edges():
    #         marked_edges[e_in.index] = True
    #         in_out_edges = [x.index for x in graph_entry.vs[e_in.source].out_edges()]
    #         if all([marked_edges[x] for x in in_out_edges]):
    #             s.add(e_in.source)

    # # Debug create your own top sort!
    # sorted_least_marked = []
    # s = set([x for x, y in zip(range(graph_entry.vcount()), graph_entry.vs.outdegree()) if y == 0])
    # marked_edges = [False]*graph_entry.ecount()

    # while len(s) != 0:
    #     t_ids = [x for x in s]
    #     t = []
    #     for t_temp in t_ids:
    #         count_t = 0
    #         for e_in in graph_entry.vs[t_temp].in_edges():
    #             in_out_edges = [x.index for x in graph_entry.vs[e_in.source].out_edges()]
    #             count_t += sum([marked_edges[x] for x in in_out_edges])
    #         t.append((count_t, t_temp))

    #     n = max(t, key=lambda x: x[0])[1]
    #     s.remove(n)

    #     sorted_least_marked.append(n)
    #     for e_in in graph_entry.vs[n].in_edges():
    #         marked_edges[e_in.index] = True
    #         in_out_edges = [x.index for x in graph_entry.vs[e_in.source].out_edges()]
    #         if all([marked_edges[x] for x in in_out_edges]):
    #             s.add(e_in.source)




    # ####DEBUG
    # iteration_info_n = []
    # true_list = [False]*len(sorted_nodes)

    # true_list[sorted_nodes[0]] = True
    # iteration_info_n.append(sum(true_list))
    # for v in sorted_nodes:
    #     true_list[v] = False
    #     for v_prev in graph_entry.vs[v].in_edges():
    #         true_list[v_prev.source] = True
    #     iteration_info_n.append(sum(true_list))
    # ####DEBUG
    # iteration_info_min = []
    # true_list = [False]*len(sorted_min)

    # true_list[sorted_min[0]] = True
    # iteration_info_min.append(sum(true_list))
    # for v in sorted_min:
    #     true_list[v] = False
    #     for v_prev in graph_entry.vs[v].in_edges():
    #         true_list[v_prev.source] = True
    #     iteration_info_min.append(sum(true_list))
    # ####DEBUG
    # iteration_info_max = []
    # true_list = [False]*len(sorted_max)

    # true_list[sorted_max[0]] = True
    # iteration_info_max.append(sum(true_list))
    # for v in sorted_max:
    #     true_list[v] = False
    #     for v_prev in graph_entry.vs[v].in_edges():
    #         true_list[v_prev.source] = True
    #     iteration_info_max.append(sum(true_list))
    # ####DEBUG
    # iteration_info_high = []
    # true_list = [False]*len(sorted_high)

    # true_list[sorted_high[0]] = True
    # iteration_info_high.append(sum(true_list))
    # for v in sorted_high:
    #     true_list[v] = False
    #     for v_prev in graph_entry.vs[v].in_edges():
    #         true_list[v_prev.source] = True
    #     iteration_info_high.append(sum(true_list))
    # ####DEBUG
    # iteration_info_low = []
    # true_list = [False]*len(sorted_low)

    # true_list[sorted_low[0]] = True
    # iteration_info_low.append(sum(true_list))
    # for v in sorted_low:
    #     true_list[v] = False
    #     for v_prev in graph_entry.vs[v].in_edges():
    #         true_list[v_prev.source] = True
    #     iteration_info_low.append(sum(true_list))
    # ####DEBUG
    # iteration_info_most_m = []
    # true_list = [False]*len(sorted_most_marked)

    # true_list[sorted_most_marked[0]] = True
    # iteration_info_most_m.append(sum(true_list))
    # for v in sorted_most_marked:
    #     true_list[v] = False
    #     for v_prev in graph_entry.vs[v].in_edges():
    #         true_list[v_prev.source] = True
    #     iteration_info_most_m.append(sum(true_list))
    # ####DEBUG
    # iteration_info_least_m = []
    # true_list = [False]*len(sorted_least_marked)

    # true_list[sorted_least_marked[0]] = True
    # iteration_info_least_m.append(sum(true_list))
    # for v in sorted_least_marked:
    #     true_list[v] = False
    #     for v_prev in graph_entry.vs[v].in_edges():
    #         true_list[v_prev.source] = True
    #     iteration_info_least_m.append(sum(true_list))
    # ###

    # import matplotlib.pyplot as plt

    # plt.plot(iteration_info_n, label="normal_igraph")
    # plt.plot(iteration_info_min, label="min_edges")
    # plt.plot(iteration_info_max, label="max_edges")
    # plt.plot(iteration_info_low, label="low_id")
    # plt.plot(iteration_info_high, label="high_id")
    # plt.plot(iteration_info_most_m, label="most_edges_marked")
    # plt.plot(iteration_info_least_m, label="least_edges_marked")
    # plt.legend()
    # plt.show()

    # sorted_nodes = sorted_most_marked
    # Create list with weights
    var_weights = defaultdict(set)

    # Initialize Path from the very last node!
    last = sorted_nodes[0]
    var_weights[last] = set([0])
    weight = _get_weight(graph_entry.vs[last]["aminoacid"], mass_dict)
    for v_prev in graph_entry.neighbors(last, mode="IN"):
        var_weights[v_prev].update([x + weight for x in var_weights[last]])
    del var_weights[last]

    # Get next node in topological sorted nodes
    for v in tqdm(sorted_nodes[1:-1]):
        weight = _get_weight(graph_entry.vs[v]["aminoacid"], mass_dict)
        res_list = [x + weight for x in var_weights[v]]
        for v_prev in tqdm(graph_entry.neighbors(v, mode="IN"), leave=False):
            var_weights[v_prev].update(res_list)
        del var_weights[v]

    return var_weights[graph_entry.vs.select(aminoacid="__start__")[0].index]


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

def _get_mass_dict(factor=1000000000, type=int):
    """
    Return a Dictionary containing the masses of each aminoacid
    We explicitly convert them by a factor of 1 000 000 000 (default) into integers

    The values are taken from: https://proteomicsresource.washington.edu/protocols06/masses.php
    """
    return dict(  # In format: AA = (MONO_MASS, AVG_MASS)
        G=(type(57.021463735 * factor), type(57.05132 * factor)),
        A=(type(71.037113805 * factor), type(71.0779 * factor)),
        S=(type(87.032028435 * factor), type(87.0773 * factor)),
        P=(type(97.052763875 * factor), type(97.11518 * factor)),
        V=(type(99.068413945 * factor), type(99.13106 * factor)),
        T=(type(101.047678505 * factor), type(101.10388 * factor)),
        C=(type(103.009184505 * factor), type(103.1429 * factor)),
        L=(type(113.084064015 * factor), type(113.15764 * factor)),
        I=(type(113.084064015 * factor), type(113.15764 * factor)),
        N=(type(114.042927470 * factor), type(114.10264 * factor)),
        D=(type(115.026943065 * factor), type(115.0874 * factor)),
        Q=(type(128.058577540 * factor), type(128.12922 * factor)),
        K=(type(128.094963050 * factor), type(128.17228 * factor)),
        E=(type(129.042593135 * factor), type(129.11398 * factor)),
        M=(type(131.040484645 * factor), type(131.19606 * factor)),
        H=(type(137.058911875 * factor), type(137.13928 * factor)),
        F=(type(147.068413945 * factor), type(147.17386 * factor)),
        U=(type(150.953633405 * factor), type(150.3079 * factor)),
        R=(type(156.101111050 * factor), type(156.18568 * factor)),
        Y=(type(163.063328575 * factor), type(163.17326 * factor)),
        W=(type(186.079312980 * factor), type(186.2099 * factor)),
        O=(type(237.147726925 * factor), type(237.29816 * factor)),
        # Special Aminoacids
        J=(type(113.084064015 * factor), type(113.1594 * factor)),
        X=(type(0.0 * factor), type(0.0 * factor)),  # Unknown Amino Acid
        Z=(type(128.55059 * factor), type(128.6231 * factor)),
        B=(type(114.53495 * factor), type(114.5962 * factor)),
        # Custom start and end points
        __start__=(type(0), type(0)),
        __end__=(type(0), type(0)),
    )


def annotate_weights(graph_entry, **kwargs):
    """
    This method annotates the graph (if explicitly set) with the
    average and monoisotopic weight. It also has the possibility to set
    the end weight of the nodes. This value tells from a specific node
    how much weight it has at least to go to get to the end.

    Such information could e.g. be used to reduce the number of paths early on,
    when using a customized DFS or BFS


    kwargs arguments:
    1. annotate_mono_weights     : If True, then mono weigths are added
    2. annotate_mono_end_weights : If True, then 1. and the lowest weight to end is added

    3. annotate_avrg_weights     : If True, then mono weigths are added
    4. annotate_avrg_end_weights : If True, then 3. and the lowest weight to end is added

    NOTE: This transforms the graph without returning it!

    Following Keys maybe set here:
    Nodes: "mono_weight", "avrg_weight", "mono_weight_to_end", "avrg_weight_to_end"
    Edges: <None>
    """
    # Get the mass dictionary, which should be happen instantly.
    mass_dict = _get_mass_dict(factor=kwargs["mass_dict_factor"], type=kwargs["mass_dict_type"])

    # If mono or mono_end is set, annotate the graph with the mono weights
    if kwargs["annotate_mono_weights"] or kwargs["annotate_mono_weight_to_end"]:
        _add_masses(graph_entry, "mono_weight", mass_dict, kwargs["mass_dict_type"], 0)

    # If avrg or avrg_end is set, annotate the graph with the average weights
    if kwargs["annotate_avrg_weights"] or kwargs["annotate_avrg_weight_to_end"]:
        _add_masses(graph_entry, "avrg_weight", mass_dict, kwargs["mass_dict_type"], 1)

    # If one of the end weights is set:
    if kwargs["annotate_mono_weight_to_end"] or kwargs["annotate_avrg_weight_to_end"]:
        # Then get the (reverse) topological sort of the graph ("all pair shortest path manner")
        top_sorted_nodes = graph_entry.topological_sorting(mode="IN")

        # If mono_end is set then annotate it
        if kwargs["annotate_avrg_weight_to_end"]:
            _add_end_masses(graph_entry, "avrg_weight_to_end", "avrg_weight", top_sorted_nodes)

        # If avrg_end is set then annotate it
        if kwargs["annotate_mono_weight_to_end"]:
            _add_end_masses(graph_entry, "mono_weight_to_end", "mono_weight", top_sorted_nodes)


def _add_masses(graph_entry, weight_name, mass_dict, mass_dict_type, idx):
    """
    Here we simply iterate over each edge, get their target node (sum up, if it contains multiple aminoacid entries
    due to merging) and set the masses in the graph.

    weight_name: name of the weight
    mass_dict: The masses dictionary (as in this file)
    idx: Which entry of the list-value from the dictionary should be taken
    """
    mono_masses = [
        mass_dict_type(
            sum(
                [
                    mass_dict[y][idx]
                    for y in x["aminoacid"]
                    .replace("__start__", "")
                    .replace("__end__", "")
                    # 2. Get the aminoacid (-chain) from target and sum it up
                ]
            ) + (x["delta_mass"] if "delta_mass" in x.attributes() and x["delta_mass"] is not None else 0)
        )
        for x in graph_entry.vs[:]
        # 1. For each edge in the graph
    ]
    # Then set the masses for each edge
    graph_entry.vs[:][weight_name] = mono_masses
    if "delta_mass" in graph_entry.vs[0].attributes():
        del graph_entry.vs["delta_mass"]  # Delete deltamass if available


def _add_end_masses(graph_entry, weight_end_name, weight_name, sorted_nodes):
    """
    By using the sorted nodes (reverse topologically sorted), we can simply iterate over this list
    and do a second iteration over each edge (which targets this node). Those edges can only
    target nodes further in the graph. We simply check, wheather we can reach it with a lower
    weight and set it appropiately. Runtime should be: 0(n + e)  (n=Nodes, e=Edges).

    This can be done, since we are working with DAGs

    weigth_end_name: The name for the end weight
    weight_name: An ALREADY set weight in the graph
    sorted_nodes: Sorted nodes in REVERSE topological order.
    """
    # First set all weights to Infinity
    graph_entry.vs[:][weight_end_name] = float("inf")
    # The target Node has end weight as 0 (Initialization point)
    graph_entry.vs.select(aminoacid="__end__")[weight_end_name] = 0

    # Complexity O (n + e)
    # For each node in sorted nodes:
    for node in sorted_nodes:
        # Grab its edges, which targets the node. For each of the edge:
        for edge in graph_entry.es.select(_target=node):
            # Get the end weight of the node and the weight from the current edge
            temp_weight = edge.target_vertex[weight_end_name] + edge.source_vertex[weight_name]
            # Set this weight for the source node if smaller
            if graph_entry.vs[edge.source][weight_end_name] > temp_weight:
                graph_entry.vs[edge.source][weight_end_name] = temp_weight

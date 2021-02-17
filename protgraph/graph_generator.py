import igraph
from protgraph.aa_masses_annotation import annotate_weights
from protgraph.digestion import digest
from protgraph.export.exporters import Exporters
from protgraph.ft_execution.init_met import execute_init_met
from protgraph.ft_execution.signal import execute_signal
from protgraph.ft_execution.var_seq import _get_isoforms_of_entry, execute_var_seq
from protgraph.ft_execution.variant import execute_variant
from protgraph.graph_statistics import get_statistics
from protgraph.merge_aminoacids import merge_aminoacids
from protgraph.verify_graphs import verify_graph


def _generate_canonical_graph(sequence: str, acc: str):
    """
    Generates the canonical directed graph from a sequence.
    This simply generates a chain of nodes and edges and sets
    specific attributes for them.
    """
    # Initialization of the directed graph (DiGraph)
    graph = igraph.Graph(directed=True)

    # Initialize the graph with the length of the sequence
    graph.add_vertices(len(sequence) + 2)  # +2 -> adding start and end node here!
    graph.add_edges([(x1, x1 + 1) for x1 in range(len(sequence) + 1)])

    # Add their amino acid to the corresponding nodes
    graph.vs["aminoacid"] = ["__start__", *[x for x in sequence], "__end__"]

    # Add position attributes to nodes as well as from which accesion they originate
    graph.vs["position"] = list(range(len(sequence) + 2))  # Set position of aa on every node!
    graph.vs["accession"] = [acc, *[acc] * len(sequence), acc]  # Set accession on every node!

    return graph


def _sort_entry_features(entry):
    """ This sorts the features according to their type into a dict. """
    sorted_features = dict()
    # For each features
    for f in entry.features:
        # Append it to a list to its corresponding key -> type
        if f.type not in sorted_features:
            sorted_features[f.type] = [f]
        else:
            sorted_features[f.type].append(f)

    # Return the dictionary
    return sorted_features


def _include_ft_information(entry, graph, kwargs):
    """ Returns num of possible isoforms and others (on the fly) """
    # Sort features of entry according to their type into a dict
    sorted_features = _sort_entry_features(entry)

    # VAR_SEQ (isoforms) need to be executed at once and before all other variations
    # since those can be referenced by others
    num_of_isoforms = 0 if not kwargs["skip_isoforms"] else None
    if "VAR_SEQ" in sorted_features and not kwargs["skip_isoforms"]:
        # Get isoform information of entry as a dict
        isoforms, num_of_isoforms = _get_isoforms_of_entry(entry.comments, entry.accessions[0])
        execute_var_seq(isoforms, graph, entry.sequence, sorted_features["VAR_SEQ"], entry.accessions[0])

    num_of_init_m = 0 if not kwargs["skip_init_met"] else None
    if "INIT_MET" in sorted_features and not kwargs["skip_init_met"]:
        num_of_init_m = len(sorted_features["INIT_MET"])
        for f in sorted_features["INIT_MET"]:
            execute_init_met(graph, f)

    num_of_signal = 0 if not kwargs["skip_signal"] else None
    if "SIGNAL" in sorted_features and not kwargs["skip_signal"]:
        num_of_signal = len(sorted_features["SIGNAL"])
        for f in sorted_features["SIGNAL"]:
            execute_signal(graph, f)

    num_of_variant = 0 if not kwargs["skip_variants"] else None
    if "VARIANT" in sorted_features and not kwargs["skip_variants"]:
        num_of_variant = len(sorted_features["VARIANT"])
        for f in sorted_features["VARIANT"]:
            execute_variant(graph, f)

    return num_of_isoforms, num_of_init_m, num_of_signal, num_of_variant


def generate_graph_consumer(entry_queue, graph_queue, **kwargs):
    """
    TODO
    describe kwargs and consumer until a graph is generated and digested etc ...
    """
    # Initialize the exporters for graphs
    graph_exporters = Exporters(**kwargs)

    while True:
        # Get next entry
        entry = entry_queue.get()

        # Stop if entry is None
        if entry is None:
            # --> Stop Condition of Process
            break

        # Beginning of Graph-Generation
        # We also collect interesting information here!

        # Generate canonical graph (initialization of the graph)
        graph = _generate_canonical_graph(entry.sequence, entry.accessions[0])

        # FT parsing and appending of Nodes and Edges into the graph
        # The amount of isoforms, etc.. can be retrieved on the fly
        num_isoforms, num_initm, num_signal, num_variant = _include_ft_information(entry, graph, kwargs)

        # Digest graph with enzyme (unlimited miscleavages)
        num_of_cleavages = digest(graph, kwargs["digestion"])

        # Merge (summarize) graph if wanted
        if not kwargs["no_merge"]:
            merge_aminoacids(graph)

        # Annotate weights for edges and nodes (maybe even the smallest weight possible to get to the end node)
        annotate_weights(graph, **kwargs)

        # Calculate statistics on the graph:
        num_nodes, num_edges, num_paths, num_paths_miscleavages, num_paths_hops = get_statistics(graph, **kwargs)

        # Verify graphs if wanted:
        if kwargs["verify_graph"]:
            verify_graph(graph)

        # Persist or export graphs with speicified exporters
        graph_exporters.export_graph(graph)

        # Output statistics we gathered during processing
        entry_protein_desc = entry.description.split(";", 1)[0]
        entry_protein_desc = entry_protein_desc[entry_protein_desc.index("=") + 1:]
        graph_queue.put(
            (
                entry.accessions[0],  # Protein Accesion
                entry.entry_name,  # Protein displayed name
                num_isoforms,  # Number of Isoforms
                num_initm,  # Number of Init_M (either 0 or 1)
                num_signal,  # Number of Signal Peptides used (either 0 or 1)
                num_variant,  # Number of Variants applied to this protein
                num_of_cleavages,  # Number of cleavages (marked edges) this protein has
                num_nodes,  # Number of nodes for the Protein/Peptide Graph
                num_edges,  # Number of edges for the Protein/Peptide Graph
                num_paths,  # Possible (non repeating paths) to the end of a graph. (may conatin repeating peptides)
                num_paths_miscleavages,  # As num_paths, but binned to the number of miscleavages (by list idx, at 0)
                num_paths_hops,  # As num_paths, only that we bin by hops (E.G. useful for determine DFS or BFS depths)
                entry_protein_desc,  # Description name of the Protein (can be lenghty)
            )
        )

    # Close exporters (maybe opened files, database connections, etc... )
    graph_exporters.close()

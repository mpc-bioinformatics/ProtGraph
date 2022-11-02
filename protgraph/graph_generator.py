import io
from collections import defaultdict

import igraph
from Bio import SwissProt

from protgraph.aa_masses_annotation import annotate_weights
from protgraph.aa_replacer import replace_aa
from protgraph.annotate_ptms import annotate_ptms
from protgraph.digestion import digest
from protgraph.export.exporters import Exporters
from protgraph.ft_execution.generic import (execute_conflict, execute_mutagen,
                                            execute_variant)
from protgraph.ft_execution.generic_cleaved_peptide import (execute_chain,
                                                            execute_peptide,
                                                            execute_propeptide)
from protgraph.ft_execution.init_met import execute_init_met
from protgraph.ft_execution.signal import execute_signal
from protgraph.ft_execution.var_seq import (_get_isoforms_of_entry,
                                            execute_var_seq)
from protgraph.graph_collapse_edges import collapse_parallel_edges
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

    # Add position attributes to nodes as well as from which accession they originate
    graph.vs["position"] = list(range(len(sequence) + 2))  # Set position of aa on every node!
    graph.vs["accession"] = [acc, *[acc] * len(sequence), acc]  # Set accession on every node!

    return graph


def _sort_entry_features(entry):
    """ This sorts the features according to their type into a dict. """
    sorted_features = defaultdict(list)
    # For each features
    for f in entry.features:
        # Append it to a list to its corresponding key -> type
        sorted_features[f.type].append(f)

    # Return the dictionary
    return sorted_features


def _include_spefic_ft(graph, ft_type, method, sorted_features, ft_dict):
    """ Execute features individually """
    num_of_feature_type = 0 if ft_type in ft_dict else None
    if ft_type in sorted_features and ft_type in ft_dict:
        num_of_feature_type = len(sorted_features[ft_type])
        for f in sorted_features[ft_type]:
            method(graph, f)
    return num_of_feature_type


# Features which are applied in order on the graph
FT_SINGLE_EXECUTION = [
    ("INIT_MET", execute_init_met, "num_init_met"),
    ("SIGNAL", execute_signal, "num_signal"),
    ("PROPEP", execute_propeptide, "num_propep"),
    ("PEPTIDE", execute_peptide, "num_peptide"),
    ("CHAIN", execute_chain, "num_chain"),
    ("VARIANT", execute_variant, "num_variant"),
    ("MUTAGEN", execute_mutagen, "num_mutagen"),
    ("CONFLICT", execute_conflict, "num_conflict"),
]


def _include_ft_information(entry, graph, ft_dict, entry_dict):
    """
    Apply feature (first isoform, then others as specified in order) on the graph.

    This method can count the number of applied features on the graph on the fly
    """
    # Sort features of entry according to their type into a dict
    sorted_features = _sort_entry_features(entry)

    # VAR_SEQ (isoforms) need to be executed at once and before all other variations
    # since those can be referenced by others
    if "VAR_SEQ" in ft_dict:
        entry_dict["num_var_seq"] = 0
    if "VAR_SEQ" in sorted_features and "VAR_SEQ" in ft_dict:
        # Get isoform information of entry as a dict
        isoforms, num_of_isoforms = _get_isoforms_of_entry(entry.comments, entry.accessions[0])
        entry_dict["num_var_seq"] = num_of_isoforms
        execute_var_seq(isoforms, graph, entry.sequence, sorted_features["VAR_SEQ"], entry.accessions[0])

    # Execute the other features tables, per feature
    for (feature, method, entry_dict_key) in FT_SINGLE_EXECUTION:
        entry_dict[entry_dict_key] = _include_spefic_ft(graph, feature, method, sorted_features, ft_dict)


def generate_graph_consumer(entry_queue, graph_queue, common_out_queue, proc_id, **kwargs):
    """
    TODO
    describe kwargs and consumer until a graph is generated and digested etc ...
    """
    # Set proc id
    kwargs["proc_id"] = proc_id

    # Set feature_table dict boolean table
    ft_dict = dict()
    if kwargs["feature_table"] is None or len(kwargs["feature_table"]) == 0 or "ALL" in kwargs["feature_table"]:
        ft_dict = dict(
            PEPTIDE=True, PROPEP=True, VARIANT=True, CHAIN=True,
            VAR_SEQ=True, SIGNAL=True, INIT_MET=True, MUTAGEN=True, CONFLICT=True
        )
    else:
        for i in kwargs["feature_table"]:
            ft_dict[i] = True

    # Initialize the exporters for graphs
    with Exporters(**kwargs) as graph_exporters:

        while True:
            # Get next entry
            io_entry = entry_queue.get()

            # Stop if entry is None
            if io_entry is None:
                # --> Stop Condition of Process
                break

            # Parse entry in graph generation process, so that more work is done in the consumer
            entry = SwissProt.read(io.BytesIO(io_entry))

            if kwargs["exclude_accessions"] and entry.accessions[0] in kwargs["exclude_accessions"]:
                # This effectively skips an entry at the cost to check whether to skip in EACH entry!
                continue

            # create new dict to collect information about this entry
            entry_dict = dict()

            # Beginning of Graph-Generation
            # We also collect interesting information here!

            # Generate canonical graph (initialization of the graph)
            graph = _generate_canonical_graph(entry.sequence, entry.accessions[0])

            # FT parsing and appending of nodes and edges into the graph
            # The amount of isoforms, etc.. can be retrieved on the fly
            _include_ft_information(entry, graph, ft_dict, entry_dict)

            # Replace Amino Acids based on user defined rules: E.G.: "X -> A,B,C"
            replace_aa(graph, kwargs["replace_aa"])

            # Digest graph with enzyme (unlimited miscleavages)
            digest(graph, kwargs["digestion"], entry_dict)

            # Annotate delta masses for PTMs
            annotate_ptms(graph, kwargs["variable_mod"], kwargs["fixed_mod"], kwargs["mass_dict_factor"])

            # Collapse parallel edges in a graph
            if not kwargs["no_collapsing_edges"]:
                collapse_parallel_edges(graph)

            # Merge (summarize) graph if wanted
            if not kwargs["no_merge"]:
                merge_aminoacids(graph)

            # Annotate weights for edges and nodes (maybe even the smallest weight possible to get to the end node)
            annotate_weights(graph, **kwargs)

            # Calculate statistics on the graph:
            get_statistics(graph, entry_dict, **kwargs)

            # Verify graphs if wanted:
            if kwargs["verify_graph"]:
                verify_graph(graph)

            # Persist or export graphs with speicified exporters
            graph_exporters.export_graph(graph, common_out_queue)

            # Output statistics we gathered during processing
            if kwargs["no_description"]:
                entry_dict["protein_description"] = entry_protein_desc = None
            else:
                entry_protein_desc = entry.description.split(";", 1)[0]
                entry_dict["protein_description"] = entry_protein_desc[entry_protein_desc.index("=") + 1:]

            # Set accession and gene for csv
            entry_dict["accession"] = entry.accessions[0]
            entry_dict["gene_id"] = entry.entry_name

            # Return the statistics which were retrieved
            graph_queue.put(
                [entry_dict[x] if x in entry_dict else None for x in kwargs["output_csv_layout"]]
            )

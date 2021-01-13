import igraph
from Bio.SeqFeature import UncertainPosition, UnknownPosition

from aa_masses_annotation import annotate_weights
# TODO comment and beautify this file!!!
from digestion import digest
from export.exporters import Exporters
from graph_statistics import get_statistics
from merge_aminoacids import merge_aminoacids

from ft_execution.var_seq import _get_isoforms_of_entry, execute_var_seq


def _execute_init_met(graph, init_met_feature):
    """ DONE? NOTE This may contain flaws! """

    # Get start
    [start] = graph.vs.select(aminoacid="__start__")

    # Get all neigbors of start
    pos_first_aas = graph.neighbors(start, mode="out")

    # Filtering is important, since we may skip an aminoacid of an already cleaved protein
    met_aas = [
        x 
        for x in pos_first_aas 
        if graph.vs[x]["aminoacid"] == "M" and graph.vs[x]["position"] == 1
    ]

    # Get possible features 
    features = [
        graph.es.select(_source=start, _target=x)["qualifiers"] for x in met_aas
    ]
    # Get the next nodes which should be after them
    targets = [graph.neighbors(x, mode="out") for x in met_aas]

    # Get current number of edges
    cur_count = graph.ecount()

    # Generate the edges list
    edge_list = [(start, y) for x in targets for y in x]
    qualifiers_list = [
        [init_met_feature, *features[x_idx][y_idx]]
        for x_idx, x in enumerate(targets)
        for y_idx, _ in enumerate(x)
    ]

    # Add new edges (bulk)
    graph.add_edges(edge_list)
    graph.es[cur_count:]["qualifiers"] = qualifiers_list 
    



def _generate_canonical_graph(sequence: str, acc: str):
    '''
        Generates the canonical directed graph from a sequence.
        This simply generates a chain of nodes and edges and sets
        specific attributes for them.
    '''
    # Initialization of the directed graph (DiGraph)
    graph = igraph.Graph(directed=True)

    # Initialize the graph with the length of the sequence
    graph.add_vertices(len(sequence) + 2)  # +2 -> adding start and end node here!
    graph.add_edges([(x1, x1+1) for x1 in range(len(sequence) + 1)])

    # Add their amino acid to the corresponding nodes
    graph.vs["aminoacid"] = ["__start__", *[x for x in sequence], "__end__"]

    # Add position attributes to nodes as well as from which accesion they originate
    graph.vs["position"] = list(range(len(sequence) + 2))  # Set position of aa on every node!
    graph.vs["accession"] = [acc, *[acc]*len(sequence), acc]  # Set accession on every node!

    # TODO add somewhere else, only if needed!
    # graph.vs["isoform_accession"] = [None, *[None]*len(sequence), None]
    # graph.vs["isoform_position"] = [None, *[None]*len(sequence), None] # Position acording to the isoform
    # Add Information about the Features used
    # graph.es["qualifiers"] = [[] for _ in range(len(sequence) + 2)]
    # graph.es["cleaved"] = False # And each added edge is currently not cleaved!

    return graph



def _combine_vertices(list_a, list_b):
    out_d = {}

    for a in list_a:
        if a["isoform_accession"] not in out_d:
            out_d[a["isoform_accession"]] = dict(inn=[a])
        else: 
            # TODO we cannot simply associate via accession!?!?
            raise Exception("Key is multiple times in dict") 
    for b in list_b:
        if b["isoform_accession"] not in out_d:
            out_d[b["isoform_accession"]] = dict(out=[b])
        else:
            if "out" not in out_d[b["isoform_accession"]]:
                out_d[b["isoform_accession"]]["out"] = [b]
            else:
                out_d[b["isoform_accession"]]["out"].append(b)
            if len(out_d[b["isoform_accession"]]) > 2:
                # TODO we cannot simply associate via accession!?!?
                raise Exception("multiple Entries Found")
    # TODO Does these two cases ever happen? 

    k = list(out_d.items())
    a_s = [x[1]["inn"] if "inn" in x[1] else [] for x in k]
    b_s = [x[1]["out"] if "out" in x[1] else [] for x in k]

    return a_s, b_s


def _execute_variant(graph, variant_feature):
    # Check if missing or if we add another node
    text = variant_feature.qualifiers["note"] 

    # TODO what about nodes which are sequentielly replaced
    # add missing node here use select for all available previous nodes (using acc AND position)!!!

    if text.lower().startswith("missing"):
        # A sequence is missing! Just append an edge and its information
        # NOTE: Shifted by 1 due to the __start__ node at 0
        aa_before = variant_feature.location.start + 1
        aa_after = variant_feature.location.end + 0

        if variant_feature.ref is None:
            vertices_before_raw = list(graph.vs.select(position=aa_before))
            vertices_after_raw = list(graph.vs.select(position=aa_after))
        else: 
            vertices_before_raw = list(
                graph.vs.select(isoform_position=aa_before, isoform_accession=variant_feature.ref)
            )
            vertices_after_raw = list(
                graph.vs.select(isoform_position=aa_after, isoform_accession=variant_feature.ref)
            )

        if len(vertices_before_raw) == 0 or len(vertices_after_raw) == 0:
            # TODO does this happen?
            raise Exception("TODO this should never happen!")

        vertices_before, vertices_after = _combine_vertices(vertices_before_raw, vertices_after_raw)

        # Association via combine Vertices OK ?? TODO
        edge_list = []
        for aa_in, aa_out in zip(vertices_before, vertices_after):
            for aa_in_in in aa_in:
                for aa_in_in_in in list(graph.es.select(_target=aa_in_in)) :
                    for aa_out_out in aa_out:
                        for aa_out_out_out in list(graph.es.select(_source=aa_out_out)) :
                            edge_list.append( ((aa_in_in_in.source, aa_out_out_out.target), [*aa_in_in_in["qualifiers"], variant_feature]) )

        # Bulk add of edges
        cur_edges = graph.ecount()
        graph.add_edges([x[0] for x in edge_list])
        graph.es[cur_edges:]["qualifiers"] = [x[1] for x in edge_list]

    else:
        # Get   X -> Y   Information
        idx = text.find("(")
        if idx != -1:
            text = text[:idx]
        xy = text.split("->")
        assert len(xy) == 2
        y = xy[1].strip().replace(" ", "")
        
        # Get start and end position
        # NOTE: Shifted by 1 due to the __start__ node at 0
        aa_before = variant_feature.location.start + 1
        aa_after = variant_feature.location.end + 0

        # Get list of all aa before and after
        if variant_feature.ref is None:
            vertices_before_raw = list(graph.vs.select(position=aa_before))
            vertices_after_raw = list(graph.vs.select(position=aa_after))
        else: 
            vertices_before_raw = list(
                graph.vs.select(isoform_position=aa_before, isoform_accession=variant_feature.ref)
            )
            vertices_after_raw = list(
                graph.vs.select(isoform_position=aa_after, isoform_accession=variant_feature.ref)
            )

        if len(vertices_before_raw) == 0 or len(vertices_after_raw) == 0:
            # TODO does this happen?
            raise Exception("TODO this should never happen!")

        vertices_before, vertices_after = _combine_vertices(
            vertices_before_raw, vertices_after_raw
        )

        # Association via combine Vertices OK ?? TODO
        for aa_in, aa_out in zip(vertices_before, vertices_after):
            if len(aa_out) == 0 or len(aa_in) == 0:
                # Skip this entry, since we do not have complete information
                # -> Not possible to link either start or end or both
                continue

            # Append new node in graph and mapping
            # Case for one or multiple amino acids TODO
            # Add each individual amino acid as a node
            y_idcs = []
            for entry in y:
                vertex = graph.add_vertex()
                graph.vs[vertex.index]["aminoacid"] = entry
                graph.vs[vertex.index]["accession"] = graph.vs[1]["accession"]
                # TODO add position??
                y_idcs.append(vertex.index)

            # Add edges between them:
            for idx, n in enumerate(y_idcs[:-1]):
                graph.add_edges([(n, y_idcs[idx+1])])

            # Get first node and last node
            first_node, last_node = y_idcs[0], y_idcs[-1]

            # Add edges to the original protein
            edge_list = []

            for aa_in_in in aa_in:
                for aa_in_in_in in list(graph.es.select(_target=aa_in_in)) :
                    edge_list.append( ((aa_in_in_in.source, first_node), [*aa_in_in_in["qualifiers"], variant_feature]) )
            for aa_out_out in aa_out:
                for aa_out_out_out in list(graph.es.select(_source=aa_out_out)) :
                    edge_list.append( ((last_node, aa_out_out_out.target), []) )

            # Bulk add of edges
            cur_edges = graph.ecount()
            graph.add_edges([x[0] for x in edge_list])
            graph.es[cur_edges:]["qualifiers"] = [x[1] for x in edge_list]


def _execute_signal(graph, signal_feature):

    if isinstance(signal_feature.location.end, UnknownPosition):
        # The Position of the end is not known. Therefore we skip
        # this entry simply 
        return


    # Get end node
    [__stop_node__]  = graph.vs.select(aminoacid="__end__")

    # get start and end position of signal peptide
    start_position, end_position = signal_feature.location.start + 1, signal_feature.location.end + 0



    # Get all nodes with position start_position
    all_start_points = list(graph.vs.select(position=start_position))

    all_end_points = [list(graph.vs.select(isoform_accession=x["isoform_accession"], position=end_position)) for x in all_start_points]
    # Check if only one end point exists for each start
    for x in all_end_points:
        if len(x) > 1:
            print("WARNING, there are multiple defined ENDPOINTS!!")


    # TODO should the signal peptide to be exactly the same as in canonical? Or can we leave it as is?
    # Should we check this? If so: do this here!!
    # TODO can we associate the end points with the start points, according to its index?

    # create information list
    all_edges = []
    for s, e in zip(all_start_points, all_end_points):
        for ee in e:
        
            edges_in = list(graph.es.select(_target=s))
            edges_out = list(graph.es.select(_source=ee))

            for ei in edges_in:
                for eo in edges_out:
                    all_edges.append(( (ei.source_vertex, eo.target_vertex), [*ei["qualifiers"], signal_feature] ) )

            
            # Special case the end can go directly to the end node
            all_edges.append(( (ee, __stop_node__), [signal_feature] ) )


    # Bulk adding edges
    cur_edges = graph.ecount()
    graph.add_edges([x[0] for x in all_edges])
    graph.es[cur_edges:]["qualifiers"] = [x[1] for x in all_edges]



def _sort_entry_features(entry):
    ''' This sorts the features according to their type into a dict. '''
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
    """ Returns num of possible isoforms and others (on the fly) TODO """ 
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
            _execute_init_met(graph, f)

    num_of_signal = 0 if not kwargs["skip_signal"] else None
    if "SIGNAL" in sorted_features and not kwargs["skip_signal"]:
        num_of_signal = len(sorted_features["SIGNAL"])
        for f in sorted_features["SIGNAL"]:
            _execute_signal(graph, f)

    num_of_variant = 0 if not kwargs["skip_variants"] else None
    if "VARIANT" in sorted_features and not kwargs["skip_variants"]:
        num_of_variant = len(sorted_features["VARIANT"])
        for f in sorted_features["VARIANT"]:
            _execute_variant(graph, f)

    return num_of_isoforms, num_of_init_m, num_of_signal, num_of_variant


# TODO parse note Missing or replace information! via method!
def generate_graph_consumer(entry_queue, graph_queue, **kwargs):
    """ TODO
        
        describe kwargs and consumer until a graph is generated and digested etc ...  TODO
    """

    # Initialize the exporters for graphs
    graph_exporters = Exporters(**kwargs)


    while True:
        # Get next entry
        entry = entry_queue.get()

        # Stop if entry is None
        if entry == None:
            # --> Stop Condition of Process
            break

        ### Beginning of Graph-Generation
        ### We also collect interesting information here!

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
        num_nodes, num_edges, num_paths = get_statistics(graph, **kwargs)

        # Persist or export graphs with speicified exporters
        graph_exporters.export_graph(graph)

        # Output statistics we gathered during processing
        entry_protein_desc = entry.description.split(";", 1)[0]
        entry_protein_desc = entry_protein_desc[entry_protein_desc.index("=")+1:]
        graph_queue.put(
            (
                entry.accessions[0], # Protein Accesion
                entry.entry_name,    # Protein displayed name
                num_isoforms,        # Number of Isoforms
                num_initm,           # Number of Init_M (either 0 or 1)
                num_signal,          # Number of Signal Peptides used (either 0 or 1)
                num_variant,         # Number of Variants applied to this protein
                num_of_cleavages,    # Number of cleavages (marked edges) this protein has
                num_nodes,           # Number of nodes for the Protein/Peptide Graph
                num_edges,           # Number of edges for the Protein/Peptide Graph
                num_paths,           # Number of possible (non repeating paths) to the end of the Protein/Peptide. NOTE: This can contain repeating Peptides!
                entry_protein_desc   # Description name of the Protein (which can be lenghty and is therefore added to the end)
            )
        )

    # Close exporters (maybe opened files, database connections, etc... )
    graph_exporters.close() 

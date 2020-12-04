
import igraph
from Bio.SeqFeature import UncertainPosition, UnknownPosition


def _execute_init_met(graph, init_met_feature):
    ''' DONE? NOTE This may contain flaws! '''


    # Get start
    [start] = graph.vs.select(aminoacid="__start__")

    # Get all neigbors of start
    pos_first_aas = graph.neighbors(start, mode="out")

    # Filtering is important, since we may skip an aminoacid of an already cleaved protein
    met_aas = [x for x in pos_first_aas if graph.vs[x]["aminoacid"] == "M" and graph.vs[x]["position"] == 1]

    # Get possible features 
    features = [graph.es.select(_source=start, _target=x)["qualifiers"] for x in met_aas]
    # Get the next nodes which should be after them
    targets = [graph.neighbors(x, mode="out") for x in met_aas]

    # Get current number of edges
    cur_count = graph.ecount()


    # Generate the edges list
    edge_list = [(start, y) for x in targets for y in x]
    qualifiers_list = [[init_met_feature, *features[x_idx][y_idx]] for x_idx, x in enumerate(targets) for y_idx, _ in enumerate(x)]

    # Add new edges (bulk)
    graph.add_edges(edge_list)


    graph.es[cur_count:]["qualifiers"] = qualifiers_list 
    

def get_isoform_dict(comments, accession):
    ''' Get all Isoforms from Comment Section '''
    # TODO comment make nice and quick
    d = {}
    num_of_isoforms = 0

    for isoforms in [x for x in comments if x.startswith("ALTERNATIVE PRODUCTS:")]:
        isoforms = isoforms[len("ALTERNATIVE PRODUCTS:"):]

        entries = isoforms.split(";")
        for e_idx, entry in enumerate(entries):
            entry = entry.strip()
            if entry:
                key, value = entry.split("=", 1)

                # Parse Information
                if key.lower() == "named isoforms":
                    num_of_isoforms = int(value)

                reference_info = ""
                if "{" in value:
                    idx = value.index("{")
                    value = value[:idx].strip()
                    reference_info = value[idx:]

                offset = 1
                if key.lower() == "name":
                    iso_key, iso_value = entries[e_idx+offset].strip().split("=", 1) # TODO There are comma seperated Synonyms!? see ADAM22=G07, 22g(D26D27)+29.3
                    if not iso_key.lower() == "isoid":
                        offset = 2
                        iso_key, iso_value = entries[e_idx+offset].strip().split("=")
                        assert iso_key.lower() == "isoid"
                    seq_key, seq_value = entries[e_idx+offset + 1].strip().split("=")
                    assert seq_key.lower() == "sequence"


                    if "," in iso_value: # BUG/Feature in embl.  Some IDs are not unique (see e.g. P12821-3, P22966-1)
                        # We simply take the first occurence
                        iso_value = iso_value.split(",", 1)[0].strip()

                    d[value] = (iso_value, [x.strip() for x in seq_value.split(",")], reference_info)
                    # if "Displayed" in d[value][1]: # Copy original accession also into isoforms
                    #     d[accession] = (iso_value, [x.strip() for x in seq_value.split(",")], reference_info)
                    #     # TODO can we then simply delete the QXXXXX-Y value in dict?


    return d, num_of_isoforms



def _generate_canonical_graph(sequence: str, acc: str):

    graph = igraph.Graph(directed=True) 

    # Initialize the graph with the sequence
    graph.add_vertices(len(sequence) + 2) 
    graph.add_edges([(x1, x1+1) for x1 in range(len(sequence) + 1)])

    # Add their amino acid to 
    graph.vs["aminoacid"] = ["__start__", *[x for x in sequence], "__end__"]

    # Add position attributes to nodes as well as from which accesion they originate
    graph.vs["position"] = list(range(len(sequence) + 2))
    graph.vs["accession"] = [acc, *[acc]*len(sequence), acc] # TODO is this okay for every entry? to set them to acc
    graph.vs["isoform_accession"] = [None, *[None]*len(sequence), None]
    graph.vs["isoform_position"] = [None, *[None]*len(sequence), None] # Position acording to the isoform


    # Add Information about the Features used
    graph.es["qualifiers"] = [[] for _ in range(len(sequence) + 2)]

    return graph



def _combine_vertices(list_a, list_b):
    out_d = {}

    for a in list_a:
        if a["isoform_accession"] not in out_d:
            out_d[a["isoform_accession"]] = dict(inn=[a])
        else: 
            print("Key is multiple times in it") # TODO we cannot simply associate via accession!?!?
    for b in list_b:
        if b["isoform_accession"] not in out_d:
            out_d[b["isoform_accession"]] = dict(out=[b])
        else:
            if "out" not in out_d[b["isoform_accession"]]:
                out_d[b["isoform_accession"]]["out"] = [b]
            else:    
                out_d[b["isoform_accession"]]["out"].append(b)
            if len(out_d[b["isoform_accession"]]) > 2:
                print("multiple time??")

    k = list(out_d.items())
    a_s = [x[1]["inn"] if "inn" in x[1] else [] for x in k]
    b_s = [x[1]["out"] if "out" in x[1] else [] for x in k]

    return a_s, b_s




def _execute_variant(graph, variant_feature):
    text = variant_feature.qualifiers["note"] # check if missing or if we add another node



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
            vertices_before_raw = list(graph.vs.select(isoform_position=aa_before, isoform_accession=variant_feature.ref))
            vertices_after_raw = list(graph.vs.select(isoform_position=aa_after, isoform_accession=variant_feature.ref))

        if len(vertices_before_raw) == 0 or len(vertices_after_raw) == 0:
            print("fdfd")

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
        graph.add_edges( [x[0] for x in edge_list] )
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
            vertices_before_raw = list(graph.vs.select(isoform_position=aa_before, isoform_accession=variant_feature.ref))
            vertices_after_raw = list(graph.vs.select(isoform_position=aa_after, isoform_accession=variant_feature.ref))

        if len(vertices_before_raw) == 0 or len(vertices_after_raw) == 0:
            print("fdfd")

        vertices_before, vertices_after = _combine_vertices(vertices_before_raw, vertices_after_raw)

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
                graph.add_edges( [(n, y_idcs[idx+1])] )

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
            graph.add_edges( [x[0] for x in edge_list] )
            graph.es[cur_edges:]["qualifiers"] = [x[1] for x in edge_list]

            




def _execute_var_seq(isoforms, graph, sequence: str, var_seqs_features, displayed_accession): # Isoforms

    execute_isoforms = {}

    for f in var_seqs_features:
        
        if "note" not in f.qualifiers:
            print("Some feature tables do not contain information about all isoforms for {}".format(displayed_accession))
            return

        # get isoform information
        note = f.qualifiers["note"]
        isoform_isoids = note[note.index("(")+1 + 3:note.rfind(")")] # +3 to remove "in "
        isoform_isoids = isoform_isoids.replace("isoform", "").replace(" and ", ",")
        for isoid in isoform_isoids.split(","):
            isoid = isoid.strip()
            if isoid not in isoforms:
                print("Isoform not found in specification. Skipping all isoforms for: {}".format(displayed_accession))
                return

            if "Displayed" in isoforms[isoid]:
                print("Isoform modification for canonical sequence found. Skipping all isoforms for: {}".format(displayed_accession))
                return
            
            if isoid not in execute_isoforms:
                execute_isoforms[isoid] = [f]
            else:
                execute_isoforms[isoid].append(f)



    [__start_node__] = graph.vs.select(aminoacid="__start__")
    [__stop_node__]  = graph.vs.select(aminoacid="__end__")
    for key in execute_isoforms.keys():
        iso_sequence, iso_orig_pos, iso_pos = _create_isoform_lists(isoforms[key][0], execute_isoforms[key], sequence)

        # Add the new sequence to the graph
        cur_nodes = graph.vcount()
        cur_edges = graph.ecount()
        graph.add_vertices(len(iso_sequence))
        nodes_indices = graph.vs[cur_nodes:].indices
        graph.add_edges([(nodes_indices[idx], nodes_indices[idx+1]) for idx, _ in enumerate(nodes_indices[:-1])])
        edges_indices = graph.es[cur_edges:].indices

        # Bulk add of all information about nodes!
        graph.vs[nodes_indices]["aminoacid"] = [x for x in iso_sequence] # Adding the Amino Acid!!!
        graph.vs[nodes_indices]["position"] = iso_orig_pos # Here the value None may be present, indicating new sequences!
        graph.vs[nodes_indices]["isoform_position"] = iso_pos # Position acording to the isoform
        graph.vs[nodes_indices]["isoform_accession"] = [isoforms[key][0]]*len(iso_sequence) # Isoform accession like "QXXXXX-3"
        graph.vs[nodes_indices]["accession"] = [displayed_accession]*len(iso_sequence) # The original accession (from canonical)

        # Bulk add Information about edges 
        graph.es[edges_indices]["qualifiers"] = [[]]*(len(edges_indices)-1)

        # Add two special edges between start and end of sequence
        graph.add_edges([(__start_node__, graph.vs[cur_nodes]), (graph.vs[-1], __stop_node__)])
        qualifiers = execute_isoforms[key]
        graph.es[-2:]["qualifiers"] = [qualifiers, []]


def _create_isoform_lists(isoform_accession, feature_list, sequence: str):

    sorted_features = sorted(feature_list, key= lambda x: x.location.start)

    for idx, _ in enumerate(feature_list[:-1]):
        if feature_list[idx].location.end > feature_list[idx+1].location.start:
            print("Isoform information for accession {} overlap! Returning no sequence!".format(isoform_accession))
            return "", [], []

    orig_positions = list(range(1, len(sequence)+1))

    for f in sorted_features[::-1]:
        text = f.qualifiers["note"]
        if text.lower().startswith("missing"):
            # Missing is set, it needs to be removed
            sequence = sequence[:f.location.start] + sequence[f.location.end:]
            orig_positions = orig_positions[:f.location.start] + orig_positions[f.location.end:]
        
        else:
            # Get   X -> Y   Information
            idx = text.find("(")
            if idx != -1:
                text = text[:idx]
            xy = text.split("->")
            assert len(xy) == 2
            y = xy[1].strip().replace(" ", "")
            # Replacing sequence!
            sequence = sequence[:f.location.start] + y + sequence[f.location.end:]
            orig_positions = orig_positions[:f.location.start] + [None]*len(y) + orig_positions[f.location.end:]

    # Return the new sequence, original positions of NOT replaced amino acids and the isoform positions
    return sequence, orig_positions, list(range(1, len(sequence)+1))


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


    # TODO should the signal peptide be exactly the same as in canonical? Or can we leave it as is?
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




# TODO parse note Missing or replace information! via method!

def generate_graph(entry_queue, graph_queue):
    while True:
        # Get next item
        try: 
            entry = entry_queue.get(timeout=1)
        except Exception:
            continue


        # Sort features according to their type
        sorted_features = {}
        for f in entry.features:
            if f.type not in sorted_features:
                sorted_features[f.type] = [f]
            else: 
                sorted_features[f.type].append(f)


        # Get Isoform information
        isoforms, num_of_isoforms = get_isoform_dict(entry.comments, entry.accessions[0])

        # Generate canonical graph (Initialization)
        graph = _generate_canonical_graph(entry.sequence, entry.accessions[0])


        # VAR_SEQ (isoforms) need to be executed at once and before all other variations
        if "VAR_SEQ" in sorted_features:
            _execute_var_seq(isoforms, graph, entry.sequence, sorted_features["VAR_SEQ"], entry.accessions[0])

        if "INIT_MET" in sorted_features: # NOTE: TODO DONE? 
            for f in sorted_features["INIT_MET"]:
                _execute_init_met(graph, f)

        if "SIGNAL" in sorted_features: # NOTE: TODO DONE? 
            for f in sorted_features["SIGNAL"]:
                _execute_signal(graph, f)

        if "VARIANT" in sorted_features: # NOTE: TODO DONE? 
            for f in sorted_features["VARIANT"]:
                _execute_variant(graph, f)


        # Add generated Graph into the next Queue
        graph_queue.put(graph)



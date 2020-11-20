
import igraph



def _execute_init_met(graph, init_met_feature):
    ''' DONE? NOTE This may contain flaws! '''


    # Get start
    [start] = graph.vs.select(aminoacid="__start__")

    # Get all neigbors of start
    pos_first_aas = graph.neighbors(start, mode="out")

    # Filtering is important, since we may skip an aminoacid of an already cleaved protein
    met_aas = [x for x in pos_first_aas if graph.vs[x]["aminoacid"] == "M" and graph.vs[x]["position"] == 1]

    if len(met_aas) > 1:
        # TODO  CHECK if it works
        print("TODO Check Working of INI_MET with multiple entries")

    # Get possible features 
    features = [graph.es.select(_source=start, _target=x)["qualifiers"] for x in met_aas]
    # Get the next nodes which should be after them
    targets = [graph.neighbors(x, mode="out") for x in met_aas]

    # Get current number of edges
    cur_count = graph.ecount()

    # Get qualifier info
    qualifiers = init_met_feature.qualifiers
    if hasattr(init_met_feature, "id"):
        qualifiers["id"] = init_met_feature.id


    # Generate the edges list
    edge_list = [(start, y) for x in targets for y in x]
    qualifiers_list = [[qualifiers, *features[x_idx][y_idx]] for x_idx, x in enumerate(targets) for y_idx, _ in enumerate(x)]

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
                    iso_key, iso_value = entries[e_idx+offset].strip().split("=")
                    if not iso_key.lower() == "isoid":
                        offset = 2
                        iso_key, iso_value = entries[e_idx+offset].strip().split("=")
                        assert iso_key.lower() == "isoid"
                    seq_key, seq_value = entries[e_idx+offset + 1].strip().split("=")
                    assert seq_key.lower() == "sequence"
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
    graph.vs["accession"] = [None, *[acc]*len(sequence), None]

    # Add Information about the Features used
    graph.es["qualifiers"] = [[] for _ in range(len(sequence) + 2)]

    return graph


def _execute_variant(graph, variant_feature):
    text = variant_feature.qualifiers["note"] # check if missing or if we add another node


    if variant_feature.ref is not None:
        print("TODO: Protein Isoform {} with VARIANT: {} is currently ignored".format(variant_feature.ref, variant_feature.id))
        return # TODO we need to handle isoforms somehow here!
        # Luckily it gives us such information by setting the ref attribute
        


    # TODO what about nodes which are sequentielly replaced
    # add missing node here use select for all available previous nodes (using acc AND position)!!!

    if text.lower().startswith("missing"):
        # A sequence is missing! Just append an edge and its information
        # NOTE: Shifted by 1 due to the __start__ node at 0
        aa_before = variant_feature.location.start + 0
        aa_after = variant_feature.location.end + 1
        graph.add_edges( [(aa_before, aa_after)] )

        # Add Information about the new edge (we are here using a variant)
        qualifiers = variant_feature.qualifiers # Get information
        if hasattr(variant_feature, "id"):
            qualifiers["id"] = variant_feature.id
        graph.es[-1]["qualifiers"] = qualifiers # add them to the last added edge

    else:

        # Get   X -> Y   Information
        idx = text.find("(")
        if idx != -1:
            text = text[:idx]
        xy = text.split("->")
        assert len(xy) == 2
        y = xy[1].strip()
        

        # Get start and end position
        # NOTE: Shifted by 1 due to the __start__ node at 0
        aa_before = variant_feature.location.start + 0
        aa_after = variant_feature.location.end + 1


        # Append new node in graph and mapping

        # Case for one or multiple amino acids TODO
        # Add each individual amino acid as a node
        y_idcs = []
        for entry in y:
            vertex = graph.add_vertex()
            graph.vs[vertex.index]["aminoacid"] = entry
            # TODO add position and acc
            y_idcs.append(vertex.index)

        # Add edges between them:
        for idx, n in enumerate(y_idcs[:-1]):
            graph.add_edges( [(n, y_idcs[idx+1])] )





        first_node, last_node = y_idcs[0], y_idcs[-1]
        graph.add_edges( [(aa_before, first_node), (last_node, aa_after)] )


        qualifiers = variant_feature.qualifiers
        if hasattr(variant_feature, "id"):
            qualifiers["id"] = variant_feature.id
        graph.es[-2]["qualifiers"] = qualifiers 



def _execute_var_seq(isoforms, graph, sequence: str, var_seqs_features, displayed_accession): # Isoforms

    execute_isoforms = {}

    # if len(var_seqs_features) < 3:
    #     return

    for f in var_seqs_features:
        
        if "note" not in f.qualifiers:
            print("Some feature tables do not contain information about the all isoforms for {}".format(displayed_accession))
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
            y = xy[1].strip()

            # Replacing sequence!
            sequence = sequence[:f.location.start] + y + sequence[f.location.end:]
            orig_positions = orig_positions[:f.location.start] + [None]*len(y) + orig_positions[f.location.end:]

    # Return the new sequence, original positions of NOT replaced amino acids and the isoform positions
    return sequence, orig_positions, list(range(1, len(sequence)+1))








# TODO parse note Missing or replace information! via method!

def generate_graph(entry_queue, graph_queue):

    note_dict = {}
    iso_dict = {}

    while True:

        try: 
            entry = entry_queue.get(timeout=1)
        except Exception:
            pass


        isoforms, num_of_isoforms = get_isoform_dict(entry.comments, entry.accessions[0])

        note_dict[entry.accessions[0]] = num_of_isoforms
        iso_dict[entry.accessions[0]] = isoforms


        # TODO lists?!?
        graph = _generate_canonical_graph(entry.sequence, entry.accessions[0])


        # TODO var seqs (isoforms) need to be executed at once!
        var_seqs = []
        for f in entry.features:
            if f.type == "VAR_SEQ":
                var_seqs.append(f)
        _execute_var_seq(isoforms, graph, entry.sequence, var_seqs, entry.accessions[0])




        for f in entry.features:

            if f.type == "INIT_MET" and False: # NOTE: TODO DONE? 
                _execute_init_met(graph, f)
                pass
            



                pass

            if f.type == "VARIANT" and False: # TODO Debugging skip
                # TODO ISOFORMS BEFORE VARIANTS! (needed!)

                # TODO make it easier than executing each one seperately?
                _execute_variant(graph, f)

        graph_queue.put(graph)



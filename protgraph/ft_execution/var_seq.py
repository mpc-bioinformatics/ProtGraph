from protgraph.ft_execution import get_content


def execute_var_seq(
    isoforms, graph, sequence: str, var_seqs_features, displayed_accession
):
    """
    Executes the Feature Table Information (with parsed comments, retrieved by the dict isoforms)
    VAR_SEQ to generate chains of nodes and edges for the corresponding isoforms.

    NOTE: This transforms the graph without returning it.

    Following Keys are set here:
    Nodes: "isoform_accession", "isoform_position"
    Edges: "qualifiers" ( -> adds VAR_SEQ)
    """
    # First sort all isoforms
    execute_isoforms = {}
    # For each feature
    for f in var_seqs_features:
        # Skip the complete isoform generation if no information is available
        if "note" not in f.qualifiers:
            print(
                "Some feature tables do not contain information "
                "about all isoforms for {}".format(displayed_accession)
            )
            return

        # Get isoform information
        note = f.qualifiers["note"]
        isoform_isoids = note[
            note.index("(") + 1 + 3: note.rfind(")")  # +3 to remove "in "
        ]
        # Replace the words "isoform" and "and"
        isoform_isoids = isoform_isoids.replace("isoform", "").replace(" and ", ", ")
        # For each comma entry:
        for isoid in isoform_isoids.split(", "):
            isoid = isoid.strip()
            # Skip if we did not found a specification
            if isoid not in isoforms:
                print(
                    "Isoform not found in specification. "
                    "Skipping all isoforms for: {}".format(displayed_accession)
                )
                return
            # Skip if we found an isoform as canonical!? Maybe the entry is corrupted?
            if "Displayed" in isoforms[isoid]:
                print(
                    "Isoform modification for canonical sequence "
                    "found. Skipping all isoforms for: {}".format(displayed_accession)
                )
                return

            # We add (or append) it to our dictionary with its isoid
            if isoid not in execute_isoforms:
                execute_isoforms[isoid] = [f]
            else:
                execute_isoforms[isoid].append(f)

    # Get exactly the ONLY start and stop node here
    [__start_node__] = graph.vs.select(aminoacid="__start__")
    [__stop_node__] = graph.vs.select(aminoacid="__end__")

    # Execute for each isoform, all information at once!
    for key in execute_isoforms.keys():
        # Get the isoform sequence (similar as to the canonical form), its
        # position as isofrom as well as its original position (none isoform)
        iso_sequence, iso_orig_pos, iso_pos = _create_isoform_lists(
            isoforms[key][0], execute_isoforms[key], sequence
        )

        # Bulk add the new sequence to the graph (similar to the code in canonical)
        cur_nodes = graph.vcount()
        graph.add_vertices(len(iso_sequence))
        nodes_indices = graph.vs[cur_nodes:].indices
        graph.add_edges([(nodes_indices[idx], nodes_indices[idx + 1])for idx, _ in enumerate(nodes_indices[:-1])])

        # Bulk add of all information about these new nodes!
        graph.vs[nodes_indices]["aminoacid"] = [x for x in iso_sequence]  # Adding the aminoacid
        graph.vs[nodes_indices]["position"] = iso_orig_pos  # Position according to  (none indicates new seqs)
        graph.vs[nodes_indices]["isoform_position"] = iso_pos  # Position according to the isoform
        graph.vs[nodes_indices]["isoform_accession"] = [isoforms[key][0]] * len(iso_sequence)  # Iso-Acc. ("PXXXXX-3")
        graph.vs[nodes_indices]["accession"] = [displayed_accession] * len(iso_sequence)  # The original accession

        # Add two special edges between start and end of sequence
        graph.add_edges([(__start_node__, graph.vs[cur_nodes]), (graph.vs[-1], __stop_node__)])
        qualifiers = execute_isoforms[key]
        # Add qualifier information here for the FT: VAR_SEQ
        graph.es[-2:]["qualifiers"] = [qualifiers, []]  # need to be set like this for igraph!?


def _create_isoform_lists(isoform_accession, feature_list, sequence: str):
    """ TODO comments """

    sorted_features = sorted(feature_list, key=lambda x: x.location.start)

    for idx, _ in enumerate(feature_list[:-1]):
        if feature_list[idx].location.end > feature_list[idx + 1].location.start:
            print(
                "Isoform information for accession {} overlap! "
                "Returning no sequence!".format(isoform_accession)
            )
            return "", [], []

    orig_positions = list(range(1, len(sequence) + 1))

    for f in sorted_features[::-1]:
        text = f.qualifiers["note"]
        if text.lower().startswith("missing"):
            # Missing is set, it needs to be removed
            sequence = sequence[: f.location.start] + sequence[f.location.end:]
            orig_positions = (orig_positions[: f.location.start] + orig_positions[f.location.end:])

        else:
            # Get to be replaced amino_acids
            y = get_content(text, "(", "->")
            # Replacing sequence!
            sequence = sequence[: f.location.start] + y + sequence[f.location.end:]
            orig_positions = orig_positions[: f.location.start] + [None] * len(y) + orig_positions[f.location.end:]

    # Return the new sequence, original positions of NOT replaced amino acids and the isoform positions
    return sequence, orig_positions, list(range(1, len(sequence) + 1))


def _get_isoforms_of_entry(comments, accession):
    """ Get and parse all isoforms from the comment section """
    # TODO make nice and quick
    d = {}
    num_of_isoforms = 0

    # for each comment line starting with: "ALTERNATIVE PRODUCTS:"
    for isoforms in [x for x in comments if x.startswith("ALTERNATIVE PRODUCTS:")]:
        # Get its text
        isoforms = isoforms[len("ALTERNATIVE PRODUCTS:"):]
        # and split each entry
        entries = isoforms.split(";")
        # Each entry consist of "Name", "Synonyms", "IsoId", "Sequence"
        # Synonyms may be missing
        # for each of the entries:
        e_idx = 0
        # for e_idx, entry in enumerate(entries):
        while e_idx < len(entries):
            entry = entries[e_idx].strip()
            # First check if it is not empty
            if entry:
                # And then parse its information
                key, value = entry.split("=", 1)

                # Here we retrieve the num of isoforms from the entry
                if key.lower() == "named isoforms":
                    num_of_isoforms = int(value)

                # Some VAR_SEQ contain reference info.
                # We retrieve them here
                # TODO this is currently not used anywhere
                reference_info = ""
                if "{" in value:
                    idx = value.index("{")
                    value = value[:idx].strip()
                    reference_info = value[idx:]

                # Parse the information here
                if key.lower() == "name":
                    # We found an entry
                    # We now try to retrieve IsoId and sequence
                    # TODO There are comma seperated Synonyms!? see ADAM22=G07, 22g(D26D27)+29.3
                    e_idx += 1  # Increase by one to get the IsoID information
                    iso_key, iso_value = entries[e_idx].strip().split("=", 1)
                    if not iso_key.lower() == "isoid":
                        e_idx += 1  # if IsoId is still not next we simply increase it again
                        iso_key, iso_value = entries[e_idx].strip().split("=")
                        assert iso_key.lower() == "isoid"
                    e_idx += 1  # Increase it again to retrieve the sequence information
                    # TODO the sequence information is not used anywhere
                    seq_key, seq_value = entries[e_idx].strip().split("=")
                    assert seq_key.lower() == "sequence"

                    if "," in iso_value:
                        # BUG/Feature in embl. Some IDs are not unique (see e.g. P12821-3, P22966-1)
                        # We simply take the first occurence
                        iso_value = iso_value.split(",", 1)[0].strip()

                    # At last add the information to the dict
                    d[value] = (
                        iso_value,
                        [x.strip() for x in seq_value.split(",")],
                        reference_info,
                    )
            e_idx += 1  # Increase for the next iteration

    return d, num_of_isoforms

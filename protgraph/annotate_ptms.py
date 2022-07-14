from Bio.SwissProt import FeatureLocation, FeatureTable, UnknownPosition


def annotate_ptms(graph_entry, var_mods: list, fix_mods: list, mass_factor: int):
    """
    Annotates a delta_mass which should be considered during the generation of weights

    E.G: var_mods contains [("M",  "+15.994915")] (via CLI "-vb 'M:+15.994915'")
    will duplicate all M's, including and ONLY add to the duplicated M's the delta mass (and a new Qualifier)

    Multiple PTMs can be added, however, NOTE that modifications on the same aminoacid are superseeded by
    the newest delta mass.

    Analogously for fix_mods, where all aminoacids get the fixed modifcation (without duplicating)

    Following Keys are set here:
    Nodes: "delta_mass" ( -> the already in mass-factor unit given delta mass )
    Edges: "qualifiers" ( -> adds Var/Fixed Modification)
    """
    var_mods = {x: y for x, y in var_mods} if var_mods else {}
    fix_mods = {x: y for x, y in fix_mods} if fix_mods else {}

    # FIXMOD Handling
    _apply_fixmod(graph_entry, fix_mods, mass_factor)


    # VARMOD Handling
    _apply_varmod(graph_entry, var_mods, mass_factor)


def _apply_fixmod(graph_entry, fix_mods, mass_factor):
    # Apply for each aminoa acid fixed modification
    for aa, delta in fix_mods.items():
        if aa not in ["NTERM", "CTERM"]:
            # CASE AAs:
            # Simply search for every aminoacid and set its delta_mass statically!
            graph_entry.vs.select(aminoacid=aa)["delta_mass"] = delta*mass_factor

            # Append Qualifier Information to ingoing edges
            for n in graph_entry.vs.select(aminoacid=aa):
                for in_edge in n.in_edges():
                    spos, epos = __get_node_pos(n)

                    # Append qualifier
                    if "qualifiers" not in in_edge.attributes() or in_edge["qualifiers"] is None:
                        in_edge["qualifiers"] = [
                            FeatureTable(
                                location=FeatureLocation(spos, epos), type="FIXMOD", strand=None, id=None,
                                qualifiers=dict(note="{}:{}".format(aa, delta))
                            )
                        ]
                    else:
                        in_edge["qualifiers"].append(
                            FeatureTable(
                                location=FeatureLocation(spos, epos), type="FIXMOD", strand=None, id=None,
                                qualifiers=dict(note="{}:{}".format(aa, delta)))
                            )
        elif aa == "NTERM":
            # CASE N-TERMINAL:
            [start_node] = graph_entry.vs.select(aminoacid="__start__")

            new_start = graph_entry.add_vertex()
            # Copy values to new start node
            for key, val in start_node.attributes().items():
                new_start[key] = val
            
            # Set amino acid of old start node to nothing (dereference it as start_node)
            start_node["aminoacid"] = ""
            start_node["delta_mass"] = delta*mass_factor

            # Append qualifier
            new_edge = graph_entry.add_edge(new_start.index, start_node.index)
            new_edge["qualifiers"] = [
                FeatureTable(
                    location=FeatureLocation(0, 1), type="FIXMOD", strand=None, id=None,
                    qualifiers=dict(note="{}:{}".format(aa, delta))
                )
            ]
        elif aa == "CTERM":
            # CASE C-TERMINAL:
            [end_node] = graph_entry.vs.select(aminoacid="__end__")

            new_end = graph_entry.add_vertex()
            # Copy values to new end node
            for key, val in end_node.attributes().items():
                new_end[key] = val
            
            # Set amino acid of old end node to nothing (dereference it as end_node)
            end_node["aminoacid"] = ""
            end_node["delta_mass"] = delta*mass_factor

            # Append new edge
            graph_entry.add_edge(end_node.index, new_end.index)

            # For all ingoing edges add this new qualifier
            for in_edge in end_node.in_edges():
                if "qualifiers" not in in_edge.attributes() or in_edge["qualifiers"] is None:
                    in_edge["qualifiers"] = ([
                        FeatureTable(
                            location=FeatureLocation(end_node["position"], end_node["position"]+1), type="FIXMOD", strand=None, id=None,
                            qualifiers=dict(note="{}:{}".format(aa, delta))
                        )
                    ])
                else:
                    in_edge["qualifiers"].append(
                        FeatureTable(
                            location=FeatureLocation(end_node["position"], end_node["position"]+1), type="FIXMOD", strand=None, id=None,
                            qualifiers=dict(note="{}:{}".format(aa, delta))
                        )
                    )


def _apply_varmod(graph_entry, var_mods, mass_factor):
    for aa, delta in var_mods.items():
        if aa not in ["NTERM", "CTERM"]:
            # Search for every aminoacid and set for the new nodes the delta_mass!
            nodes_to_clone = graph_entry.vs.select(aminoacid=aa)

            # Clone Nodes
            vc = graph_entry.vcount()
            graph_entry.add_vertices(len(nodes_to_clone))
            graph_entry.vs[vc:]["aminoacid"] = nodes_to_clone["aminoacid"]
            graph_entry.vs[vc:]["position"] = nodes_to_clone["position"]
            graph_entry.vs[vc:]["accession"] = nodes_to_clone["accession"]
            if "isoform_accession" in graph_entry.vs[0].attributes():
                graph_entry.vs[vc:]["isoform_accession"] = nodes_to_clone["isoform_accession"]
            if "isoform_position" in graph_entry.vs[0].attributes():
                graph_entry.vs[vc:]["isoform_position"] = nodes_to_clone["isoform_position"]

            # Add Delta Mass to cloned nodes
            graph_entry.vs[vc:]["delta_mass"] = delta*mass_factor

            # Clone the edges
            edges_to_add = []
            qualifier_info = []
            cleaved_info = []
            for offset, node in enumerate(nodes_to_clone):
                spos, epos = __get_node_pos(node)

                for in_edge in node.in_edges():
                    edges_to_add.append((in_edge.source, offset + vc))
                    if "qualifiers" not in in_edge.attributes() or in_edge["qualifiers"] is None:
                        qualifier_info.append([
                            FeatureTable(
                                location=FeatureLocation(spos, epos), type="VARMOD", strand=None, id=None,
                                qualifiers=dict(note="{}:{}".format(aa, delta))
                            )
                        ])
                    else:
                        qualifier_info.append(in_edge["qualifiers"].copy())
                        qualifier_info[-1].append(
                            FeatureTable(
                                location=FeatureLocation(spos, epos), type="VARMOD", strand=None, id=None,
                                qualifiers=dict(note="{}:{}".format(aa, delta))
                            )
                        )
                    if "cleaved" in graph_entry.es[0].attributes():
                        cleaved_info.append(in_edge["cleaved"])

                for out_edge in node.out_edges():
                    edges_to_add.append((offset + vc, out_edge.target))
                    if "qualifiers" in out_edge.attributes():
                        qualifier_info.append(out_edge["qualifiers"])
                    else:
                        qualifier_info.append(None)
                    if "cleaved" in graph_entry.es[0].attributes():
                        cleaved_info.append(out_edge["cleaved"])
            # Add edges in bulk
            ec = graph_entry.ecount()
            graph_entry.add_edges(edges_to_add)
            graph_entry.es[ec:]["qualifiers"] = qualifier_info
            if "cleaved" in graph_entry.es[0].attributes():
                graph_entry.es[ec:]["cleaved"] = cleaved_info

        elif aa == "NTERM":
            # CASE N-TERMINAL:
            # Clone the start node and create a new start node
            [start_node] = graph_entry.vs.select(aminoacid="__start__")
            cloned_start_node = graph_entry.add_vertex()
            new_start_node = graph_entry.add_vertex()

            # Copy values to new start node
            for key, val in start_node.attributes().items():
                cloned_start_node[key] = val
                new_start_node[key] = val

            # De-reference the old start_node and cloned start_node, set deltamass
            cloned_start_node["aminoacid"] = ""
            cloned_start_node["delta_mass"] = delta*mass_factor
            start_node["aminoacid"] = ""

            # Clone edges
            edges_to_add = []
            qualifier_info = []
            cleaved_info = []
            for e in start_node.out_edges():
                edges_to_add.append((cloned_start_node.index, e.target))
                qualifier_info.append(e["qualifiers"])
                if "cleaved" in graph_entry.es[0].attributes():
                    cleaved_info.append(e["cleaved"])

            # connect the new starts together:
            # new_start --> start_node
            # new_start --> cloned_start_node
            edges_to_add.append((new_start_node.index, start_node.index))
            edges_to_add.append((new_start_node.index, cloned_start_node.index))
            qualifier_info.extend([ 
                None,
                [FeatureTable(location=FeatureLocation(0, 1), type="VARMOD", strand=None, id=None,
                    qualifiers=dict(note="{}:{}".format(aa, delta)))]
            ])
            if "cleaved" in graph_entry.es[0].attributes():
                cleaved_info.extend([None, None])

            # Add edges in bulk
            ec = graph_entry.ecount()
            graph_entry.add_edges(edges_to_add)
            graph_entry.es[ec:]["qualifiers"] = qualifier_info
            if "cleaved" in graph_entry.es[0].attributes():
                graph_entry.es[ec:]["cleaved"] = cleaved_info

        elif aa == "CTERM":
            # CASE N-TERMINAL:
            # Clone the end node and create a new end node
            [end_node] = graph_entry.vs.select(aminoacid="__end__")
            cloned_end_node = graph_entry.add_vertex()
            new_end_node = graph_entry.add_vertex()

            # Copy values to new end node
            for key, val in end_node.attributes().items():
                cloned_end_node[key] = val
                new_end_node[key] = val

            # De-reference the old end node and cloned end node, set deltamass
            cloned_end_node["aminoacid"] = ""
            cloned_end_node["delta_mass"] = delta*mass_factor
            end_node["aminoacid"] = ""

            # Clone edges
            edges_to_add = []
            qualifier_info = []
            cleaved_info = []
            for e in end_node.in_edges():
                edges_to_add.append((e.source, cloned_end_node.index))
                if "qualifiers" not in e.attributes() or e["qualifiers"] is None:
                    qualifier_info.append([FeatureTable(location=FeatureLocation(end_node["position"], end_node["position"]+1), type="VARMOD", strand=None, id=None,
                    qualifiers=dict(note="{}:{}".format(aa, delta)))])
                else:
                    qualifier_info.append(e["qualifiers"].copy())
                    qualifier_info[-1].append(FeatureTable(location=FeatureLocation(end_node["position"], end_node["position"]+1), type="VARMOD", strand=None, id=None,
                    qualifiers=dict(note="{}:{}".format(aa, delta))))
                if "cleaved" in graph_entry.es[0].attributes():
                    cleaved_info.append(e["cleaved"])

            # connect the new ends together:
            # end_node --> new_end_node
            # cloned_end_node --> new_end_node
            edges_to_add.append((end_node.index, new_end_node.index))
            edges_to_add.append((cloned_end_node.index, new_end_node.index))
            qualifier_info.extend([ 
                None, None
            ])
            if "cleaved" in graph_entry.es[0].attributes():
                cleaved_info.extend([None, None])

            # Add edges in bulk
            ec = graph_entry.ecount()
            graph_entry.add_edges(edges_to_add)
            graph_entry.es[ec:]["qualifiers"] = qualifier_info
            if "cleaved" in graph_entry.es[0].attributes():
                graph_entry.es[ec:]["cleaved"] = cleaved_info


def __get_node_pos(n):
    """ Parse Position (depending wheather we know it from the isoform or canonical sequence) """
    if "isoform_position" in n.attributes() and n["isoform_position"] is not None:
        if n["isoform_position"] is None:
            return UnknownPosition(), UnknownPosition()
        else:
            return n["isoform_position"] - 1, n["isoform_position"]
    else:
        if n["position"] is None:
            return UnknownPosition(), UnknownPosition()
        else:
            return n["position"] - 1, n["position"]

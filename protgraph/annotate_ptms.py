from collections import defaultdict

import igraph
from Bio.SeqFeature import FeatureLocation, UnknownPosition
from Bio.SwissProt import FeatureTable


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
    d_var_mods, d_fix_mods = defaultdict(list), defaultdict(list)
    if var_mods is not None:
        for k, v in var_mods:
            d_var_mods[k].append(v)
    if fix_mods is not None:
        for k, v in fix_mods:
            d_fix_mods[k].append(v)

    # FIXMOD Handling
    _apply_fixmod(graph_entry, d_fix_mods, mass_factor)

    # VARMOD Handling
    _apply_varmod(graph_entry, d_var_mods, mass_factor)


def _flatten_list(double_list: list):
    return [x for sl in double_list for x in sl]


def _create_feature(start_pos, end_pos, mod_type, amino_acid, delta):
    ''' A wrapper to return a VARMOD or FIXMOD feature '''
    return FeatureTable(
        location=FeatureLocation(start_pos, end_pos), type=mod_type, id=None,
        qualifiers=dict(note="{}:{}".format(amino_acid, delta))
    )


def _append_feature_in_edge(attribute, edge, start_pos, end_pos, mod_type, amino_acid, delta):
    ''' Wrapper function to append an feature '''
    if attribute not in edge.attributes() or edge[attribute] is None:
        edge[attribute] = [_create_feature(start_pos, end_pos, mod_type, amino_acid, delta)]
    else:
        edge[attribute].append(_create_feature(start_pos, end_pos, mod_type, amino_acid, delta))


def _apply_fixmod(graph_entry, fix_mods, mass_factor):
    # Apply for each aminoa acid fixed modification
    for aa, delta in fix_mods.items():
        delta = delta[0]  # We only take the first fix modification. Only a single fixed modification can be applied
        if aa not in ["NPEPTERM", "CPEPTERM", "NPROTERM", "CPROTERM"] \
                and not (aa.startswith("NPEP") or aa.startswith("CPEP")):
            # CASE AAs:
            # Simply search for every aminoacid and set its delta_mass statically!
            graph_entry.vs.select(aminoacid=aa)["delta_mass"] = delta*mass_factor

            # Append Qualifier Information to ingoing edges
            for n in graph_entry.vs.select(aminoacid=aa):
                for in_edge in n.in_edges():
                    spos, epos = __get_node_pos(n)

                    # Append qualifier
                    _append_feature_in_edge("qualifiers", in_edge, spos, epos, "FIXMOD", aa, delta)

    if "NPEPTERM" in fix_mods:
        deltas = fix_mods["NPEPTERM"]
        delta = deltas[0]
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
        new_edge["qualifiers"] = [_create_feature(0, 1, "FIXMOD", "NPEPTERM", delta)]

    if "CPEPTERM" in fix_mods:
        deltas = fix_mods["CPEPTERM"]
        delta = deltas[0]
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
            _append_feature_in_edge(
                "qualifiers", in_edge, end_node["position"], end_node["position"] + 1, "FIXMOD", "CPEPTERM", delta
            )

    if "NPROTERM" in fix_mods:
        deltas = fix_mods["NPROTERM"]
        delta = deltas[0]
        # CASE N-TERMINAL for Protein-N

        # Add new nodes
        nodes_to_update = get_first_last_nodes_from_protein(graph_entry, start_end_aa="__start__", min_max=min)

        # Updata delta masses
        for n in nodes_to_update:
            if "delta_mass" in graph_entry.vs[0].attributes() and graph_entry.vs[n.index]["delta_mass"] is not None:
                graph_entry.vs[n.index]["delta_mass"] += delta*mass_factor
            else:
                graph_entry.vs[n.index]["delta_mass"] = delta*mass_factor

        for n in nodes_to_update:
            for e in n.in_edges():
                spos, epos = __get_node_pos(n)
                for in_edge in n.in_edges():
                    _append_feature_in_edge("qualifiers", e, spos, epos, "FIXMOD", "NPROTERM", delta)

    if "CPROTERM" in fix_mods:
        deltas = fix_mods["CPROTERM"]
        delta = deltas[0]
        # CASE C-TERMINAL for Protein-C

        # Add new nodes
        nodes_to_update = get_first_last_nodes_from_protein(graph_entry, start_end_aa="__end__", min_max=max)

        # Updata delta masses
        for n in nodes_to_update:
            if "delta_mass" in graph_entry.vs[0].attributes() and graph_entry.vs[n.index]["delta_mass"] is not None:
                graph_entry.vs[n.index]["delta_mass"] += delta*mass_factor
            else:
                graph_entry.vs[n.index]["delta_mass"] = delta*mass_factor

        for n in nodes_to_update:
            for e in n.in_edges():
                spos, epos = __get_node_pos(n)
                for in_edge in n.in_edges():
                    _append_feature_in_edge("qualifiers", e, spos, epos, "FIXMOD", "CPROTERM", delta)


def _apply_varmod(graph_entry, var_mods, mass_factor):
    for aa, deltas in var_mods.items():
        if aa not in ["NPEPTERM", "CPEPTERM", "NPROTERM", "CPROTERM"] \
                and not (aa.startswith("NPEP") or aa.startswith("CPEP")):
            # Search for every aminoacid and set for the new nodes the delta_mass!
            nodes_to_clone = graph_entry.vs.select(aminoacid=aa)

            # Clone Nodes
            vc = graph_entry.vcount()
            graph_entry.add_vertices(len(nodes_to_clone)*len(deltas))
            graph_entry.vs[vc:]["aminoacid"] = nodes_to_clone["aminoacid"]*len(deltas)
            graph_entry.vs[vc:]["position"] = nodes_to_clone["position"]*len(deltas)
            graph_entry.vs[vc:]["accession"] = nodes_to_clone["accession"]*len(deltas)
            if "isoform_accession" in graph_entry.vs[0].attributes():
                graph_entry.vs[vc:]["isoform_accession"] = nodes_to_clone["isoform_accession"]*len(deltas)
            if "isoform_position" in graph_entry.vs[0].attributes():
                graph_entry.vs[vc:]["isoform_position"] = nodes_to_clone["isoform_position"]*len(deltas)

            # Add Delta Mass to cloned nodes
            graph_entry.vs[vc:]["delta_mass"] = _flatten_list(
                [[delta*mass_factor]*len(nodes_to_clone) for delta in deltas]
            )

            # Clone the edges
            edges_to_add = []
            qualifier_info = []
            cleaved_info = []
            for offset, node in enumerate(nodes_to_clone):
                spos, epos = __get_node_pos(node)

                for in_edge in node.in_edges():
                    for delta_idx, delta in enumerate(deltas):
                        edges_to_add.append((in_edge.source, offset + len(nodes_to_clone)*delta_idx + vc))
                        qualifier_info.append([_create_feature(spos, epos, "VARMOD", aa, delta)])

                        if "cleaved" in graph_entry.es[0].attributes():
                            cleaved_info.append(in_edge["cleaved"])

                for out_edge in node.out_edges():
                    for delta_idx, delta in enumerate(deltas):
                        edges_to_add.append((offset + len(nodes_to_clone)*delta_idx + vc, out_edge.target))
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

        if aa.startswith("NPEP") or aa.startswith("CPEP"):
            # CASE AAs on C-Term or N-Term
            [start_node] = graph_entry.vs.select(aminoacid="__start__")
            [end_node] = graph_entry.vs.select(aminoacid="__end__")

            if aa.startswith("NPEP"):
                filter_expression_for_cloning = lambda node: any([x["position"] == 0 for x in node.neighbors(mode="in")])  # noqa: E731,E501
                filter_ingoing_edges = lambda x: x == 0  # x --> position of Node # noqa: E731,E501
                filter_outgoing_edges = lambda _: True  # No Filtering of outgoing edges # noqa: E731,E501
            else:
                filter_expression_for_cloning = lambda node: any([x["position"] == end_node["position"] for x in node.neighbors(mode="out")])  # noqa: E731,E501
                filter_ingoing_edges = lambda _: True  # No Filtering of ingoing edges # noqa: E731,E501
                filter_outgoing_edges = lambda x: x == end_node["position"]  # x --> position of Node # noqa: E731,E501

            # Filter Nodes which should be cloned
            nodes_to_check = graph_entry.vs.select(aminoacid=aa[4:])
            nodes_to_clone = []
            for n in nodes_to_check:
                if filter_expression_for_cloning(n):
                    nodes_to_clone.append(n)
            nodes_to_clone = igraph.VertexSeq(graph_entry, [n.index for n in nodes_to_clone])

            # Clone Nodes
            vc = graph_entry.vcount()
            graph_entry.add_vertices(len(nodes_to_clone)*len(deltas))
            graph_entry.vs[vc:]["aminoacid"] = nodes_to_clone["aminoacid"]*len(deltas)
            graph_entry.vs[vc:]["position"] = nodes_to_clone["position"]*len(deltas)
            graph_entry.vs[vc:]["accession"] = nodes_to_clone["accession"]*len(deltas)
            if "isoform_accession" in graph_entry.vs[0].attributes():
                graph_entry.vs[vc:]["isoform_accession"] = nodes_to_clone["isoform_accession"]*len(deltas)
            if "isoform_position" in graph_entry.vs[0].attributes():
                graph_entry.vs[vc:]["isoform_position"] = nodes_to_clone["isoform_position"]*len(deltas)

            # Add Delta Mass to cloned nodes
            graph_entry.vs[vc:]["delta_mass"] = _flatten_list(
                [[delta*mass_factor]*len(nodes_to_clone) for delta in deltas]
            )

            # Clone the edges
            edges_to_add = []
            qualifier_info = []
            cleaved_info = []
            for offset, node in enumerate(nodes_to_clone):
                spos, epos = __get_node_pos(node)

                # Only a connection from the "start_node"
                for in_edge in node.in_edges():
                    if filter_ingoing_edges(in_edge.source_vertex["position"]):
                        for delta_idx, delta in enumerate(deltas):
                            edges_to_add.append((in_edge.source, offset + len(nodes_to_clone)*delta_idx + vc))
                            qualifier_info.append([_create_feature(spos, epos, "VARMOD", aa, delta)])

                            if "cleaved" in graph_entry.es[0].attributes():
                                cleaved_info.append(in_edge["cleaved"])

                # connection to all outgoing edges
                for out_edge in node.out_edges():
                    if filter_outgoing_edges(out_edge.target_vertex["position"]):
                        for delta_idx, delta in enumerate(deltas):
                            edges_to_add.append((offset + len(nodes_to_clone)*delta_idx + vc, out_edge.target))
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

    if "NPEPTERM" in var_mods:
        deltas = var_mods["NPEPTERM"]
        # CASE N-TERMINAL:
        # Clone the start node and create new start nodes
        [start_node] = graph_entry.vs.select(aminoacid="__start__")
        vc = graph_entry.vcount()
        graph_entry.add_vertices(len(deltas) + 1)

        # Copy values to new start nodes
        for key, val in start_node.attributes().items():
            for offset in range(len(deltas) + 1):
                graph_entry.vs[vc + offset][key] = val

        # De-reference the old start_node and all variable modification nodes. Set the delta masses
        start_node["aminoacid"] = ""
        for offset, delta in enumerate(deltas):
            graph_entry.vs[vc + offset]["aminoacid"] = ""
            graph_entry.vs[vc + offset]["delta_mass"] = delta*mass_factor

        # Clone edges
        edges_to_add = []
        qualifier_info = []
        cleaved_info = []
        for e in start_node.out_edges():
            for offset in range(len(deltas)):
                edges_to_add.append((vc + offset, e.target))
                if "qualifiers" in graph_entry.es[0].attributes():
                    qualifier_info.append(e["qualifiers"])
                if "cleaved" in graph_entry.es[0].attributes():
                    cleaved_info.append(e["cleaved"])

        # connect the new starts together:
        # new_start --> cloned_start_node
        # new_start --> second_cloned_start_nodes etc..
        for offset in range(len(deltas)):
            edges_to_add.append((vc + len(deltas), vc + offset))
            qualifier_info.append([_create_feature(0, 1, "VARMOD", "NPEPTERM", delta)])
            if "cleaved" in graph_entry.es[0].attributes():
                cleaved_info.append(None)

        # Also do this for the old start node:
        # new_start --> start_node
        edges_to_add.append((vc + len(deltas), start_node.index))
        qualifier_info.append(None)
        if "cleaved" in graph_entry.es[0].attributes():
            cleaved_info.append(None)

        # Add edges in bulk
        ec = graph_entry.ecount()
        graph_entry.add_edges(edges_to_add)
        graph_entry.es[ec:]["qualifiers"] = qualifier_info
        if "cleaved" in graph_entry.es[0].attributes():
            graph_entry.es[ec:]["cleaved"] = cleaved_info

    if "CPEPTERM" in var_mods:
        deltas = var_mods["CPEPTERM"]
        # CASE N-TERMINAL:
        # Clone the end node and create a new end node
        [end_node] = graph_entry.vs.select(aminoacid="__end__")
        vc = graph_entry.vcount()
        graph_entry.add_vertices(len(deltas) + 1)

        # Copy values to new end nodes
        for key, val in end_node.attributes().items():
            for offset in range(len(deltas) + 1):
                graph_entry.vs[vc + offset][key] = val

        # De-reference the old end:node and all variable modification nodes. Set the delta masses
        end_node["aminoacid"] = ""
        for offset, delta in enumerate(deltas):
            graph_entry.vs[vc + offset]["aminoacid"] = ""
            graph_entry.vs[vc + offset]["delta_mass"] = delta*mass_factor

        # Clone edges
        edges_to_add = []
        qualifier_info = []
        cleaved_info = []
        for e in end_node.in_edges():
            for offset in range(len(deltas)):
                edges_to_add.append((e.source, vc + offset))
                qualifier_info.append([
                    _create_feature(end_node["position"], end_node["position"]+1, "VARMOD", "CPEPTERM", delta)
                ])
                if "cleaved" in graph_entry.es[0].attributes():
                    cleaved_info.append(e["cleaved"])

        # connect the new ends together:
        # cloned_end_node --> new_end_node
        # second_cloned_end_node --> new_end_node
        # ...
        for offset in range(len(deltas)):
            edges_to_add.append((vc + offset, vc + len(deltas)))
            qualifier_info.append(None)
            if "cleaved" in graph_entry.es[0].attributes():
                cleaved_info.append(None)

        # end_node --> new_end_node
        edges_to_add.append((end_node.index, vc + len(deltas)))
        qualifier_info.append(None)
        if "cleaved" in graph_entry.es[0].attributes():
            cleaved_info.append(None)

        # Add edges in bulk
        ec = graph_entry.ecount()
        graph_entry.add_edges(edges_to_add)
        graph_entry.es[ec:]["qualifiers"] = qualifier_info
        if "cleaved" in graph_entry.es[0].attributes():
            graph_entry.es[ec:]["cleaved"] = cleaved_info

    if "NPROTERM" in var_mods:
        deltas = var_mods["NPROTERM"]
        # CASE N-TERMINAL for Protein-N

        # Add new nodes
        nodes_to_clone = get_first_last_nodes_from_protein(graph_entry, start_end_aa="__start__", min_max=min)
        vc = graph_entry.vcount()
        graph_entry.add_vertices(len(nodes_to_clone) * len(deltas))

        # Copy values to new end nodes and update the delta_mass
        offset = 0
        for n in nodes_to_clone:
            for delta in deltas:
                for key, val in n.attributes().items():
                    graph_entry.vs[vc + offset][key] = val
                if "delta_mass" in graph_entry.vs[0].attributes() \
                        and graph_entry.vs[vc + offset]["delta_mass"] is not None:
                    graph_entry.vs[vc + offset]["delta_mass"] += delta*mass_factor
                else:
                    graph_entry.vs[vc + offset]["delta_mass"] = delta*mass_factor
                offset += 1

        # Clone edges
        edges_to_add = []
        qualifier_info = []
        cleaved_info = []

        offset = 0
        for n in nodes_to_clone:
            for delta in deltas:
                for e in n.in_edges():
                    spos, epos = __get_node_pos(n)
                    for in_edge in n.in_edges():
                        edges_to_add.append((in_edge.source, vc + offset))
                        qualifier_info.append([_create_feature(spos, epos, "VARMOD", "NPROTERM", delta)])
                        if "cleaved" in graph_entry.es[0].attributes():
                            cleaved_info.append(in_edge["cleaved"])

                    for out_edge in n.out_edges():
                        edges_to_add.append((vc + offset, out_edge.target))
                        if "qualifiers" in out_edge.attributes():
                            qualifier_info.append(out_edge["qualifiers"])
                        else:
                            qualifier_info.append(None)
                        if "cleaved" in graph_entry.es[0].attributes():
                            cleaved_info.append(out_edge["cleaved"])
                offset += 1

        # Add edges in bulk
        ec = graph_entry.ecount()
        graph_entry.add_edges(edges_to_add)
        graph_entry.es[ec:]["qualifiers"] = qualifier_info
        if "cleaved" in graph_entry.es[0].attributes():
            graph_entry.es[ec:]["cleaved"] = cleaved_info

    if "CPROTERM" in var_mods:
        deltas = var_mods["CPROTERM"]
        # CASE C-TERMINAL for Protein-C

        # Add new nodes
        nodes_to_clone = get_first_last_nodes_from_protein(graph_entry, start_end_aa="__end__", min_max=max)
        vc = graph_entry.vcount()
        graph_entry.add_vertices(len(nodes_to_clone) * len(deltas))

        # Copy values to new end nodes and update the delta_mass
        offset = 0
        for n in nodes_to_clone:
            for delta in deltas:
                for key, val in n.attributes().items():
                    graph_entry.vs[vc + offset][key] = val
                if "delta_mass" in graph_entry.vs[0].attributes() \
                        and graph_entry.vs[vc + offset]["delta_mass"] is not None:
                    graph_entry.vs[vc + offset]["delta_mass"] += delta*mass_factor
                else:
                    graph_entry.vs[vc + offset]["delta_mass"] = delta*mass_factor
                offset += 1

        # Clone edges
        edges_to_add = []
        qualifier_info = []
        cleaved_info = []

        offset = 0
        for n in nodes_to_clone:
            for delta in deltas:
                for e in n.in_edges():
                    spos, epos = __get_node_pos(n)
                    for in_edge in n.in_edges():
                        edges_to_add.append((in_edge.source, vc + offset))
                        qualifier_info.append([_create_feature(spos, epos, "VARMOD", "CPROTERM", delta)])
                        if "cleaved" in graph_entry.es[0].attributes():
                            cleaved_info.append(in_edge["cleaved"])

                    for out_edge in n.out_edges():
                        edges_to_add.append((vc + offset, out_edge.target))
                        if "qualifiers" in out_edge.attributes():
                            qualifier_info.append(out_edge["qualifiers"])
                        else:
                            qualifier_info.append(None)
                        if "cleaved" in graph_entry.es[0].attributes():
                            cleaved_info.append(out_edge["cleaved"])
                offset += 1

        # Add edges in bulk
        ec = graph_entry.ecount()
        graph_entry.add_edges(edges_to_add)
        graph_entry.es[ec:]["qualifiers"] = qualifier_info
        if "cleaved" in graph_entry.es[0].attributes():
            graph_entry.es[ec:]["cleaved"] = cleaved_info


def get_first_last_nodes_from_protein(graph_entry, start_end_aa="__start__", min_max=min):
    ''' Get the beginning of the protein '''
    if "isoform_position" in graph_entry.vs[0].attributes():
        # Get maximum for each isoform end
        isos = list(set(graph_entry.vs["isoform_accession"]))
        nodes = set()
        for iso in isos:
            vseq = graph_entry.vs.select(isoform_accession=iso)
            pos = "position" if iso is None else "isoform_position"  # Decide which position to use
            highetst_num = min_max(
                [
                    x
                    for x, aa in zip(vseq[pos], vseq["aminoacid"])
                    if x and aa != "" and aa != start_end_aa
                ]
            )
            nodes.update(vseq.select(**{pos: highetst_num}))
        return list(nodes)
    else:
        # Get maximum position and just report this
        highetst_num = min_max(
            [
                x
                for x, aa in zip(graph_entry.vs["position"], graph_entry.vs["aminoacid"])
                if x and aa != "" and aa != start_end_aa
            ]
        )
        return list(graph_entry.vs.select(position=highetst_num))


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

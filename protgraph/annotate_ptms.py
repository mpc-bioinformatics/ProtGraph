from Bio.SwissProt import FeatureTable, FeatureLocation, UnknownPosition

def annotate_ptms(graph_entry, var_mods: list, fix_mods: list, mass_factor: int):
    """
    Annotates a delta_mass which should be considered during the generation of weights

    E.G: var_mods contains [("M",  "+15.994915")] (via CLI "-vb 'M:+15.994915'")
    will duplicate all M's, including and ONLY add to the duplicated M's the delta mass (and a new Qualifier)

    Multiple PTMs can be added, however, note that modifications on the same aminoacid are superseeded by 
    the newest delta mass.

    Analogously for fix_mods, where all aminoacids get the fixed modifcation (without duplicating)

    Following Keys are set here:
    Nodes: "delta_mass" ( -> the already in mass-factor unit given delta mass )
    Edges: "qualifiers" ( -> adds Var/Fixed Modification)

    TODO we currently CANNOT handly Cterm or Nterm PTMs!
    """
    var_mods = {x: y for x,y in var_mods} if var_mods else {}
    fix_mods = {x: y for x,y in fix_mods} if fix_mods else {}

    # FIXMOD Handling
    for aa, delta in fix_mods.items():
        # Simply search for every aminoacid and set its delta_mass statically!
        graph_entry.vs.select(aminoacid=aa)["delta_mass"] = delta*mass_factor
        
        # Append Qualifier Information to ingoing edges
        for n in graph_entry.vs.select(aminoacid=aa):
            for in_edge in n.in_edges():
                spos, epos = __get_node_pos(n)

                # Append qualifier
                if "qualifiers" not in in_edge.attributes() or in_edge["qualifiers"] is None:
                    in_edge["qualifiers"] = [FeatureTable(location=FeatureLocation(spos, epos), type="FIXMOD", strand=None, id=None, qualifiers=dict(note="{}:{}".format(aa, delta)))]
                else:
                    in_edge["qualifiers"].append(FeatureTable(location=FeatureLocation(spos, epos), type="FIXMOD", strand=None, id=None, qualifiers=dict(note="{}:{}".format(aa, delta))))


    # VARMOD Handling
    for aa, delta in var_mods.items():
        # Simply search for every aminoacid and set its delta_mass statically!
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
                    qualifier_info.append([FeatureTable(location=FeatureLocation(spos, epos), type="VARMOD", strand=None, id=None, qualifiers=dict(note="{}:{}".format(aa, delta)))])
                else:
                    qualifier_info.append(in_edge["qualifiers"])
                    qualifier_info[-1].append(FeatureTable(location=FeatureLocation(spos, epos), type="VARMOD", strand=None, id=None, qualifiers=dict(note="{}:{}".format(aa, delta))))
                if "cleaved" in graph_entry.es[0].attributes(): cleaved_info.append(in_edge["cleaved"])

        # Add edges in bulk
        ec = graph_entry.ecount()
        graph_entry.add_edges(edges_to_add)
        graph_entry.es[ec:]["qualifiers"] = qualifier_info
        if "cleaved" in graph_entry.es[0].attributes(): graph_entry.es[ec:]["cleaved"] = cleaved_info


def __get_node_pos(n):
    """ Parse Position (depending wheather we know it from the isoform or canonical sequence) """
    if "isoform_position" in n.attributes() and n["isoform_position"] is not None:
        if n["isoform_position"] is None:
            return UnknownPosition(), UnknownPosition() 
        else:
            return n["isoform_position"] -1, n["isoform_position"]
    else:
        if n["position"] is None:
            return UnknownPosition(), UnknownPosition() 
        else:
            return n["position"] - 1, n["position"]

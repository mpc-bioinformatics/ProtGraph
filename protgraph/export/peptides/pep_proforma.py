from decimal import Decimal
from math import isfinite

from protgraph.export.peptides.pep_fasta import PepFasta
from protgraph.graph_collapse_edges import Or


class PepProForma(PepFasta):
    """
    Peptide exporter emitting ProForma 2.0 strings (HUPO-PSI standard notation
    for proteoforms and peptidoforms, LeDuc et al., J. Proteome Res. 2022).

    Writes one tab-separated line per exported peptide path with columns:
        accession, proforma, start, end, misscleavages, qualifiers

    ProtGraph knows modifications only by their delta mass (-fm/-vm), so
    modifications are emitted as mass-delta tags by default ("S[+79.966]"),
    which is valid ProForma. A user-declared mapping (--pep_proforma_mod_names
    "79.966=UNIMOD:21") replaces a delta by a named tag ("S[UNIMOD:21]").
    Terminal modifications (NPEPTERM etc.) use ProForma terminal syntax
    ("[+42.010565]-PEPTIDE" or "[UNIMOD:1]-PEPTIDE").

    A modification is placed on the residue its edge points at: FIXMOD/VARMOD
    features are attached to the in-edges of the (cloned) modified node when
    the graph is generated, so the first residue of the traversed edge's
    target node is the modified one. This also holds for residues introduced
    by VARIANT/MUTAGEN/CONFLICT, whose features carry no usable reference
    location, and for nodes merged after annotation.

    NOTE: as with the FASTA exporter, this exports all possible paths, so make
    sure the traversal can terminate in forseeable future!
    """

    def start_up(self, **kwargs):
        super(PepProForma, self).start_up(**kwargs)
        self.output_file = kwargs["pep_proforma_out"]
        self.mod_names = dict(kwargs["pep_proforma_mod_names"] or [])
        self.warned_or_once = False
        self.write_qualifiers = kwargs["pep_proforma_write_qualifiers"]

    def export(self, prot_graph, queue):
        tsv_header = [
            "accession", "proforma", "start", "end", "misscleavages"
        ] + (["qualifiers"] if self.write_qualifiers else [])

        queue.put((self.output_file, "\t".join(tsv_header) + "\n", True, "w"))
        super().export(prot_graph, queue)

    def export_peptides(self, prot_graph, l_path_nodes, l_path_edges, l_peptide, l_miscleavages, queue):
        entries = ""
        for nodes, edges, misses in zip(l_path_nodes, l_path_edges, l_miscleavages):
            # the first/last node carrying sequence; with terminal modifications
            # applied, nodes[1] can be an empty helper node (see annotate_ptms)
            acc = self._get_accession_or_isoform(prot_graph.vs[nodes[1]])
            start_pos = self._get_position_or_isoform_position(prot_graph.vs[nodes[1]])
            end_pos = self._get_position_or_isoform_position(prot_graph.vs[nodes[-2]], end=True)
            proforma = self._build_proforma(prot_graph, nodes[1:-1], edges)
            l_str_qualifiers = self._get_qualifiers(prot_graph, edges)
            entries += "\t".join(
                [acc, proforma, str(start_pos), str(end_pos), str(misses)]
                + (l_str_qualifiers if self.write_qualifiers else [])
            ) + "\n"

        # w: ensures to have no stale entries, overwriting the file, if it exists.
        queue.put((self.output_file, entries, False, "w"))

    def _build_proforma(self, prot_graph, inner_nodes, edges):
        """ ProForma string for one peptide path. """
        # Get all modifications by iterating over the edges
        n_term, c_term, by_residue = [], [], {}
        for edge_id in edges:
            edge = prot_graph.es[edge_id]
            for key, delta in self._edge_mods(edge):
                tag = self.mod_names.get(float(delta)) or self._delta_tag(delta)
                if key in ("NPEPTERM", "NPROTERM"):
                    n_term.append(tag)
                elif key in ("CPEPTERM", "CPROTERM"):
                    c_term.append(tag)
                else:
                    by_residue.setdefault(edge.target, []).append((key, tag))

        # Build the final proforma sequence
        parts = []
        if n_term:
            parts.append("".join("[{}]".format(t) for t in n_term) + "-")
        for ie in inner_nodes:
            aas = prot_graph.vs[ie]["aminoacid"]
            if len(aas) == 1:
                # Simply append UniMod
                parts.append(aas)
                parts.extend("[{}]".format(t[1]) for t in by_residue.get(ie, ()))
            else:
                # We need to add PTM on correct position in merged node:
                # E.G.  Merged node "EIN", Modification I->+23.123
                # should yield "EI[+23.123]N"  and NOT "EIN[+23.123]"!
                by_aa = {}
                for mods in by_residue.values():
                    for aa, mod in mods:
                        by_aa.setdefault(aa, []).append(mod)
                for aa, mod in by_aa.items():
                    # We can savely use replace as this case only can happen in FIXMODs
                    aas = aas.replace(aa, aa + "".join("[{}]".format(t) for t in mod))
                parts.append(aas)

        if c_term:
            parts.append("-" + "".join("[{}]".format(t) for t in c_term))
        return "".join(parts)

    def _edge_mods(self, edge):
        """ Deduplicated (key, delta) modifications of one traversed edge.

        Or-wrapped qualifiers (collapsed parallel edges) are alternatives, not
        a conjunction: only modifications present in EVERY branch are certain
        for this traversal and get exported; branch-specific ones are dropped
        with a one-time warning rather than fabricated onto every path.
        """
        if "qualifiers" not in edge.attributes():
            return []
        mods = []
        for f in edge["qualifiers"] or []:
            if isinstance(f, Or):
                branches = [self._collect_mods(branch) for branch in f]
                common = set.intersection(*map(set, branches)) if branches else set()
                if any(set(b) - common for b in branches) and not self.warned_or_once:
                    self.warned_or_once = True
                    print(
                        "Warning: modifications differ between collapsed edge alternatives; "
                        "only modifications shared by all alternatives are exported. "
                        "Re-run with --no_collapsing_edges for exact per-path modifications."
                    )
                mods.extend(m for m in dict.fromkeys(x for b in branches for x in b) if m in common)
            elif getattr(f, "type", None) in ("FIXMOD", "VARMOD"):
                mods.append(self._feature_mod(f))
        return list(dict.fromkeys(mods))

    def _collect_mods(self, qualifier):
        out = []
        for f in qualifier or []:
            if isinstance(f, Or):
                for branch in f:
                    out.extend(self._collect_mods(branch))
            elif getattr(f, "type", None) in ("FIXMOD", "VARMOD"):
                out.append(self._feature_mod(f))
        return out

    def _feature_mod(self, feature):
        key, _, delta = feature.qualifiers["note"].rpartition(":")
        return key, delta

    def _delta_tag(self, delta):
        """ '79.966' -> '+79.966': explicitly signed, plain decimal notation.

        ProForma's mass grammar has no exponent and no non-finite forms, so
        '1e-05' is rewritten as '+0.00001' and nan/inf are refused.
        """
        if not isfinite(float(delta)):
            raise ValueError(
                "Cannot express the delta mass '{}' as a ProForma mass tag".format(delta)
            )
        if "e" in delta.lower():
            delta = format(Decimal(delta), "f")
        return delta if delta.startswith("-") else "+" + delta

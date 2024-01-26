from Bio.SeqFeature import UnknownPosition

from protgraph.export.peptides.abstract_peptide_exporter import \
    APeptideExporter
from protgraph.graph_collapse_edges import Or


class PepFasta(APeptideExporter):
    """
    This is a rather 'simple' Peptide-Fasta-File-Exporter into a single file.

    NOTE: It exports all possible paths from start to end, so before executing this
    exporter make sure that it can terminate in forseable future!
    """

    def start_up(self, **kwargs):
        super(PepFasta, self).start_up(**kwargs)
        self.id_gen = self.unique_id_gen(**kwargs)

        # Get output file from configuration
        self.output_file = kwargs["pep_fasta_out"]

    def export_peptides(self, prot_graph, l_path_nodes, l_path_edges, l_peptide, l_miscleavages, queue):
        # Export a batch of peptides at onces
        # Build up the entries for the batch
        entries = ""
        for peptide, nodes, edges, misses in zip(l_peptide, l_path_nodes, l_path_edges, l_miscleavages):
            acc = self._get_accession_or_isoform(prot_graph.vs[nodes[1]])
            start_pos = self._get_position_or_isoform_position(prot_graph.vs[nodes[1]])
            end_pos = self._get_position_or_isoform_position(prot_graph.vs[nodes[-2]], end=True)
            l_str_qualifiers = self._get_qualifiers(prot_graph, edges)
            quali_entries = ",".join(l_str_qualifiers)
            if quali_entries:
                quali_entries = "," + quali_entries
            # Generate header: >pg|ENTRY_ID|UNIPROT_ISO(START:END, mssclvg:MISS; QUALIFIERS)
            entries += "|".join(
                [
                    ">pg",
                    "ID_" + str(next(self.id_gen)),
                    "".join(
                        [
                            acc, "(", str(start_pos), ":", str(end_pos), ",",
                            "mssclvg:", str(misses), quali_entries,  ")"
                        ]
                    )
                ]
            )
            entries += "\n" + '\n'.join(peptide[i:i+60] for i in range(0, len(peptide), 60)) + "\n"

        # put the entries into the queue, so that the main thread can write the actual content.
        queue.put((self.output_file, entries, False, "a"))

    def _get_qualifiers(self, prot_graph, edges):
        """ Generate a list of strings for all avilable qualifiers. """
        if "qualifiers" not in prot_graph.es[0].attributes():
            return []
        str_qualifiers = []
        for qualifier in prot_graph.es[edges]["qualifiers"]:
            if qualifier is None:
                continue
            if len(qualifier) == 0:
                continue

            str_qualifiers.extend(self._map_qualifier_to_string(qualifier))
        return str_qualifiers

    def _map_qualifier_to_string(self, qualifier):
        str_qualifiers = []
        if qualifier is None:
            # This case only happens if there are multiple ways from
            # the start or end node due to the feature beeing there
            # E.G. PROPEP and PEPTIDE have this in common
            return ["None"]
        for f in qualifier:
            if isinstance(f, Or):
                f_qs = []
                for f_q in f:
                    f_qs.append(self._map_qualifier_to_string(f_q))

                str_qualifiers.append(
                    "OR[" + "|".join([",".join(x) for x in f_qs]) + "]"
                )

            elif f.type == "VAR_SEQ":
                continue  # Skip since we encode this differently

            elif f.type == "VARIANT":
                str_qualifiers.append(
                    "VARIANT[" + self._get_location(f) + "," + self._get_variant_mutagen_qualifier(f) + "]"
                )

            elif f.type == "MUTAGEN":
                str_qualifiers.append(
                    "MUTAGEN[" + self._get_location(f) + "," +
                    self._get_variant_mutagen_qualifier(f, stop_codon=":", offset=0) + "]"
                )

            elif f.type == "CONFLICT":
                str_qualifiers.append(
                    "CONFLICT[" + self._get_location(f) + "," + self._get_variant_mutagen_qualifier(f) + "]"
                )

            elif f.type == "SIGNAL":
                str_qualifiers.append(
                    "SIGNAL[" + self._get_location(f) + self._get_pep_prop_chain_sig_id(f) + "]"
                )

            elif f.type == "INIT_MET":
                str_qualifiers.append(
                    "INIT_MET[" + self._get_location(f) + self._get_pep_prop_chain_sig_id(f) + "]"
                )
            elif f.type == "PROPEP":
                str_qualifiers.append(
                    "PROPEP[" + self._get_location(f) + self._get_pep_prop_chain_sig_id(f) + "]"
                )
            elif f.type == "PEPTIDE":
                str_qualifiers.append(
                    "PEPTIDE[" + self._get_location(f) + self._get_pep_prop_chain_sig_id(f) + "]"
                )
            elif f.type == "CHAIN":
                str_qualifiers.append(
                    "CHAIN[" + self._get_location(f) + self._get_pep_prop_chain_sig_id(f) + "]"
                )
            elif f.type == "FIXMOD":
                str_qualifiers.append(
                    "FIXMOD[" + self._get_location(f) + "," + f.qualifiers["note"] + "]"
                )
            elif f.type == "VARMOD":
                str_qualifiers.append(
                    "VARMOD[" + self._get_location(f) + "," + f.qualifiers["note"] + "]"
                )

            else:
                print("Warning: FASTA-Export-Case is not covered: {}".format(f))

        return str_qualifiers

    def _get_location(self, feature):
        if type(feature.location.start) is UnknownPosition:
            spos = "?"
        else:
            spos = str(feature.location.start + 1)
        if type(feature.location.end) is UnknownPosition:
            epos = "?"
        else:
            epos = str(feature.location.end)
        return spos + ":" + epos

    def _get_pep_prop_chain_sig_id(self, feature):
        if feature.id is not None:
            return "," + feature.id
        else:
            return ""

    def _get_variant_mutagen_qualifier(self, feature, stop_codon="(", offset=1):
        """ Get x -> y or missing """
        message = feature.qualifiers["note"]
        message = message[:message.index(stop_codon)-offset] if stop_codon in message else message

        if feature.id is not None:
            message += ", " + feature.id
        return message

    def _get_accession_or_isoform(self, node):
        """ get accession in every case """
        attrs = node.attributes()
        if "isoform_accession" in attrs and attrs["isoform_accession"] is not None:
            return attrs["isoform_accession"]
        else:
            return attrs["accession"]

    def _get_position_or_isoform_position(self, node, end=False):
        """ get position, '?' if not specified"""
        attrs = node.attributes()
        val = -1
        if "isoform_position" in attrs and attrs["isoform_position"] is not None:
            val = attrs["isoform_position"]
        elif attrs["position"] is not None:
            val = attrs["position"]
        else:
            val = "?"

        if type(val) is not str and end:
            return val + len(node["aminoacid"]) - 1
        else:
            return val

    def tear_down(self):
        # We do not need to do something here, since we output it into
        # the queue which handles the tear down of opened files.
        pass

from protgraph.export.peptides.abstract_peptide_exporter import \
    APeptideExporter


class PepFasta(APeptideExporter):
    """
    This is a rather 'simple' Peptide-Fasta-File-Exporter into a single file.

    NOTE: It exports all possible paths from start to end, so before executing this
    exporter make sure that it can terminate in forseable future!
    """

    @property
    def skip_x(self) -> bool:
        return self.get_postgres_skip_x

    @property
    def peptide_min_length(self) -> int:
        return self.get_peptide_min_length

    @property
    def max_miscleavages(self) -> int:
        return self.get_miscleavages

    @property
    def use_igraph(self) -> bool:
        return self.get_use_igraph

    @property
    def peptide_max_length(self) -> int:
        return self.get_peptide_length

    @property
    def batch_size(self) -> int:
        return self.get_batch_size

    def start_up(self, **kwargs):
        # Traversal parameters:
        self.get_peptide_length = kwargs["pep_fasta_hops"]  # Number of hops. E.G. 2: s -> h_1 -> h_2 -> e
        self.get_miscleavages = kwargs["pep_fasta_miscleavages"]  # A filter criterion how many miscleavages?
        self.get_peptide_min_length = kwargs["pep_fasta_min_pep_length"]  # Peptide minimum length
        self.get_postgres_skip_x = kwargs["pep_fasta_skip_x"]
        self.get_use_igraph = kwargs["pep_fasta_use_igraph"]
        self.get_batch_size = kwargs["pep_fasta_batch_size"]

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
            # Generate header: >pg|ENTRY_ID|UNIPROT_ISO(START, END, miss:MISS; QUALIFIERS)
            entries += "|".join(
                [
                    ">pg",
                    "ID_" + str(next(self.id_gen)),
                    "".join(
                        [
                            acc, "(", str(start_pos), ":", str(end_pos), ",",
                            "mssclvg:", str(misses),
                            quali_entries,  ")"
                        ]
                    )
                ]
            )
            entries += "\n" + '\n'.join(peptide[i:i+60] for i in range(0, len(peptide), 60)) + "\n"

        # put the entries into the queue, so that the main thread can write the actual content.
        queue.put((self.output_file, entries, False))

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

            for f in qualifier:
                if f.type == "VAR_SEQ":
                    continue  # Skip since we encode this differently

                elif f.type == "VARIANT":
                    str_qualifiers.append(
                        "VARIANT[" + str(f.location.start + 1) + ":"
                        + str(f.location.end) + "," + self._get_variant_qualifier(f) + "]"
                    )

                elif f.type == "SIGNAL":
                    str_qualifiers.append(
                        "SIGNAL[" + str(f.location.start + 1) + ":" + str(f.location.end) + "]"
                    )

                elif f.type == "INIT_MET":
                    str_qualifiers.append(
                        "INIT_MET[" + str(f.location.start + 1) + ":" + str(f.location.end) + "]"
                    )

        return str_qualifiers

    def _get_variant_qualifier(self, feature):
        """ Get x -> y or missing """
        message = feature.qualifiers["note"]
        message = message[:message.index("(")-1]

        if feature.id is not None:
            message + ", " + feature.id
        return message

    def _get_accession_or_isoform(self, node):
        """ get accesion in every case """
        attrs = node.attributes()
        if "isoform_accession" in attrs and attrs["isoform_accession"] is not None:
            return attrs["isoform_accession"]
        else:
            return attrs["accession"]

    def _get_position_or_isoform_position(self, node, end=False):
        """ get position, '?' if not specified"""
        attrs = node.attributes()
        val = -1
        if "isoform_accession" in attrs and attrs["isoform_accession"] is not None:
            val = attrs["isoform_position"]
        elif attrs["position"] is not None:
            val = attrs["position"]
        else:
            val = "?"

        if type(val) != str and end:
            return val + len(node["aminoacid"]) - 1
        else:
            return val

    def tear_down(self):
        # We do not need to do something here, since we output it into
        # the queue which handles the tear down of opened files.
        pass

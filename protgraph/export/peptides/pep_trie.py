import os

from protgraph.export.peptides.pep_fasta import PepFasta


class PepTrie(PepFasta):
    """
    This is a trie exporter for peptides (which directly saves them on the filesystem).

    NOTE: It is recommended to put the export folder in a filesystem, which can handle large number
    of folders and files. XFS would be one option.
    """

    def start_up(self, **kwargs):
        super(PepTrie, self).start_up(**kwargs)

        self.id_gen = self.unique_id_gen(**kwargs)

        # Get output file from configuration
        self.output_folder = kwargs["pep_trie_folder_out"]

    def export_peptides(self, prot_graph, l_path_nodes, l_path_edges, l_peptide, l_miscleavages, queue):
        # Export a batch of peptides at onces
        # Build up the entries for the batch
        for peptide, nodes, edges, misses in zip(l_peptide, l_path_nodes, l_path_edges, l_miscleavages):
            # Get Peptide information
            acc = self._get_accession_or_isoform(prot_graph.vs[nodes[1]])
            start_pos = self._get_position_or_isoform_position(prot_graph.vs[nodes[1]])
            end_pos = self._get_position_or_isoform_position(prot_graph.vs[nodes[-2]], end=True)
            l_str_qualifiers = self._get_qualifiers(prot_graph, edges)
            quali_entries = ",".join(l_str_qualifiers)
            if quali_entries:
                quali_entries = "," + quali_entries

            # Generate only specific entry ,UNIPROT_ISO(START:END, mssclvg:MISS; QUALIFIERS)
            entry = "," + "".join([
                acc, "(", str(start_pos), ":", str(end_pos), ",",
                "mssclvg:", str(misses), quali_entries,  ")"
            ])

            # Create Trie on Filesystem (using OS)
            p = os.path.join(self.output_folder, *[x for x in peptide[:]])
            os.makedirs(p, exist_ok=True)

            # Write actual Peptide
            queue.put((
                os.path.join(p, ".proteins"),
                entry, False, "ac"
            ))

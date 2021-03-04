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

        # Get output file from configuration
        self.output_file = kwargs["pep_fasta_out"]

    def export_peptides(self, prot_graph, l_path_nodes, l_path_edges, l_peptide, l_miscleavages, queue):
        # Export a batch of peptides at once
        accession = prot_graph.vs[0]["accession"]

        # Build up the entries for the batch
        entries = ""
        if "qualifiers" in prot_graph.es[0].attributes():
            for peptide, edges in zip(l_peptide, l_path_edges):
                qualifiers = prot_graph.es[edges]["qualifiers"]
                # replace empty lists with None and get type
                qualifiers = [None if x is None or len(x) == 0 else [y.type for y in x] for x in qualifiers]

                entries += (
                    ">lcl|ACCESSION=" + accession +
                    "|QUALIFIERS=" + ",".join(";".join(x) if x is not None else "" for x in qualifiers)
                )
                entries += "\n" + '\n'.join(peptide[i:i+60] for i in range(0, len(peptide), 60)) + "\n"
        else:
            for peptide, edges in zip(l_peptide, l_path_edges):
                entries += ">lcl|ACCESSION=" + accession + "|QUALIFIERS="
                entries += "\n" + '\n'.join(peptide[i:i+60] for i in range(0, len(peptide), 60)) + "\n"

        # put the entries into the queue, so that the main thread can write the actual content.
        queue.put((self.output_file, entries))

    def tear_down(self):
        # We do not need to do something here, since we output it into
        # the queue which handles the tear down of opened files.
        pass

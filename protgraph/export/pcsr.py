import csv
import os

from protgraph.export.abstract_exporter import AExporter


class PCsr(AExporter):
    """ A simple CSV Exporter. This export is compatible with Gephi """

    def start_up(self, **kwargs):
        # Here we simply create the folder if it does not exist
        self.out_folder = kwargs["export_output_folder"]
        os.makedirs(self.out_folder, exist_ok=True)

        self.flat = not kwargs["export_in_directories"]

    def export(self, prot_graph, _):
        if self.flat:
            accession = prot_graph.vs["accession"][0]
            self._write_csr_attrs(os.path.join(self.out_folder, accession + ".pcsv"), prot_graph)
        else:
            accession = prot_graph.vs["accession"][0]
            out_dir = os.path.join(self.out_folder, *[x for x in accession[:-1]])
            # Create outfolders if needed
            os.makedirs(out_dir, exist_ok=True)
            self._write_csr_attrs(os.path.join(os.path.join(out_dir, accession[-1:] + ".pcsv")), prot_graph)

    def _write_csr_attrs(self, out_file, graph):

        AC = graph.vs[0]["accession"]
        IA = None  # Get List of Isos (Accessions)
        NO = None  # Get List of Nodes
        ED = None  # Get List of Edges
        SQ = None  # Get List of Node[Sequence]
        IS = None  # Get List of Node[Isos] substituted
        MW = None  # Get List of Node[Mono Weight]
        AW = None  # Get List of Node[Average Weight]
        TO = None  # The Topological Order
        CL = None  # Get List of Edges[Cleaved]
        QU = None  # Get List of Edges[Qualifiers] (simplified as in FASTA?)
        VC = None  # Get List of Edges[Var-Count]
        PD = None  # List of Lists of PDB-Entries



        with open(out_file, "w") as csvfile:
            writer = csv.writer(csvfile)
            # Write Header
            header = ["Id", *graph.vs.attributes()]
            writer.writerow(header)

            # Write "Body"
            for n in graph.vs:
                writer.writerow([n.index, *n.attributes().values()])

    def tear_down(self):
        # We do not need to tear down a graph export to dot files
        pass

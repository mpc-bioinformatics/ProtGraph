import csv
import os

from protgraph.export.abstract_exporter import AExporter


class CSV(AExporter):
    """ A simple CSV Exporter. This export is compatible with Gephi """

    def start_up(self, **kwargs):
        # Here we simply create the folder if it does not exist
        self.out_folder = kwargs["export_output_folder"]
        os.makedirs(self.out_folder, exist_ok=True)

        self.flat = not kwargs["export_in_directories"]

    def export(self, prot_graph, _):
        if self.flat:
            accession = prot_graph.vs["accession"][0]
            self._write_nodes_data(os.path.join(self.out_folder, accession + "_nodes.csv"), prot_graph)
            self._write_edges_data(os.path.join(self.out_folder, accession + "_edges.csv"), prot_graph)
        else:
            accession = prot_graph.vs["accession"][0]
            out_dir = os.path.join(self.out_folder, *[x for x in accession[:-1]])
            # Create outfolders if needed
            os.makedirs(out_dir, exist_ok=True)

            self._write_nodes_data(os.path.join(os.path.join(out_dir, accession[-1:] + "_nodes.csv")), prot_graph)
            self._write_edges_data(os.path.join(os.path.join(out_dir, accession[-1:] + "_edges.csv")), prot_graph)

    def _write_nodes_data(self, out_file, graph):
        with open(out_file, "w") as csvfile:
            writer = csv.writer(csvfile)
            # Write Header
            header = ["Id", *graph.vs.attributes()]
            writer.writerow(header)

            # Write "Body"
            for n in graph.vs:
                writer.writerow([n.index, *n.attributes().values()])

    def _write_edges_data(self, out_file, graph):
        with open(out_file, "w") as csvfile:
            writer = csv.writer(csvfile)
            # Write Header
            header = ["Source", "Target", *graph.es.attributes()]
            writer.writerow(header)

            # Write "Body"
            for e in graph.es:
                writer.writerow([e.source, e.target, *e.attributes().values()])

    def tear_down(self):
        # We do not need to tear down a graph export to dot files
        pass

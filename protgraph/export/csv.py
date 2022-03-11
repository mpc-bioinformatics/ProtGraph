import csv

from protgraph.export.generic_file_exporter import GenericFileExporter


class CSV(GenericFileExporter):
    """ A simple CSV exporter. This export is compatible with Gephi."""

    def __init__(self):
        super(CSV, self).__init__(
            self._lambda_export
        )

    def _lambda_export(self, pg, path):
        self.write_nodes_data(path + "_nodes.csv", pg)
        self.write_edges_data(path + "_edges.csv", pg)

    def write_nodes_data(self, out_file, graph):
        with open(out_file, "w") as csvfile:
            writer = csv.writer(csvfile)
            # Write Header
            header = ["Id", *graph.vs.attributes()]
            writer.writerow(header)

            # Write "Body"
            for n in graph.vs:
                writer.writerow([n.index, *n.attributes().values()])

    def write_edges_data(self, out_file, graph):
        with open(out_file, "w") as csvfile:
            writer = csv.writer(csvfile)
            # Write Header
            header = ["Source", "Target", *graph.es.attributes()]
            writer.writerow(header)

            # Write "Body"
            for e in graph.es:
                writer.writerow([e.source, e.target, *e.attributes().values()])

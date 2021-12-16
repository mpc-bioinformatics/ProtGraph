import csv
import io
import os

from protgraph.export.abstract_exporter import AExporter


class LargeCSV(AExporter):
    """ A CSV Exporter to export to one/two large files. This export is compatible with Gephi """

    def start_up(self, **kwargs):
        # Here we simply set the folder
        self.out_folder = kwargs["export_output_folder"]
        self.nodes_out = os.path.join(self.out_folder, "nodes.csv")
        self.edges_out = os.path.join(self.out_folder, "edges.csv")
        self.id_get = self.unique_id_gen(**kwargs)

    def export(self, prot_graph, queue):
        # Write Header Once!
        n_header = [
            "Id", "accession", "aminoacid", "position", "isoform_accession", "isoform_position",
            "mono_weight", "mono_weight_to_end", "avrg_weight", "avrg_weight_to_end"
        ]
        e_header = [
            "Source", "Target", "cleaved", "qualifiers"
        ]

        queue.put((self.nodes_out, ",".join(n_header) + "\n", True, "a"))
        queue.put((self.edges_out, ",".join(e_header) + "\n", True, "a"))

        # Get Id-Mapper
        node_mapping = [next(self.id_get) for _ in range(prot_graph.vcount())]

        # Write Nodes and edges
        self._write_nodes_data(node_mapping, n_header[1:], queue, prot_graph)
        self._write_edges_data(node_mapping, e_header[2:], queue, prot_graph)

    def _write_nodes_data(self, node_mapping, attrs, queue, graph):
        # Build up string
        str_out = io.StringIO()
        writer = csv.writer(str_out)
        for n in graph.vs:
            writer.writerow(
                [node_mapping[n.index], *self.getattrs(n, attrs)]
            )

        # Add to queue
        queue.put(
            (
                self.nodes_out,
                str_out.getvalue(),
                False, "a"
            )
        )

    def _write_edges_data(self, node_mapping, attrs, queue, graph):
        # Build up string
        str_out = io.StringIO()
        writer = csv.writer(str_out)
        for e in graph.es:
            writer.writerow(
                [node_mapping[e.source], node_mapping[e.target], *self.getattrs(e, attrs)]
            )

        # Add to queue
        queue.put(
            (
                self.edges_out,
                str_out.getvalue(),
                False, "a"
            )
        )

    def getattrs(self, node_edge, attrs):
        neattr = node_edge.attributes()
        return [
            neattr[a] if a in neattr else None
            for a in attrs
        ]

    def tear_down(self):
        # We do not need to tear down a graph export to dot files
        pass

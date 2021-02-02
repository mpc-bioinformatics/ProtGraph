import os

from export.abstract_exporter import AExporter


class GraphML(AExporter):
    """
    A simple GraphML (XML-format) exporter. This one is
    recommended since it can save all attributes used in the graph!
    """

    def start_up(self, **kwargs):
        # Here wer simply create the folder if it does not exist
        self.out_folder = kwargs["export_output_folder"]
        os.makedirs(self.out_folder, exist_ok=True)

        self.flat = not kwargs["export_in_directories"]

    def export(self, prot_graph):
        if self.flat:
            accession = prot_graph.vs["accession"][0]
            prot_graph.write_graphml(os.path.join(self.out_folder, accession + ".graphml"))
        else:
            accession = prot_graph.vs["accession"][0]
            out_dir = os.path.join(self.out_folder, *[x for x in accession[:-1]])
            # Create outfolders if needed
            os.makedirs(out_dir, exist_ok=True)

            prot_graph.write_graphml(os.path.join(out_dir, accession[-1:] + ".graphml"))

    def tear_down(self):
        # We do not need to tear down a graph export to graphml files
        pass

import os

from protgraph.export.abstract_exporter import AExporter


class Pickle(AExporter):
    """ A simple Pickle Exporter using igraph """

    def start_up(self, **kwargs):
        # Here we simply create the folder if it does not exist
        self.out_folder = kwargs["export_output_folder"]
        os.makedirs(self.out_folder, exist_ok=True)

        self.flat = not kwargs["export_in_directories"]

    def export(self, prot_graph, _):
        if self.flat:
            accession = prot_graph.vs["accession"][0]
            prot_graph.write_pickle(os.path.join(self.out_folder, accession + ".pickle"))
        else:
            accession = prot_graph.vs["accession"][0]
            out_dir = os.path.join(self.out_folder, *[x for x in accession[:-1]])
            # Create outfolders if needed
            os.makedirs(out_dir, exist_ok=True)

            prot_graph.write_pickle(os.path.join(out_dir, accession[-1:] + ".pickle"))

    def tear_down(self):
        # We do not need to tear down a graph export to dot files
        pass

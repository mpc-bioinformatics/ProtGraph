import os

from export.abstract_exporter import AExporter


class Pickle(AExporter):
    """ A simple Pickle Exporter using igraph """

    def start_up(self, **kwargs):
        # Here we simply create the folder if it does not exist
        self.out_folder = kwargs["export_output_folder"]
        if not os.path.exists(self.out_folder):
            os.makedirs(self.out_folder)

    def export(self, prot_graph):
        accession = prot_graph.vs["accession"][0]
        out_dir = os.path.join(self.out_folder, *[x for x in accession[:-1]])
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        prot_graph.write_pickle(os.path.join(out_dir, accession[-1:] + ".pickle"))

    def tear_down(self):
        # We do not need to tear down a graph export to dot files
        pass

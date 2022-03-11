import os

from protgraph.export.abstract_exporter import AExporter


class GenericFileExporter(AExporter):
    """ A simple Generic File Exporter """

    def __init__(self, export_function):
        """
        Insert Export-Function as Lambda, containing the following inputs:
        Input: (Graph, Output-Folder)
        Output: <None>
        """
        self.export_function = export_function

    def start_up(self, **kwargs):
        # Here we simply create the folder if it does not exist
        self.out_folder = kwargs["export_output_folder"]
        os.makedirs(self.out_folder, exist_ok=True)

        self.flat = not kwargs["export_in_directories"]

    def export(self, prot_graph, _):
        if self.flat:
            accession = prot_graph.vs["accession"][0]
            self.export_function(prot_graph, os.path.join(self.out_folder, accession))
        else:
            accession = prot_graph.vs["accession"][0]
            out_dir = os.path.join(self.out_folder, *[x for x in accession[:-1]])
            # Create outfolders if needed
            os.makedirs(out_dir, exist_ok=True)

            self.export_function(prot_graph, os.path.join(out_dir, accession[-1:]))

    def tear_down(self):
        # We do not need to tear down a file export
        pass

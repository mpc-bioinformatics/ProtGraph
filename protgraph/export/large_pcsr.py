import os

from protgraph.export.abstract_exporter import AExporter
from protgraph.export.pcsr import PCsr

class LargePCsr(AExporter):
    """ A simple Protein Compressed Sparse Row Exporter, exporting to a single file """

    def start_up(self, **kwargs):
        self.pcsr_exporter = PCsr()
        self.pcsr_exporter.pdb_count = kwargs["export_large_pcsr_pdb_entries"]
        self.out_folder = kwargs["export_output_folder"]
        self.out_file = os.path.join(self.out_folder, "database.pcsv")


    def export(self, prot_graph, queue):
        queue.put((
            self.out_file, 
            self.pcsr_exporter._build_csr_entry(prot_graph),
            False
        ))

    def tear_down(self):
        # We do not need to tear down a simple file export
        pass

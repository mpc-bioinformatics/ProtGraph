import os

from protgraph.export.abstract_exporter import AExporter
from protgraph.export.pcsr import PCSR


class LargePCSR(AExporter):
    """ A simple Protein Compressed Sparse Row Exporter, exporting to a single file """

    def start_up(self, **kwargs):
        self.pcsr_exporter = PCSR()
        self.pcsr_exporter.start_up(**kwargs)
        self.pcsr_exporter.pdb_count = kwargs["export_large_pcsr_pdb_entries"]
        self.out_folder = kwargs["export_output_folder"]
        self.out_file = os.path.join(self.out_folder, "database.pcsr")

    def export(self, prot_graph, queue):
        queue.put((
            self.out_file,
            self.pcsr_exporter._build_csr_string(
                self.pcsr_exporter._build_csr_entry(prot_graph)
            ),
            False, "a"
        ))

    def tear_down(self):
        # We do not need to tear down a simple file export
        pass

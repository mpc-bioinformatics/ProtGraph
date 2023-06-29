import os

from protgraph.export.abstract_exporter import AExporter
from protgraph.export.binary_pcsr import BinaryPCSR


class LargeBinaryPCSR(AExporter):
    """ A simple Protein Compressed Sparse Row Exporter, exporting to a single file """

    def start_up(self, **kwargs):
        self.bpcsr_exporter = BinaryPCSR()
        self.bpcsr_exporter.start_up(**kwargs)
        self.bpcsr_exporter.pdb_count = kwargs["export_large_binary_pcsr_pdb_entries"]
        self.bpcsr_exporter.ctype_mapping = self.bpcsr_exporter._mapping()
        self.out_folder = kwargs["export_output_folder"]
        self.out_file = os.path.join(self.out_folder, "database.bpcsr")

    def export(self, prot_graph, queue):
        queue.put((
            self.out_file,
            self.bpcsr_exporter._pcsr_to_bytes(
                self.bpcsr_exporter._build_csr_entry(prot_graph)
            ),
            False, "ab"
        ))

    def tear_down(self):
        # We do not need to tear down a simple file export
        pass

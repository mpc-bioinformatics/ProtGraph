import numpy as np

from protgraph.export.pcsr import PCSR


class BinaryPCSR(PCSR):
    """ A simple Protein Compressed Sparse Row Exporter, exporting to a single file """

    def _mapping(self):
        return dict(
            AC=str,
            CN=32,  # 4 Numbers always!
            NO=32,
            ED=32,
            SQ=str,
            PO=16,
            IS=8,
            IP=16,
            MW=64,
            CL=bool,
            QU=str,
            VC=8,
            PD=64,
        )

    def start_up(self, **kwargs):
        self.ctype_mapping = self._mapping()
        super(BinaryPCSR, self).start_up(**kwargs)
        self.pdb_count = kwargs["export_binary_pcsr_pdb_entries"]

    def write_pcsr(self, pg, path):
        out_list = self._build_csr_entry(pg)
        out_bytes = self._pcsr_to_bytes(out_list)

        with open(path + ".bpcsr", "wb") as out:
            out.write(out_bytes)

    def _pcsr_to_bytes(self, b_list):
        b = bytearray()
        for key, values in b_list:
            if self.ctype_mapping[key] == str:
                b.extend(b"".join([x.encode() + bytes(1) for x in values]))

            elif key == "PD":
                uint_bytes = b"".join([
                    (int(k)).to_bytes(int(self.ctype_mapping[key]/8), byteorder="big", signed=False)
                    if k is not np.nan
                    else b"\xff"*int(self.ctype_mapping[key]/8)
                    for x in values for y in x for k in y
                ])
                b.extend(uint_bytes)

            elif self.ctype_mapping[key] == bool:
                hot_encoded = [255 if x == "t" else 0 for x in values]
                uint_bytes = b"".join([
                    (x).to_bytes(1, byteorder="big", signed=False)
                    for x in hot_encoded
                ])
                b.extend(uint_bytes)

            elif type(self.ctype_mapping[key]) is int:
                uint_bytes = b"".join([
                    (x).to_bytes(int(self.ctype_mapping[key]/8), byteorder="big", signed=False)
                    if x >= 0
                    else (2**self.ctype_mapping[key]-1).to_bytes(
                        int(self.ctype_mapping[key]/8), byteorder="big", signed=False
                    )
                    for x in values
                ])
                b.extend(uint_bytes)

        return b

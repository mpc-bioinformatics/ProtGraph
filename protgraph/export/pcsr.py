import numpy as np

from protgraph.export.generic_file_exporter import GenericFileExporter
from protgraph.export.peptides.pep_fasta import PepFasta
from protgraph.graph_statistics import _count_feature


class PCSR(GenericFileExporter):
    """ A simple Protein Compressed Sparse Row  Exporter """

    def __init__(self):
        super(PCSR, self).__init__(
            self.write_pcsr
        )
        self._get_qualifier = PepFasta()._map_qualifier_to_string

    def start_up(self, **kwargs):
        self.pdb_count = kwargs["export_pcsr_pdb_entries"]

        if kwargs["pcsr_feature_to_count"] is None:
            self.features_to_check = ["VARIANT"]
        else:
            self.features_to_check = []
            for i in kwargs["pcsr_feature_to_count"]:
                self.features_to_check.append(i)

        super(PCSR, self).start_up(**kwargs)

    def write_pcsr(self, pg, path):
        out_str = self._build_csr_string(self._build_csr_entry(pg))
        with open(path + ".pcsr", "w") as out:
            out.write(out_str)

    def _build_csr_entry(self, graph):
        # The CSR is reordered in Top Order in favor on going through the RAM sequentially. To be tested...

        # get the topological order
        TO = self.__get_protein_graph_specific_top_order(graph)  # The Topological Order
        TO_EDGES = []

        # get Accessions (including isoforms)
        AC = graph.vs[0]["accession"]
        IA = \
            list(set(graph.vs["isoform_accession"]).difference(set([None]))) \
            if "isoform_accession" in graph.vs[0].attributes() \
            else []  # Get List of Isos (Accessions)

        NO = []  # Get List of Nodes
        ED = []  # Get List of Edges
        TO_EDGES = []

        out_edges = graph.vs[TO[0]].out_edges()
        NO.append(len(out_edges))
        sorted_edges = sorted([(TO.index(e.target), e.index) for e in out_edges], key=lambda x: x[0])
        ED.extend(x[0] for x in sorted_edges)
        TO_EDGES.extend(x[1] for x in sorted_edges)
        for n_idx in TO[1:]:
            out_edges = graph.vs[n_idx].out_edges()
            NO.append(len(out_edges) + NO[-1])
            sorted_edges = sorted([(TO.index(e.target), e.index) for e in out_edges], key=lambda x: x[0])
            ED.extend(x[0] for x in sorted_edges)
            TO_EDGES.extend(x[1] for x in sorted_edges)

        CN = [len([AC, *IA]), graph.vcount(), graph.ecount(), self.pdb_count]  # Every Count
        # All these attributes need to be reordered if NO or ED gets modified!

        # Get List of Node[Sequence]
        SQ = graph.vs[TO]["aminoacid"]

        # Get List of Node[Isos] substituted (idx + 1 due to AC)
        IS = [
            IA.index(x) + 1 if x is not None else 0
            for x in graph.vs[TO]["isoform_accession"]
        ] if "isoform_accession" in graph.vs[0].attributes() else [0]*len(NO)

        # Get List of Node[Mono Weight]
        MW = graph.vs[TO]["mono_weight"] if "mono_weight" in graph.vs[0].attributes() else [0]*len(NO)

        # Add Position-Information
        PO = [x if x is not None else -1 for x in graph.vs[TO]["position"]]
        if "isoform_position" in graph.vs[0].attributes():
            IP = [x if x is not None else -1 for x in graph.vs[TO]["isoform_position"]]
        else:
            IP = [-1]*len(NO)

        # Attribute for Edges
        CL = ["t" if x else "f" for x in graph.es[TO_EDGES]["cleaved"]]  # Get List of Edges[Cleaved]

        if "qualifiers" in graph.es[0].attributes():
            VC = [
                sum([_count_feature(x, ftc, min) for ftc in self.features_to_check])
                for x in graph.es[TO_EDGES]["qualifiers"]
            ]
            QU = []  # Get List of Edges[Qualifiers] (simplified as in FASTA)
            for qualifier in graph.es[TO_EDGES]["qualifiers"]:
                if qualifier is None or len(qualifier) == 0:
                    QU.append([])
                else:
                    QU.append(self._get_qualifier(qualifier))
            QU = ["&".join(x) for x in QU]

        else:
            VC = [0]*len(ED)
            QU = [""]*len(ED)

        # Generate additional information like PDBs
        if "mono_weight" in graph.vs[0].attributes():
            # Build PDB
            self.__build_pdb(graph, k=self.pdb_count)
            # Fill entries with nans
            pdb_entries = [x + [[np.nan, np.nan]]*(self.pdb_count - len(x)) for x in graph.vs["pdb"]]

            # Delete entries in graph itself
            del graph.vs["pdb"]

            PD = pdb_entries  # List of Lists of PDB-Entries
            PD = [PD[x] for x in TO]

        else:
            PD = [[[0, 0]]*self.pdb_count]*len(NO)

        # Order to write:
        build_list = [
            ("CN", CN),
            ("AC", [AC, *IA]),
            ("NO", NO),
            ("ED", ED),
            ("SQ", SQ),
            ("PO", PO),
            ("IS", IS),
            ("IP", IP),
            ("MW", MW),
            ("CL", CL),
            ("QU", QU),
            ("VC", VC),
            ("PD", PD),
        ]

        return build_list

    def _build_csr_string(self, build_list):
        build_str = [
            (build_list[0][0]+"   " + ";".join([str(x) for x in build_list[0][1]])),
            (build_list[1][0]+"   " + ";".join(build_list[1][1])),
            (build_list[2][0]+"   " + ";".join([str(x) for x in build_list[2][1]])),
            (build_list[3][0]+"   " + ";".join([str(x) for x in build_list[3][1]])),
            (build_list[4][0]+"   " + ";".join(build_list[4][1])),
            (build_list[5][0]+"   " + ";".join([str(x) for x in build_list[5][1]])),
            (build_list[6][0]+"   " + ";".join([str(x) for x in build_list[6][1]])),
            (build_list[7][0]+"   " + ";".join([str(x) for x in build_list[7][1]])),
            (build_list[8][0]+"   " + ";".join([str(x) for x in build_list[8][1]])),
            (build_list[9][0]+"   " + ";".join(build_list[9][1])),
            (build_list[10][0]+"   " + ";".join(build_list[10][1])),
            (build_list[11][0]+"   " + ";".join([str(x) for x in build_list[11][1]])),
            (build_list[12][0]+"   " + ";".join(["^".join(["^".join(str(z) for z in y) for y in x]) for x in build_list[12][1]])), # noqa E501
        ]

        return "\n".join(build_str) + "\n\n"

    def __get_protein_graph_specific_top_order(self, _graph):
        """
        TODO Explanation!
        """
        sorted_by_position_attr = []
        s = set([x for x, y in zip(range(_graph.vcount()), _graph.vs.indegree()) if y == 0])
        marked_edges = [False]*_graph.ecount()

        while len(s) != 0:
            t = []  # (isoform_name, iso_pos, pos, n)
            for x in s:
                node_attrs = _graph.vs[x].attributes()

                if "isoform_accession" in node_attrs:
                    t1 = node_attrs["isoform_accession"] if node_attrs["isoform_accession"] else node_attrs["accession"]
                else:
                    t1 = node_attrs["accession"]

                if "isoform_position" in node_attrs:
                    t2 = node_attrs["isoform_position"] if node_attrs["isoform_position"] else float("-inf")
                else:
                    t2 = float("-inf")

                t3 = node_attrs["position"] if node_attrs["position"] else float("-inf")
                t.append((t1, t2, t3, x))
            # sorted up down down
            res = sorted(t, key=lambda x: (-len(x[0]), [-ord(c) for c in x[0]], x[1], x[2]))
            n = res[0][3]
            s.remove(n)

            sorted_by_position_attr.append(n)
            for e_out in _graph.vs[n].out_edges():
                marked_edges[e_out.index] = True
                in_out_edges = [x.index for x in _graph.vs[e_out.target].in_edges()]
                if all([marked_edges[x] for x in in_out_edges]):
                    s.add(e_out.target)

        return sorted_by_position_attr

    def __shift_interval_by(self, intervals, weight):
        """ Shift the intervals by weight """
        return [[x + weight, y + weight] for [x, y] in intervals]

    def __merge_overlapping_intervals(self, intervals):
        """ Get overlapping intervals and merge them """
        intervals = np.array(intervals)
        starts = intervals[:, 0]
        ends = np.maximum.accumulate(intervals[:, 1])
        valid = np.zeros(len(intervals) + 1, dtype=bool)
        valid[0] = True
        valid[-1] = True
        valid[1:-1] = starts[1:] >= ends[:-1]
        return [list(x) for x in np.vstack((starts[:][valid[:-1]], ends[:][valid[1:]])).T]

    def __merge_closest_intervals(self, intervals):
        """ Get the closest interval and merge those two """
        diff = [y[0]-x[1] for x, y in zip(intervals, intervals[1:])]
        argmin = diff.index(min(diff))

        new_interval = [intervals[argmin][0], intervals[argmin+1][1]]
        return intervals[:argmin] + [new_interval] + intervals[argmin+2:]

    def __build_pdb(self, graph, k=5):
        """
        Generates the pdb (intervals). Each node will have up to k many intervals
        We build it up via the rev. top. sort. Intervals are merged they overlap.
        If too many intervals are present the closest ones will be merged.

        This uses the mono_weight ONLY, but could be extended to use the avrg_weight (TODO DL?)
        """
        rev_top_sort = graph.topological_sorting(mode="IN")

        # Initial attribute values:
        graph.vs[rev_top_sort[0]]["pdb"] = [[0, 0]]

        # iterate
        for node in rev_top_sort[1:]:
            intervals = []
            for out_edge in graph.incident(node, mode="OUT"):
                intervals.extend(
                    self.__shift_interval_by(
                        graph.vs[graph.es[out_edge].target]["pdb"],
                        graph.vs[graph.es[out_edge].target]["mono_weight"]
                    )
                )

            sorted_intervals = self.__merge_overlapping_intervals(sorted(intervals, key=lambda x: x[0]))

            while True:
                if len(sorted_intervals) <= k:
                    break
                else:
                    sorted_intervals = self.__merge_closest_intervals(sorted_intervals)

            graph.vs[node]["pdb"] = sorted_intervals

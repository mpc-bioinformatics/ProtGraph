import os
import unittest

import protgraph


class FunctionalTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """ Set base Example Folder """
        main_file_path = os.path.dirname(os.path.abspath(protgraph.__file__))
        cls.examples_path = os.path.join(main_file_path, "..", "examples")
        cls.example_files = [
            os.path.abspath(os.path.join(cls.examples_path, "e_coli.dat")),
            os.path.abspath(os.path.join(cls.examples_path, "p53_human.txt"))
        ]
        cls.procs_num = ["-n", "1"]

    def test_minimal(self):
        args = protgraph.parse_args([] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_none(self):
        args = protgraph.parse_args(["-ft", "NoNE"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_all(self):
        args = protgraph.parse_args(["-ft", "ALl"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_isoforms(self):
        args = protgraph.parse_args(["-ft", "VAR_SeQ"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_variants(self):
        args = protgraph.parse_args(["-ft", "VARIAnT"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_met(self):
        args = protgraph.parse_args(["-ft", "IniT_MET"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_signal(self):
        args = protgraph.parse_args(["-ft", "SIGnaL"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_digestion_skip(self):
        args = protgraph.parse_args(["-d", "skip"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_digestion_trypsin(self):
        args = protgraph.parse_args(["-d", "trypsin"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_digestion_full(self):
        args = protgraph.parse_args(["-d", "full"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_no_merge(self):
        args = protgraph.parse_args(["-nm"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_annotate_weights(self):
        args = protgraph.parse_args(["-aawe", "-amwe"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_replacement(self):
        args = protgraph.parse_args(["-raa", "A->b,C,d"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_statistics_possibilites(self):
        args = protgraph.parse_args(["-cnp"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_statistics_miscleavages(self):
        args = protgraph.parse_args(["-cnpm"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_statistics_hops(self):
        args = protgraph.parse_args(["-cnph"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_export_dot(self):
        args = protgraph.parse_args(["-edot"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_export_csv(self):
        args = protgraph.parse_args(["-ecsv"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_export_lcsv(self):
        args = protgraph.parse_args(["-elcsv"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_export_graphml(self):
        args = protgraph.parse_args(["-egraphml"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_export_gml(self):
        args = protgraph.parse_args(["-epickle"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_export_pickle(self):
        args = protgraph.parse_args(["-egml"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_export_pep_fasta(self):
        args = protgraph.parse_args(["-epepfasta", "--pep_fasta_hops", "2"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_issue8(self):
        args = protgraph.parse_args(["-n", "1", os.path.join(self.examples_path, "Q9QXS1.txt")])
        protgraph.prot_graph(**args)

    def test_issue13(self):
        args = protgraph.parse_args(["-n", "1", os.path.join(self.examples_path, "F1SN05.txt")])
        protgraph.prot_graph(**args)

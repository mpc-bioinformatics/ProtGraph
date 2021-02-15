import os
import unittest

import prot_graph


class FunctionalTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """ Set base Example Folder """
        main_file_path = os.path.dirname(os.path.abspath(prot_graph.__file__))
        cls.examples_path = os.path.join(main_file_path, "examples")
        cls.example_files = [
            os.path.join(cls.examples_path, "e_coli.dat"),
            os.path.join(cls.examples_path, "p53_human.txt")
        ]

    def test_minimal(self):
        args = prot_graph.parse_args([] + self.example_files)
        prot_graph.main(args)

    def test_skip_isoforms(self):
        args = prot_graph.parse_args(["-si"] + self.example_files)
        prot_graph.main(args)

    def test_skip_variants(self):
        args = prot_graph.parse_args(["-sv"] + self.example_files)
        prot_graph.main(args)

    def test_skip_met(self):
        args = prot_graph.parse_args(["-sm"] + self.example_files)
        prot_graph.main(args)

    def test_skip_signal(self):
        args = prot_graph.parse_args(["-ss"] + self.example_files)
        prot_graph.main(args)

    def test_digestion_skip(self):
        args = prot_graph.parse_args(["-d", "skip"] + self.example_files)
        prot_graph.main(args)

    def test_digestion_trypsin(self):
        args = prot_graph.parse_args(["-d", "trypsin"] + self.example_files)
        prot_graph.main(args)

    def test_digestion_full(self):
        args = prot_graph.parse_args(["-d", "full"] + self.example_files)
        prot_graph.main(args)

    def test_no_merge(self):
        args = prot_graph.parse_args(["-nm"] + self.example_files)
        prot_graph.main(args)

    def test_annotate_weights(self):
        args = prot_graph.parse_args(["-aawe", "-amwe"] + self.example_files)
        prot_graph.main(args)

    def test_statistics_possibilites(self):
        args = prot_graph.parse_args(["-cnp"] + self.example_files)
        prot_graph.main(args)

    def test_statistics_miscleavages(self):
        args = prot_graph.parse_args(["-cnpm"] + self.example_files)
        prot_graph.main(args)

    def test_statistics_hops(self):
        args = prot_graph.parse_args(["-cnph"] + self.example_files)
        prot_graph.main(args)

    def test_export_dot(self):
        args = prot_graph.parse_args(["-edot"] + self.example_files)
        prot_graph.main(args)

    def test_export_graphml(self):
        args = prot_graph.parse_args(["-egraphml"] + self.example_files)
        prot_graph.main(args)

    def test_export_gml(self):
        args = prot_graph.parse_args(["-epickle"] + self.example_files)
        prot_graph.main(args)

    def test_export_pickel(self):
        args = prot_graph.parse_args(["-egml"] + self.example_files)
        prot_graph.main(args)

    def test_issue8(self):
        args = prot_graph.parse_args(["-n", "1", os.path.join(self.examples_path, "Q9QXS1.txt")])
        prot_graph.main(args)

    def test_issue13(self):
        args = prot_graph.parse_args(["-n", "1", os.path.join(self.examples_path, "F1SN05.txt")])
        prot_graph.main(args)

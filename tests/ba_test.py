import filecmp
import os
import shutil
import unittest

import protgraph

from tests.test_helper import get_file_from_name


class FunctionalTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        base_dir = os.path.dirname(__file__)
        cls.input_dir = os.path.join(base_dir, "test_data", "input")
        cls.expected_dir = os.path.join(base_dir, "test_data", "expected")
        cls.output_dir = os.path.join(base_dir, "temp")
        cls.procs_num = ["-n", "1"]

        cls.test_cases = [
            {
                "name": "Canonical",
                "input_files": ["JKBAAA.txt"],
                "expected": "expected_JKBAAA.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Canonical-Merge",
                "input_files": ["JKBAAB.txt"],
                "expected": "expected_JKBAAB.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Variant",
                "input_files": ["JKBAAC.txt"],
                "expected": "expected_JKBAAC.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Mutation",
                "input_files": ["JKBAAD.txt"],
                "expected": "expected_JKBAAD.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Conflict",
                "input_files": ["JKBAAE.txt"],
                "expected": "expected_JKBAAE.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Variant-Missing",
                "input_files": ["JKBAAF.txt"],
                "expected": "expected_JKBAAF.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Isoform",
                "input_files": ["JKBAIA.txt"],
                "expected": "expected_JKBAIA.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Isoform-Canonical-Merge",
                "input_files": ["JKBAIB.txt"],
                "expected": "expected_JKBAIB.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Shared-VARSEQ",
                "input_files": ["JKBAIC.txt"],
                "expected": "expected_JKBAIC.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Isoform-into-Variant",
                "input_files": ["JKBAID.txt"],
                "expected": "expected_JKBAID.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Partially-Shared-VARSEQs",
                "input_files": ["JKBAIE.txt"],
                "expected": "expected_JKBAIE.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Isoform-and-Variant-Simultaneous",
                "input_files": ["JKBAIF.txt"],
                "expected": "expected_JKBAIF.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Isoform-Missing",
                "input_files": ["JKBAIG.txt"],
                "expected": "expected_JKBAIG.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Canonical-Peptide",
                "input_files": ["JKBAAA.txt"],
                "output_file": "JKBAAA_peptide.graphml",
                "expected": "expected_JKBAAA_peptide.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pep", "A", "-of", "JKBAAA_peptide"],
            },
            {
                "name": "Canonical-Merge-Peptide-Single-Cut",
                "input_files": ["JKBAAB.txt"],
                "output_file": "JKBAAB_single.graphml",
                "expected": "expected_JKBAAB_peptide_single.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pep", "B", "-of", "JKBAAB_single"],
            },
            {
                "name": "Canonical-Merge-Peptide-Overlapping",
                "input_files": ["JKBAAB.txt"],
                "output_file": "JKBAAB_overlapping.graphml",
                "expected": "expected_JKBAAB_peptide_overlapping.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pep", "AB", "-pep", "BC", "-of", "JKBAAB_overlapping"],
            },
            {
                "name": "Canonical-Merge-Peptide-Children",
                "input_files": ["JKBAAB.txt"],
                "output_file": "JKBAAB_children.graphml",
                "expected": "expected_JKBAAB_peptide_children.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pep", "B", "-pep", "ABC", "-of", "JKBAAB_children"],
            },
            {
                "name": "Variant-Spanning-Peptide",
                "input_files": ["JKBAAC.txt"],
                "output_file": "JKBAAC_spanning.graphml",
                "expected": "expected_JKBAAC_peptide_spanning.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pep", "AVC", "-of", "JKBAAC_spanning"],
            },
            {
                "name": "Peptide-CSV",
                "input_files": ["JKBAAA.txt"],
                "output_file": "JKBAAA_csv.graphml",
                "expected": "expected_JKBAAA_peptide.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pf", "tests/test_data/input/peptide.csv", "-of", "JKBAAA_csv"],
            },
            {
                "name": "Peptide-CSV-Multiple",
                "input_files": ["JKBAAB.txt"],
                "output_file": "JKBAAB_csv.graphml",
                "expected": "expected_JKBAAB_peptide_overlapping.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pf", "tests/test_data/input/multiple_peptide.csv", "-of", "JKBAAB_csv"],
            },
            {
                "name": "Peptide-CSV-Isoform",
                "input_files": ["JKBAAA.txt"],
                "output_file": "JKBAAA_csv.graphml",
                "expected": "expected_JKBAAA_peptide.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pf", "tests/test_data/input/peptide_isoform.csv", "-of", "JKBAAA_csv"],
            },
            {
                "name": "Peptide-Intensity-CSV",
                "input_files": ["JKBAAA.txt"],
                "output_file": "JKBAAA_intensity_csv.graphml",
                "expected": "expected_JKBAAA_intensity.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pf", "tests/test_data/input/peptide_intensity.csv", "-int", "-of", "JKBAAA_intensity_csv"],
            },
            {
                "name": "Isoform-Spanning-Peptide",
                "input_files": ["JKBAIA.txt"],
                "output_file": "JKBAIA_spanning.graphml",
                "expected": "expected_JKBAIA_peptide_spanning_isoform.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pep", "AXYZC", "-of", "JKBAIA_spanning"],
            },
            {
                "name": "Variant-into-Isoform",
                "input_files": ["JKXXIH.txt"],
                "expected": "expected_JKXXIH.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
        ]

    def tearDown(self):
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)

    def run_test_case(self, case):
        input_paths = [os.path.join(self.input_dir, fname) for fname in case["input_files"]]
        expected_path = os.path.join(self.expected_dir, case["expected"])
        output_path = os.path.join(self.output_dir, get_file_from_name(case["input_files" if not "output_file" in case.keys() else "output_file"]))

        args = protgraph.parse_args(case["extra_args"] + self.procs_num + input_paths)
        protgraph.prot_graph(**args)

        print(output_path)
        self.assertTrue(os.path.exists(output_path), "Output file not created.")
        self.assertTrue(
            filecmp.cmp(output_path, expected_path, shallow=False),
            f"Output mismatch in test '{case['name']}'"
        )

    def test_all_cases(self):
        for case in self.test_cases:
            with self.subTest(case=case["name"]):
                self.run_test_case(case)

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
                "input_files": ["JKXXAA.txt"],
                "expected": "expected_JKXXAA.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Canonical-Merge",
                "input_files": ["JKXXAB.txt"],
                "expected": "expected_JKXXAB.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Variant",
                "input_files": ["JKXXAC.txt"],
                "expected": "expected_JKXXAC.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Mutation",
                "input_files": ["JKXXAD.txt"],
                "expected": "expected_JKXXAD.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Conflict",
                "input_files": ["JKXXAE.txt"],
                "expected": "expected_JKXXAE.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Variant-Missing",
                "input_files": ["JKXXAF.txt"],
                "expected": "expected_JKXXAF.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Isoform",
                "input_files": ["JKXXIA.txt"],
                "expected": "expected_JKXXIA.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Isoform-Canonical-Merge",
                "input_files": ["JKXXIB.txt"],
                "expected": "expected_JKXXIB.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Shared-VARSEQ",
                "input_files": ["JKXXIC.txt"],
                "expected": "expected_JKXXIC.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Isoform-into-Variant",
                "input_files": ["JKXXID.txt"],
                "expected": "expected_JKXXID.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Partially-Shared-VARSEQs",
                "input_files": ["JKXXIE.txt"],
                "expected": "expected_JKXXIE.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Isoform-and-Variant-Simultaneous",
                "input_files": ["JKXXIF.txt"],
                "expected": "expected_JKXXIF.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Isoform-Missing",
                "input_files": ["JKXXIG.txt"],
                "expected": "expected_JKXXIG.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Canonical-Peptide",
                "input_files": ["JKXXAA.txt"],
                "output_file": "JKXXAA_peptide.graphml",
                "expected": "expected_JKXXAA_peptide.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pep", "A", "-of", "JKXXAA_peptide"],
            },
            {
                "name": "Canonical-Merge-Peptide-Single-Cut",
                "input_files": ["JKXXAB.txt"],
                "output_file": "JKXXAB_single.graphml",
                "expected": "expected_JKXXAB_peptide_single.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pep", "B", "-of", "JKXXAB_single"],
            },
            {
                "name": "Canonical-Merge-Peptide-Overlapping",
                "input_files": ["JKXXAB.txt"],
                "output_file": "JKXXAB_overlapping.graphml",
                "expected": "expected_JKXXAB_peptide_overlapping.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pep", "AB", "-pep", "BC", "-of", "JKXXAB_overlapping"],
            },
            {
                "name": "Canonical-Merge-Peptide-Children-False",
                "input_files": ["JKXXAB.txt"],
                "output_file": "JKXXAB_children.graphml",
                "expected": "expected_JKXXAB_peptide_children.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pep", "B", "-pep", "ABC", "-of", "JKXXAB_children"],
            },
            {
                "name": "Variant-Spanning-Peptide",
                "input_files": ["JKXXAC.txt"],
                "output_file": "JKXXAC_spanning.graphml",
                "expected": "expected_JKXXAC_peptide_spanning.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pep", "AVC", "-of", "JKXXAC_spanning"],
            },
            {
                "name": "Peptide-CSV",
                "input_files": ["JKXXAA.txt"],
                "output_file": "JKXXAA_csv.graphml",
                "expected": "expected_JKXXAA_peptide.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pf", "tests/test_data/input/peptide.csv", "-of", "JKXXAA_csv"],
            },
            {
                "name": "Peptide-CSV-Multiple",
                "input_files": ["JKXXAB.txt"],
                "output_file": "JKXXAB_csv.graphml",
                "expected": "expected_JKXXAB_peptide_overlapping.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pf", "tests/test_data/input/multiple_peptide.csv", "-of", "JKXXAB_csv"],
            },
            {
                "name": "Peptide-CSV-Isoform",
                "input_files": ["JKXXAA.txt"],
                "output_file": "JKXXAA_csv.graphml",
                "expected": "expected_JKXXAA_peptide.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pf", "tests/test_data/input/peptide_isoform.csv", "-of", "JKXXAA_csv"],
            },
            {
                "name": "Peptide-Intensity-CSV",
                "input_files": ["JKXXAA.txt"],
                "output_file": "JKXXAA_intensity_csv.graphml",
                "expected": "expected_JKXXAA_intensity.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pf", "tests/test_data/input/peptide_intensity.csv", "-int", "-of", "JKXXAA_intensity_csv"],
            },
            {
                "name": "Isoform-Spanning-Peptide",
                "input_files": ["JKXXIA.txt"],
                "output_file": "JKXXIA_spanning.graphml",
                "expected": "expected_JKXXIA_peptide_spanning_isoform.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pep", "AXYZC", "-of", "JKXXIA_spanning"],
            },
            {
                "name": "Variant-into-Isoform",
                "input_files": ["JKXXIH.txt"],
                "expected": "expected_JKXXIH.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip"],
            },
            {
                "name": "Canonical-Merge-Peptide-Children-True",
                "input_files": ["JKXXAB.txt"],
                "output_file": "JKXXAB_children_true.graphml",
                "expected": "expected_JKXXAB_peptide_children_true.graphml",
                "extra_args": ["-egraphml", "--export_output_folder", cls.output_dir, "--digestion", "skip", "-sg", "-pep", "B", "-pep", "ABC", "-of", "JKXXAB_children_true", "-mp"],
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

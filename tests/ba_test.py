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
        ]

    def tearDown(self):
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)

    def run_test_case(self, case):
        input_paths = [os.path.join(self.input_dir, fname) for fname in case["input_files"]]
        expected_path = os.path.join(self.expected_dir, case["expected"])
        output_path = os.path.join(self.output_dir, get_file_from_name(case["input_files"]))

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

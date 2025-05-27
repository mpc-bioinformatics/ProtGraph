import filecmp
import os
import unittest

import protgraph


class FunctionalTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        base_dir = os.path.dirname(__file__)
        cls.input_dir = os.path.join(base_dir, "data", "input")
        cls.expected_dir = os.path.join(base_dir, "data", "expected")
        cls.output_file = "output.graphml"
        cls.procs_num = ["-n", "1"]

        cls.test_cases = [
            {
                "name": "Canonical",
                "input_files": ["JKBAAA.txt"],
                "expected": "expected_JKBAAA.graphml",
                "extra_args": ["-egraphml"],
            },
            {
                "name": "Canonical-Merge",
                "input_files": ["JKBAAB.txt"],
                "expected": "expected_JKBAAB.graphml",
                "extra_args": ["-egraphml"],
            },
            {
                "name": "Variant",
                "input_files": ["JKBAAC.txt"],
                "expected": "expected_JKBAAC.graphml",
                "extra_args": ["-egraphml"],
            },
            {
                "name": "Mutation",
                "input_files": ["JKBAAD.txt"],
                "expected": "expected_fm_npepterm.graphml",
                "extra_args": ["-egraphml"],
            },
            {
                "name": "Conflict",
                "input_files": ["JKBAAE.txt"],
                "expected": "expected_fm_npepterm.graphml",
                "extra_args": ["-egraphml"],
            },
            {
                "name": "Isoform",
                "input_files": ["JKBAIA.txt"],
                "expected": "expected_fm_npepterm.graphml",
                "extra_args": ["-egraphml"],
            },
            {
                "name": "Isoform-Canonical-Merge",
                "input_files": ["JKBAIB.txt"],
                "expected": "expected_fm_npepterm.graphml",
                "extra_args": ["-egraphml"],
            },
            {
                "name": "Shared-VARSEQ",
                "input_files": ["JKBAIC.txt"],
                "expected": "expected_fm_npepterm.graphml",
                "extra_args": ["-egraphml"],
            },
            {
                "name": "Isoform-into-Variant",
                "input_files": ["JKBAID.txt"],
                "expected": "expected_fm_npepterm.graphml",
                "extra_args": ["-egraphml"],
            },
            {
                "name": "Partially-Shared-VARSEQs",
                "input_files": ["JKBAIE.txt"],
                "expected": "expected_fm_npepterm.graphml",
                "extra_args": ["-egraphml"],
            },
            {
                "name": "Isoform-and-Variant-Simultaneous",
                "input_files": ["JKBAIF.txt"],
                "expected": "expected_fm_npepterm.graphml",
                "extra_args": ["-egraphml"],
            },
        ]

    def tearDown(self):
        if os.path.exists(self.output_file):
            os.remove(self.output_file)

    def run_test_case(self, case):
        input_paths = [os.path.join(self.input_dir, fname) for fname in case["input_files"]]
        expected_path = os.path.join(self.expected_dir, case["expected_file"])

        args = protgraph.parse_args(case["extra_args"] + self.procs_num + input_paths)
        protgraph.prot_graph(**args)

        self.assertTrue(os.path.exists(self.output_file), "Output file not created.")
        self.assertTrue(
            filecmp.cmp(self.output_file, expected_path, shallow=False),
            f"Output mismatch in test '{case['name']}'"
        )

    def test_all_cases(self):
        for case in self.test_cases:
            with self.subTest(case=case["name"]):
                self.run_test_case(case)

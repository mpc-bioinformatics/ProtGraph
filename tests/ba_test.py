import filecmp
import os
import unittest

import protgraph


class FunctionalTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.test_cases_base = os.path.join(os.path.dirname(__file__), "data")
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
                "extra_args": ["-fm", "C:57.021464"],
            },
            {
                "name": "Variant",
                "input_files": ["e_coli.dat", "p53_human.txt"],
                "expected": "expected_fm_none.graphml",
                "extra_args": ["-ft", "none", "-fm", "c:57.021464"],
            },
            {
                "name": "Mutation",
                "input_files": ["e_coli.dat", "p53_human.txt"],
                "expected": "expected_fm_npepterm.graphml",
                "extra_args": ["-fm", "nPePtErm:57.021464"],
            },
            {
                "name": "Mutation",
                "input_files": ["e_coli.dat", "p53_human.txt"],
                "expected": "expected_fm_npepterm.graphml",
                "extra_args": ["-fm", "nPePtErm:57.021464"],
            },
        ]

    def tearDown(self):
        if os.path.exists(self.output_file):
            os.remove(self.output_file)

    def run_test_case(self, case):
        input_paths = [os.path.join(self.test_cases_base, fname) for fname in case["input_files"]]
        expected_path = os.path.join(self.test_cases_base, case["expected"])

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

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

    def test_fm(self):
        args = protgraph.parse_args(["-fm", "C:57.021464"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_fm_none(self):
        args = protgraph.parse_args(["-ft", "none", "-fm", "c:57.021464"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_fm_npepterm(self):
        args = protgraph.parse_args(["-fm", "nPePtErm:57.021464"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_fm_none_cpepterm(self):
        args = protgraph.parse_args(["-ft", "none", "-fm", "CPePtERM:57.021464"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_fm_none_cprotterm(self):
        args = protgraph.parse_args(["-ft", "none", "-fm", "CPRotERM:57.021464"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_vm_cpepterm(self):
        args = protgraph.parse_args(["-vm", "cPEPtErm:57.021464"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_vm_cproterm(self):
        args = protgraph.parse_args(["-vm", "cprotErm:57.021464"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_vm_none_npepterm(self):
        args = protgraph.parse_args(["-ft", "none", "-vm", "nPePtERM:57.021464"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_vm(self):
        args = protgraph.parse_args(["-vm", "C:57.021464"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_vm_none(self):
        args = protgraph.parse_args(["-ft", "none", "-vm", "c:57.021464"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_vm_nproterm(self):
        args = protgraph.parse_args(["-vm", "nPRotErm:57.021464"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_vm_npepterm(self):
        args = protgraph.parse_args(["-vm", "nPePtErm:57.021464"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_vm_none_cpepterm(self):
        args = protgraph.parse_args(["-ft", "none", "-vm", "CprotERM:57.021464"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_vm_none_cproterm(self):
        args = protgraph.parse_args(["-ft", "none", "-vm", "CpeptERM:57.021464"] + self.procs_num + self.example_files)
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

    def test_statistics_var(self):
        args = protgraph.parse_args(["-cnpvar"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_statistics_mut(self):
        args = protgraph.parse_args(["-cnpmut"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_statistics_con(self):
        args = protgraph.parse_args(["-cnpcon"] + self.procs_num + self.example_files)
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
        args = protgraph.parse_args(["-egml"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_export_pickle(self):
        args = protgraph.parse_args(["-epickle"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_export_pcsr(self):
        args = protgraph.parse_args(["-epcsr"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_export_bpcsr(self):
        args = protgraph.parse_args(["-ebpcsr"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_export_lpcsr(self):
        args = protgraph.parse_args(["-elpcsr"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_export_lbpcsr(self):
        args = protgraph.parse_args(["-elbpcsr"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_export_pep_fasta(self):
        args = protgraph.parse_args(["-epepfasta", "--pep_hops", "2"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_export_pep_trie(self):
        args = protgraph.parse_args(
            ["-epeptrie", "--pep_hops", "5"] + self.procs_num + [os.path.join(self.examples_path, "Q9QXS1.txt")]
        )
        protgraph.prot_graph(**args)

    def test_export_pep_sqlite(self):
        args = protgraph.parse_args(["-epepsqlite", "--pep_hops", "2"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_issue8(self):
        args = protgraph.parse_args(["-n", "1", os.path.join(self.examples_path, "Q9QXS1.txt")])
        protgraph.prot_graph(**args)

    def test_issue13(self):
        args = protgraph.parse_args(["-n", "1", os.path.join(self.examples_path, "F1SN05.txt")])
        protgraph.prot_graph(**args)

    def test_issue41(self):
        args = protgraph.parse_args(["-n", "1", "-epepfasta", os.path.join(self.examples_path, "P49782.txt")])
        protgraph.prot_graph(**args)

    def _run_proforma(self, extra_args, out):
        """ ProForma export of tests/fixtures/minimal.txt (MKCTMSAK, VARIANT T4M). """
        fixture = os.path.abspath(os.path.join(os.path.dirname(__file__), "fixtures", "minimal.txt"))
        args = protgraph.parse_args(
            ["-epepproforma", "--pep_proforma_out", out, "--pep_hops", "10", "-d", "trypsin"]
            + extra_args + self.procs_num + [fixture]
        )
        protgraph.prot_graph(**args)
        with open(out) as f:
            return [ln.rstrip("\n").split("\t") for ln in f if ln.strip()]

    def test_export_pep_proforma_smoke(self):
        args = protgraph.parse_args(["-epepproforma", "--pep_hops", "2"] + self.procs_num + self.example_files)
        protgraph.prot_graph(**args)

    def test_export_pep_proforma_pins_mod_placement(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            lines = self._run_proforma(
                ["-ft", "VARIANT", "-fm", "C:57.021464", "-vm", "M:15.994915",
                 "--pep_proforma_mod_names", "15.994915=UNIMOD:35"],
                os.path.join(tmp, "p.tsv"),
            )
            assert all(len(x) == 5 for x in lines)
            rows = {tuple(x) for x in lines}
            proformas = {x[1] for x in lines}
            # exact placement with coordinates: the fixed mod sits on the C
            assert ("X9TEST", "C[+57.021464]TMSAK", "3", "8", "0",) in rows
            # a variable mod on a VARIANT-introduced residue (no reference
            # position!) is placed on that residue, not dropped
            assert "C[+57.021464]M[UNIMOD:35]MSAK" in proformas
            # variable mods also export their unmodified counterparts
            assert {"MK", "M[UNIMOD:35]K"} <= proformas
            # the fixed mod is applied to every C, always signed
            assert all("C[+57.021464]" in p for p in proformas if "C" in p)

    def test_export_pep_proforma_terminal_mods_and_overwrite(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "p.tsv")
            lines = self._run_proforma(["-fm", "NPEPTERM:42.010565", "-vm", "CPEPTERM:79.966"], out)
            proformas = {x[1] for x in lines}
            proformas.discard("proforma")
            assert all(p.startswith("[+42.010565]-") for p in proformas)
            # variable C-terminal: modified and unmodified paths both exported
            assert any(p.endswith("-[+79.966]") for p in proformas)
            assert any(not p.endswith("-[+79.966]") for p in proformas)
            # a rerun overwrites instead of appending stale results
            again = self._run_proforma(["-fm", "NPEPTERM:42.010565", "-vm", "CPEPTERM:79.966"], out)
            assert len(again) == len(lines)

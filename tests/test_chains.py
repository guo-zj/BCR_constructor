"""Regression tests for chain typing and V-domain trimming.

    python3 -m unittest discover -s tests -v
"""

import pathlib
import unittest

from bcr import chains, read_fasta

FIXTURE = pathlib.Path(__file__).parent / "fixtures" / "example_input.fasta"

LEADER_H = "MGWSCIILFLVATATGVHS"
CH1 = "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSW"


class FixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.seqs = {r.name: r.seq for r in read_fasta(FIXTURE)}

    def call(self, name):
        return chains.process(name, self.seqs[name])

    def test_clean_vh_passes_unchanged(self):
        r = self.call("mabA_VH")
        self.assertTrue(r.ok)
        self.assertEqual(r.detected, "heavy")
        self.assertEqual((r.n_trimmed, r.c_trimmed), (0, 0))
        self.assertEqual(r.clean_seq, self.seqs["mabA_VH"])
        self.assertEqual(r.confidence, "high")

    def test_clean_vl_typed_as_kappa(self):
        r = self.call("mabA_VL")
        self.assertTrue(r.ok)
        self.assertEqual(r.detected, "light")
        self.assertEqual(r.subclass, "kappa")
        self.assertEqual((r.n_trimmed, r.c_trimmed), (0, 0))

    def test_leader_and_ch1_are_stripped(self):
        r = self.call("mabB_VH")
        self.assertTrue(r.ok)
        self.assertEqual(r.n_trimmed, len(LEADER_H))
        self.assertEqual(r.c_trimmed, len(CH1))
        # The trimmed domain must match the independently clean record exactly.
        self.assertEqual(r.clean_seq, self.seqs["mabA_VH"])

    def test_lambda_leader_and_cl_are_stripped(self):
        r = self.call("mabC_VL")
        self.assertTrue(r.ok)
        self.assertEqual(r.detected, "light")
        self.assertEqual(r.subclass, "lambda")
        self.assertTrue(r.clean_seq.startswith("QSVLTQ"))
        self.assertTrue(r.clean_seq.endswith("FGGGTQLTVL"))

    def test_mislabelled_chain_is_rejected(self):
        r = self.call("mabD_VH")
        self.assertFalse(r.ok)
        self.assertEqual(r.detected, "light")
        self.assertEqual(r.declared, "heavy")
        self.assertEqual(r.clean_seq, "")
        self.assertIn("declares heavy", r.errors[0])

    def test_non_antibody_is_rejected(self):
        r = self.call("mabE_VH")
        self.assertFalse(r.ok)
        self.assertIn("framework 4", r.errors[0])


class BoundaryRuleTests(unittest.TestCase):
    """The counting rules themselves: Trp + 10 and Phe + 9."""

    def test_heavy_ends_at_trp_plus_ten(self):
        seq = read_fasta(FIXTURE)
        vh = {r.name: r.seq for r in seq}["mabA_VH"]
        r = chains.process("x_VH", vh)
        trp = vh.index("WGQGTLVTVSS")
        self.assertEqual(r.v_end, trp + 10)
        self.assertEqual(vh[r.v_end], "S")

    def test_light_ends_at_phe_plus_nine(self):
        vl = {r.name: r.seq for r in read_fasta(FIXTURE)}["mabA_VL"]
        r = chains.process("x_VL", vl)
        phe = vl.index("FGQGTKVEIK")
        self.assertEqual(r.v_end, phe + 9)
        self.assertEqual(vl[r.v_end], "K")

    def test_murine_style_vtvsa_is_not_dropped(self):
        """Section 2's warning: a literal TVSS match loses ...VTVSA."""
        vh = {r.name: r.seq for r in read_fasta(FIXTURE)}["mabA_VH"]
        mutated = vh[:-1] + "A"
        r = chains.process("x_VH", mutated)
        self.assertTrue(r.ok)
        self.assertEqual(r.v_end, len(mutated) - 1)

    def test_all_human_ighj_fr4_are_matched(self):
        vh = {r.name: r.seq for r in read_fasta(FIXTURE)}["mabA_VH"]
        stem = vh[: vh.index("WGQGTLVTVSS")]
        for fr4 in (
            "WGQGTLVTVSS",
            "WGRGTLVTVSS",
            "WGQGTMVTVSS",
            "WGQGTTVTVSS",
        ):
            with self.subTest(fr4=fr4):
                r = chains.process("x_VH", stem + fr4)
                self.assertTrue(r.ok, fr4)
                self.assertEqual(r.detected, "heavy")

    def test_all_human_light_j_fr4_are_matched(self):
        vl = {r.name: r.seq for r in read_fasta(FIXTURE)}["mabA_VL"]
        stem = vl[: vl.index("FGQGTKVEIK")]
        expected = {
            "FGQGTKVEIK": "kappa",
            "FGQGTKLEIK": "kappa",
            "FGPGTKVDIK": "kappa",
            "FGGGTKVEIK": "kappa",
            "FGQGTRLEIK": "kappa",
            "FGTGTKVTVL": "lambda",
            "FGGGTKLTVL": "lambda",
            "FGSGTKVTVL": "lambda",
            "FGGGTQLTVL": "lambda",
        }
        for fr4, subclass in expected.items():
            with self.subTest(fr4=fr4):
                r = chains.process("x_VL", stem + fr4)
                self.assertTrue(r.ok, fr4)
                self.assertEqual(r.detected, "light")
                self.assertEqual(r.subclass, subclass)


class HybridTests(unittest.TestCase):
    def test_lambda_v_on_kappa_c_is_flagged(self):
        """Section 3's trap: nine genuine lambda CL residues, then kappa."""
        vl = {r.name: r.seq for r in read_fasta(FIXTURE)}["mabC_VL"]
        v_only = chains.process("x_VL", vl).clean_seq
        hybrid = v_only + "GQPKAAPSV" + "FIFPPSDEQLKSGTASVVC"
        r = chains.process("x_VL", hybrid)
        self.assertTrue(r.ok)
        self.assertEqual(r.clean_seq, v_only)
        self.assertTrue(
            any("hybrid" in w for w in r.warnings),
            f"expected a hybrid warning, got {r.warnings}",
        )


class HeaderTests(unittest.TestCase):
    def test_suffix_parsing(self):
        self.assertEqual(chains.parse_declared("mabA_VH"), ("mabA", "heavy"))
        self.assertEqual(chains.parse_declared("mab_A_VL"), ("mab_A", "light"))
        self.assertEqual(chains.parse_declared("mabA_vh"), ("mabA", "heavy"))
        self.assertEqual(chains.parse_declared("mabA"), ("mabA", None))

    def test_unlabelled_header_is_typed_but_warned(self):
        vh = {r.name: r.seq for r in read_fasta(FIXTURE)}["mabA_VH"]
        r = chains.process("mystery", vh)
        self.assertTrue(r.ok)
        self.assertEqual(r.detected, "heavy")
        self.assertIsNone(r.declared)
        self.assertTrue(any("_VH/_VL" in w for w in r.warnings))

    def test_bad_characters_are_rejected(self):
        r = chains.process("x_VH", "EVQLVESGGG123")
        self.assertFalse(r.ok)
        self.assertIn("non-amino-acid", r.errors[0])


if __name__ == "__main__":
    unittest.main()

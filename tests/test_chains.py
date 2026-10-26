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


class LongCdr3Tests(unittest.TestCase):
    """A domain longer than the expected range must not have FR1 amputated.

    `len_range` bounds which start positions are scored at all, so a cap set
    too tightly does not merely lower confidence -- the true N-terminus is
    never considered and the best start inside the window cuts into FR1. C05
    (anti-influenza, 26-residue CDR-H3) is a 138 aa VH that is already clean,
    and was being trimmed by 15 residues on a medium-confidence call.
    """

    C05_VH = (
        "EVQLQESGGGLVQPGESLRLSCVGSGSSFGESTLSYYAVSWVRQAPGKGLEWLSIINAGGGDIDYADS"
        "VEGRFTISRDNSKETLYLQMTNLRVEDTGVYYCAKHMSMQQVVSAGWERADLVGDAFDVWGQGTMVTVSA"
    )

    def test_a_long_cdr3_heavy_chain_is_not_trimmed_into_fr1(self):
        r = chains.process("C05_VH", self.C05_VH)
        self.assertTrue(r.ok)
        self.assertEqual(r.detected, "heavy")
        self.assertEqual(r.n_trimmed, 0)
        self.assertEqual(r.c_trimmed, 0)
        self.assertEqual(r.clean_seq, self.C05_VH)

    def test_the_true_n_terminus_wins_on_its_own_anchors(self):
        r = chains.process("C05_VH", self.C05_VH)
        self.assertEqual(r.confidence, "high")
        # The anchors that the amputated start lost: SG at 7 (section 1b) and
        # the conserved FR1 Cys (section 4 cross-check).
        self.assertTrue(any("SG at residue 7" in e for e in r.evidence))
        self.assertTrue(any("FR1 Cys at residue" in e for e in r.evidence))

    def test_an_unusual_length_is_still_reported(self):
        """Widening the scoring window does not hide the length itself."""
        r = chains.process("C05_VH", self.C05_VH)
        self.assertEqual(len(r.clean_seq), 138)
        self.assertTrue(
            any("outside the usual 98-128" in w for w in r.warnings), r.warnings
        )

    def test_the_v_end_rule_is_unchanged_by_the_wider_window(self):
        """Still Trp + 10: the C-terminal boundary is exact either way."""
        r = chains.process("C05_VH", self.C05_VH)
        self.assertTrue(r.clean_seq.endswith("WGQGTMVTVSA"))
        self.assertEqual(r.v_end, len(self.C05_VH) - 1)


class AcceptAsGivenTests(unittest.TestCase):
    """Overruling a failure stores the sequence as it came in, and nothing else.

    chains.accept_as_given is the only way a rejected record reaches the
    store, so what it does and refuses to do is the whole guarantee.
    """

    GFP = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTF"

    def test_the_sequence_is_stored_whole_and_untrimmed(self):
        r = chains.accept_as_given(chains.process("mab_VH", self.GFP))
        self.assertTrue(r.ok)
        self.assertTrue(r.accepted_as_given)
        self.assertEqual(r.clean_seq, self.GFP)
        self.assertEqual((r.n_trimmed, r.c_trimmed), (0, 0))
        self.assertEqual((r.v_start, r.v_end), (0, len(self.GFP) - 1))

    def test_nothing_is_derived_from_it(self):
        r = chains.accept_as_given(chains.process("mab_VH", self.GFP))
        self.assertEqual(r.confidence, "none")
        self.assertIsNone(r.subclass)

    def test_a_light_chain_has_to_say_kappa_or_lambda(self):
        """Both of section 3's checkpoints need a V end, and an overruled
        record has no boundary at all -- so the user states it."""
        failed = chains.process("mab_VL", self.GFP)
        with self.assertRaises(ValueError) as caught:
            chains.accept_as_given(failed)
        self.assertIn("kappa or lambda", str(caught.exception))
        with self.assertRaises(ValueError):
            chains.accept_as_given(failed, "igkc")

    def test_a_stated_subclass_is_recorded_as_stated(self):
        r = chains.accept_as_given(chains.process("mab_VL", self.GFP), "lambda")
        self.assertEqual(r.subclass, "lambda")
        self.assertTrue(any("stated by you, not derived" in w
                            for w in r.warnings), r.warnings)

    def test_kappa_lambda_does_not_apply_to_a_heavy_chain(self):
        with self.assertRaises(ValueError) as caught:
            chains.accept_as_given(chains.process("mab_VH", self.GFP), "kappa")
        self.assertIn("does not apply", str(caught.exception))

    def test_the_errors_survive_as_warnings(self):
        """The verdict was overruled, not withdrawn: it was still true."""
        r = chains.accept_as_given(chains.process("mab_VH", self.GFP))
        self.assertEqual(r.errors, [])
        self.assertTrue(any("did not pass the sanity check" in w
                            for w in r.warnings))
        self.assertTrue(any("framework 4" in w for w in r.warnings))

    def test_the_header_label_wins_over_the_detected_class(self):
        """A mislabelled record is stored as the chain its header claims: the
        user is insisting on their own label."""
        light = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGV"
                 "PSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
        failed = chains.process("mab_VH", light)
        self.assertEqual(failed.detected, "light")
        self.assertFalse(failed.ok)
        self.assertEqual(chains.accept_as_given(failed).detected, "heavy")

    def test_a_passing_record_is_not_a_candidate(self):
        passed = chains.process(
            "mab_VL",
            "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGV"
            "PSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK",
        )
        self.assertTrue(passed.ok)
        chain, reason = chains.can_accept_as_given(passed)
        self.assertIsNone(chain)
        self.assertIn("passed the sanity check", reason)

    def test_a_non_sequence_is_refused(self):
        for seq, expect in (("QVQL123@@", "not amino acids"), ("", "no sequence")):
            with self.subTest(seq=seq):
                failed = chains.process("mab_VH", seq)
                chain, reason = chains.can_accept_as_given(failed)
                self.assertIsNone(chain)
                self.assertIn(expect, reason)
                with self.assertRaises(ValueError):
                    chains.accept_as_given(failed)


if __name__ == "__main__":
    unittest.main()

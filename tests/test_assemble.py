"""Tests for BCR construct assembly.

    python3 -m unittest discover -s tests -t . -v
"""

import pathlib
import unittest

from bcr import assemble, chains, read_fasta

ROOT = pathlib.Path(__file__).resolve().parents[1]
PARTS = pathlib.Path(__file__).parent / "fixtures" / "parts_sc.fasta"
FIXTURE = pathlib.Path(__file__).parent / "fixtures" / "example_input.fasta"


def typed(path):
    return [chains.process(r.name, r.seq) for r in read_fasta(path)]


class PartsTests(unittest.TestCase):
    def test_real_misc_file_validates(self):
        parts, problems = assemble.load_parts(PARTS)
        self.assertEqual(problems, [])
        for name in assemble.REQUIRED_PARTS:
            self.assertIn(name, parts)

    def test_constant_parts_start_at_their_own_first_residue(self):
        """Section 6: start the constant region at its own first residue."""
        parts, _ = assemble.load_parts(PARTS)
        self.assertTrue(parts["IGKC"].startswith("RTVAAPSV"))
        self.assertTrue(parts["IGLC2"].startswith("GQPKAAPSV"))
        self.assertTrue(parts["CHandCTSC"].startswith("ASTKGPSVFPLAP"))

    def test_trailing_stop_is_stripped_from_parts(self):
        """Parts hold residues only; the layout contributes the stop codon."""
        parts, problems = assemble.load_parts(PARTS)
        self.assertEqual(problems, [])
        self.assertNotIn("*", parts["CHandCTSC"])
        self.assertEqual(len(parts["CHandCTSC"]), 399)

    def test_only_the_light_constants_are_required(self):
        """Part names come from the FASTA, so only IGKC/IGLC2 are demanded.

        A layout chooses between those two per antibody, so the code has to
        know their names; everything else is whatever the file calls it.
        """
        import tempfile

        with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as fh:
            fh.write(">SignalHC\nMGWSCIILFLVATATGVHS\n")
            path = fh.name
        _, problems = assemble.load_parts(path)
        self.assertTrue(any("IGKC" in p for p in problems), problems)
        self.assertTrue(any("IGLC2" in p for p in problems), problems)
        # No complaint about a linker or a heavy constant: nothing says a
        # parts file has to contain one, or what it would be called.
        self.assertFalse(any("Linker" in p for p in problems), problems)
        self.assertFalse(any("CHandCT" in p for p in problems), problems)

    def test_arbitrary_part_names_are_accepted(self):
        import tempfile

        kappa = "RTVAAPSVFIFPPSDEQLK"
        lam = "GQPKAAPSVTLFPPSSEELQ"
        with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as fh:
            fh.write(f">IGKC\n{kappa}\n>IGLC2\n{lam}\n>whatever_I_like\nGGGGSGGGGS\n")
            path = fh.name
        parts, problems = assemble.load_parts(path)
        self.assertEqual(problems, [])
        self.assertIn("whatever_I_like", parts)

    def test_truncated_part_with_internal_stop_is_reported(self):
        import tempfile

        with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as fh:
            fh.write(">IGKC\nRTVAAPSV*FIFPPSDEQLK\n")
            path = fh.name
        _, problems = assemble.load_parts(path)
        self.assertTrue(
            any("truncated" in p for p in problems), f"got {problems}"
        )


class BuildTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.parts, problems = assemble.load_parts(PARTS)
        assert not problems, problems

    def build_from(self, path, name):
        results = [r for r in typed(path) if r.ok]
        chain_map = assemble.pair_chains(results)[name]
        return assemble.build(
            name, chain_map.get("heavy"), chain_map.get("light"), self.parts
        )

    def test_order_is_signal_vl_cl_linker_vh_chandct(self):
        c = self.build_from(ROOT / "sequences.fasta", "GL_VRC01")
        self.assertTrue(c.ok, c.errors)
        labels = [s.label.split(":")[0] for s in c.segments]
        self.assertEqual(
            labels,
            ["SignalHC", "VL", "IGKC", "LinkerSC", "VH", "CHandCTSC", "Stop"],
        )

    def test_construct_is_exact_concatenation_of_its_parts(self):
        c = self.build_from(ROOT / "sequences.fasta", "GL_VRC01")
        self.assertEqual(c.seq, "".join(s.seq for s in c.segments))
        self.assertEqual(len(c.seq), sum(len(s.seq) for s in c.segments))

    def test_segment_coordinates_address_the_right_residues(self):
        c = self.build_from(ROOT / "sequences.fasta", "GL_VRC01")
        for seg in c.segments:
            self.assertEqual(c.seq[seg.start - 1 : seg.end], seg.seq, seg.label)

    def test_stop_codon_survives_into_the_construct(self):
        c = self.build_from(ROOT / "sequences.fasta", "GL_VRC01")
        self.assertTrue(c.seq.endswith("*"))
        self.assertTrue(c.has_stop)
        # Exactly one, and only at the very end.
        self.assertEqual(c.seq.count("*"), 1)
        self.assertNotIn("*", c.seq[:-1])

    def test_residue_count_excludes_the_stop(self):
        c = self.build_from(ROOT / "sequences.fasta", "GL_VRC01")
        self.assertEqual(c.residue_count, len(c.seq) - 1)
        self.assertIn("+ stop", c.size_note)

    def test_stop_is_appended_although_no_part_carries_one(self):
        """Parts are stop-free, so the layout must supply exactly one."""
        self.assertFalse(any("*" in s for s in self.parts.values()))
        results = [r for r in typed(ROOT / "sequences.fasta") if r.ok]
        chain_map = assemble.pair_chains(results)["GL_VRC01"]
        c = assemble.build(
            "GL_VRC01", chain_map.get("heavy"), chain_map.get("light"), self.parts
        )
        self.assertTrue(c.ok, c.errors)
        self.assertTrue(c.has_stop)
        self.assertEqual(c.seq.count("*"), 1)
        self.assertEqual(c.residue_count, len(c.seq) - 1)

    def test_kappa_vl_takes_igkc(self):
        c = self.build_from(ROOT / "sequences.fasta", "GL_VRC01")
        self.assertEqual(c.subclass, "kappa")
        self.assertIn("IGKC", [s.label for s in c.segments])
        self.assertIn(self.parts["IGKC"], c.seq)
        self.assertNotIn(self.parts["IGLC2"], c.seq)

    def test_lambda_vl_takes_iglc2(self):
        c = self.build_from(FIXTURE, "mabC")
        # mabC has no VH partner in the fixture, so assembly must fail -- but
        # the lambda routing is checked below on a complete pair.
        self.assertFalse(c.ok)

    def test_lambda_pair_routes_to_iglc2(self):
        results = [r for r in typed(FIXTURE) if r.ok]
        by_name = {r.name: r for r in results}
        vh = by_name["mabA_VH"]
        vl = by_name["mabC_VL"]
        c = assemble.build("chimera", vh, vl, self.parts)
        self.assertTrue(c.ok, c.errors)
        self.assertEqual(c.subclass, "lambda")
        self.assertIn("IGLC2", [s.label for s in c.segments])
        self.assertIn(self.parts["IGLC2"], c.seq)
        self.assertNotIn(self.parts["IGKC"], c.seq)

    def test_missing_heavy_chain_fails_cleanly(self):
        results = [r for r in typed(FIXTURE) if r.ok]
        vl = {r.name: r for r in results}["mabA_VL"]
        c = assemble.build("mabA", None, vl, self.parts)
        self.assertFalse(c.ok)
        self.assertEqual(c.seq, "")
        self.assertTrue(any("heavy" in e for e in c.errors))

    def test_missing_light_chain_fails_cleanly(self):
        results = [r for r in typed(FIXTURE) if r.ok]
        vh = {r.name: r for r in results}["mabA_VH"]
        c = assemble.build("mabA", vh, None, self.parts)
        self.assertFalse(c.ok)
        self.assertEqual(c.seq, "")
        self.assertTrue(any("light" in e for e in c.errors))


class JunctionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.parts, _ = assemble.load_parts(PARTS)
        results = [r for r in typed(ROOT / "sequences.fasta") if r.ok]
        chain_map = assemble.pair_chains(results)["GL_VRC01"]
        cls.c = assemble.build(
            "GL_VRC01", chain_map.get("heavy"), chain_map.get("light"), cls.parts
        )

    def test_signal_is_flush_against_the_mature_n_terminus(self):
        """Section 6: no GS linker between leader and mature N-terminus."""
        signal = self.parts["SignalHC"]
        self.assertTrue(self.c.seq.startswith(signal))
        self.assertEqual(self.c.seq[len(signal) : len(signal) + 6], "EIVLTQ")
        self.assertNotEqual(self.c.seq[len(signal) : len(signal) + 2], "GS")

    def test_vl_runs_straight_into_its_constant_region(self):
        self.assertIn("TKLEIK" + "RTVAAPSV", self.c.seq)

    def test_vh_runs_straight_into_ch1(self):
        self.assertIn("VTVSS" + "ASTKGPSVFPLAP", self.c.seq)

    def test_internal_stop_is_an_error(self):
        parts = dict(self.parts)
        parts["LinkerSC"] = "GGSSGS*GGSAGG"
        results = [r for r in typed(ROOT / "sequences.fasta") if r.ok]
        chain_map = assemble.pair_chains(results)["GL_VRC01"]
        c = assemble.build(
            "GL_VRC01", chain_map.get("heavy"), chain_map.get("light"), parts
        )
        self.assertFalse(c.ok)
        self.assertTrue(
            any("stop codon at residue" in e for e in c.errors), c.errors
        )
        self.assertTrue(any("LinkerSC" in e for e in c.errors), c.errors)


class RoleInferenceTests(unittest.TestCase):
    """A part's role is read from its sequence, since names are free-form."""

    def test_constant_regions_are_identified_by_their_opening(self):
        self.assertEqual(
            assemble.infer_role("RTVAAPSVFIFPPSDEQLK"), ("light_constant", "kappa")
        )
        self.assertEqual(
            assemble.infer_role("GQPKAAPSVTLFPPSSEELQ"), ("light_constant", "lambda")
        )
        self.assertEqual(
            assemble.infer_role("ASTKGPSVFPLAPSSKSTSGG"), ("heavy_constant", "")
        )

    def test_role_does_not_depend_on_the_name(self):
        parts, _ = assemble.load_parts(PARTS)
        self.assertEqual(assemble.infer_role(parts["CHandCTSC"])[0], "heavy_constant")
        self.assertEqual(assemble.infer_role(parts["SignalHC"])[0], "signal_peptide")

    def test_a_short_met_sequence_reads_as_a_signal_peptide(self):
        self.assertEqual(
            assemble.infer_role("MGWSCIILFLVATATGVHS"), ("signal_peptide", "")
        )

    def test_anything_unrecognised_is_other_rather_than_guessed(self):
        self.assertEqual(assemble.infer_role("GGGGSGGGGSGGGGS"), ("other", ""))
        # Long and Met-initiated is not enough to call a signal peptide.
        self.assertEqual(assemble.infer_role("M" + "A" * 80), ("other", ""))


if __name__ == "__main__":
    unittest.main()

"""Tests for layout presets, the build plan, and building from the CSV store.

    python3 -m unittest discover -s tests -t . -v
"""

import pathlib
import shutil
import tempfile
import unittest

from bcr import assemble, read_fasta
from bcr.store import Store
from bcr.web import create_app

ROOT = pathlib.Path(__file__).resolve().parents[1]
VRC01 = pathlib.Path(__file__).parent / "fixtures" / "vrc01.fasta"

def clean_vrc01():
    """The VRC01 pair, trimmed. Derived with the real typing code rather than
    read from sequences_clean.fasta, which is a generated working file."""
    from bcr import chains
    from bcr.fasta import read_fasta as _rf

    return {
        r.name: chains.process(r.name, r.seq).clean_seq for r in _rf(VRC01)
    }

PARTS = pathlib.Path(__file__).parent / "fixtures" / "parts_sc.fasta"
FIXTURE = pathlib.Path(__file__).parent / "fixtures" / "example_input.fasta"


class PresetTests(unittest.TestCase):
    def test_single_chain_bcr_is_the_default_preset(self):
        self.assertEqual(assemble.DEFAULT_PRESET, "single_chain_bcr")
        self.assertEqual(
            assemble.PRESETS["single_chain_bcr"]["label"], "Single chain BCR"
        )

    def test_preset_order_matches_the_agreed_layout(self):
        slots = assemble.PRESETS["single_chain_bcr"]["slots"]
        self.assertEqual(
            [assemble.slot_label(s) for s in slots],
            ["SignalHC", "VL", "IGKC / IGLC2", "LinkerSC", "VH", "CHandCTSC", "Stop"],
        )

    def test_light_constant_slot_resolves_per_subclass(self):
        slot = {"kind": "light_constant"}
        self.assertEqual(assemble.slot_label(slot, "kappa"), "IGKC")
        self.assertEqual(assemble.slot_label(slot, "lambda"), "IGLC2")
        self.assertEqual(assemble.slot_label(slot, None), "IGKC / IGLC2")

    def test_describe_layout_reports_part_lengths(self):
        parts, _ = assemble.load_parts(PARTS)
        layout = assemble.describe_layout("single_chain_bcr", parts)
        by_label = {s["label"]: s for s in layout}
        self.assertEqual(by_label["SignalHC"]["length"], 19)
        # Length excludes the stop codon, as elsewhere.
        self.assertEqual(by_label["CHandCTSC"]["length"], 399)
        self.assertEqual(by_label["VL"]["detail"], "per antibody")

    def test_describe_layout_flags_a_part_that_is_not_loaded(self):
        layout = assemble.describe_layout("single_chain_bcr", {})
        signal = next(s for s in layout if s["label"] == "SignalHC")
        self.assertEqual(signal["detail"], "not in the parts table")


class BuildFromLayoutTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.parts, problems = assemble.load_parts(PARTS)
        assert not problems, problems
        clean = clean_vrc01()
        cls.vh = clean["GL_VRC01_VH"]
        cls.vl = clean["GL_VRC01_VL"]

    def test_matches_the_fasta_driven_builder(self):
        """Building from stored sequences must equal building from FASTA."""
        from bcr import chains

        raw = {r.name: r.seq for r in read_fasta(VRC01)}
        vh = chains.process("GL_VRC01_VH", raw["GL_VRC01_VH"])
        vl = chains.process("GL_VRC01_VL", raw["GL_VRC01_VL"])
        from_fasta = assemble.build("GL_VRC01", vh, vl, self.parts)
        from_store = assemble.build_from_layout(
            "GL_VRC01", self.vh, self.vl, "kappa", self.parts
        )
        self.assertTrue(from_store.ok, from_store.errors)
        self.assertEqual(from_store.seq, from_fasta.seq)

    def test_kappa_resolves_to_igkc(self):
        c = assemble.build_from_layout(
            "x", self.vh, self.vl, "kappa", self.parts
        )
        self.assertTrue(c.ok, c.errors)
        cl = next(s for s in c.segments if s.kind == "light_constant")
        self.assertEqual(cl.label, "IGKC")
        self.assertEqual(cl.seq, self.parts["IGKC"])

    def test_lambda_resolves_to_iglc2(self):
        # A lambda V domain ends L, which is what drives the choice.
        lambda_vl = self.vl[:-1] + "L"
        c = assemble.build_from_layout(
            "x", self.vh, lambda_vl, "lambda", self.parts
        )
        self.assertTrue(c.ok, c.errors)
        cl = next(s for s in c.segments if s.kind == "light_constant")
        self.assertEqual(cl.label, "IGLC2")
        self.assertEqual(cl.seq, self.parts["IGLC2"])

    def test_unknown_subclass_is_refused(self):
        c = assemble.build_from_layout("x", self.vh, self.vl, "", self.parts)
        self.assertFalse(c.ok)
        self.assertTrue(any("neither kappa nor lambda" in e for e in c.errors))

    def test_missing_part_is_refused_rather_than_skipped(self):
        parts = dict(self.parts)
        del parts["LinkerSC"]
        c = assemble.build_from_layout("x", self.vh, self.vl, "kappa", parts)
        self.assertFalse(c.ok)
        self.assertEqual(c.seq, "")
        self.assertTrue(any("LinkerSC" in e for e in c.errors))

    def test_segments_carry_kind_and_source(self):
        c = assemble.build_from_layout(
            "GL_VRC01", self.vh, self.vl, "kappa", self.parts
        )
        kinds = [s.kind for s in c.segments]
        self.assertEqual(
            kinds,
            ["part", "vl", "light_constant", "part", "vh", "part", "stop"],
        )
        vl_seg = next(s for s in c.segments if s.kind == "vl")
        self.assertEqual(vl_seg.source, "GL_VRC01")

    def test_unknown_preset_is_refused(self):
        c = assemble.build_from_layout(
            "x", self.vh, self.vl, "kappa", self.parts, preset="nope"
        )
        self.assertFalse(c.ok)
        self.assertTrue(any("unknown layout preset" in e for e in c.errors))


class PlanAndBuildEndpointTests(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.out = pathlib.Path(self.dir) / "out.fasta"
        app = create_app(self.dir)
        app.testing = True
        self.client = app.test_client()
        self.store = Store(self.dir)
        # Populate both tables the way the UI would: parts are previewed then
        # applied, antibodies are checked then approved.
        self.client.post("/api/parts/apply", json={"path": str(PARTS)})
        self._approve(VRC01, ["GL_VRC01"])

    def _approve(self, path, ids):
        self.client.post("/api/check", json={"path": str(path)})
        return self.client.post(
            "/api/approve", json={"path": str(path), "ids": ids}
        )

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_presets_endpoint_exposes_the_default(self):
        d = self.client.get("/api/presets").get_json()
        self.assertEqual(d["default"], "single_chain_bcr")
        self.assertTrue(d["parts_loaded"])
        labels = [s["label"] for s in d["presets"]["single_chain_bcr"]["layout"]]
        self.assertEqual(
            labels,
            ["SignalHC", "VL", "IGKC / IGLC2", "LinkerSC", "VH", "CHandCTSC", "Stop"],
        )

    def test_plan_shows_the_resolved_light_constant(self):
        res = self.client.post("/api/plan", json={"ids": ["GL_VRC01"]})
        self.assertEqual(res.status_code, 200)
        d = res.get_json()
        self.assertEqual(d["buildable"], 1)
        c = d["constructs"][0]
        self.assertEqual(c["cl_part"], "IGKC")
        self.assertEqual(c["subclass"], "kappa")
        self.assertIn("+ stop", c["size_note"])
        self.assertEqual(
            [s["label"] for s in c["segments"]],
            ["SignalHC", "VL:GL_VRC01", "IGKC", "LinkerSC",
             "VH:GL_VRC01", "CHandCTSC", "Stop"],
        )

    def test_plan_writes_nothing(self):
        self.client.post("/api/plan", json={"ids": ["GL_VRC01"]})
        self.assertFalse(self.out.exists())

    def test_plan_rejects_empty_selection(self):
        res = self.client.post("/api/plan", json={"ids": []})
        self.assertEqual(res.status_code, 400)

    def test_plan_reports_an_unknown_antibody(self):
        d = self.client.post("/api/plan", json={"ids": ["nope"]}).get_json()
        self.assertEqual(d["buildable"], 0)
        self.assertIn("not in the antibody library", d["constructs"][0]["errors"][0])

    def test_build_writes_the_fasta(self):
        res = self.client.post(
            "/api/build", json={"ids": ["GL_VRC01"], "output": str(self.out)}
        )
        self.assertEqual(res.status_code, 200)
        d = res.get_json()
        self.assertEqual(d["written"], 1)
        self.assertEqual(d["replaced"], 0)

        written = read_fasta(self.out)
        self.assertEqual(len(written), 1)
        self.assertEqual(written[0].name, "GL_VRC01_BCR")
        self.assertTrue(written[0].seq.endswith("*"))

        # Identical to assembling straight from the store.
        parts = self.store.parts_as_dict()
        row = self.store.antibodies.read_map()["GL_VRC01"]
        expected = assemble.build_from_layout(
            "GL_VRC01", row["vh_seq"], row["vl_seq"], row["light_subclass"], parts
        )
        self.assertEqual(written[0].seq, expected.seq)

    def test_rebuilding_reports_what_it_replaced(self):
        self.client.post("/api/build", json={"ids": ["GL_VRC01"], "output": str(self.out)})
        d = self.client.post(
            "/api/build", json={"ids": ["GL_VRC01"], "output": str(self.out)}
        ).get_json()
        self.assertEqual(d["replaced"], 1)
        self.assertEqual(len(read_fasta(self.out)), 1)

    def test_build_refuses_when_nothing_is_buildable(self):
        res = self.client.post(
            "/api/build", json={"ids": ["nope"], "output": str(self.out)}
        )
        self.assertEqual(res.status_code, 400)
        self.assertFalse(self.out.exists())

    def test_plan_requires_parts_to_be_loaded(self):
        empty = tempfile.mkdtemp()
        try:
            app = create_app(empty)
            app.testing = True
            client = app.test_client()
            client.post("/api/check", json={"path": str(VRC01)})
            res = client.post("/api/plan", json={"ids": ["GL_VRC01"]})
            self.assertEqual(res.status_code, 400)
            self.assertIn("parts table is empty", res.get_json()["error"])
        finally:
            shutil.rmtree(empty, ignore_errors=True)

    def test_multiple_antibodies_plan_independently(self):
        """A kappa and a lambda in one selection must resolve separately."""
        from bcr import chains, read_fasta

        # mabA is kappa and has a complete pair in the fixture.
        self._approve(FIXTURE, ["mabA"])

        # The fixture's lambda VL has no VH partner, so pair it with a heavy
        # chain directly in the store to get a genuine lambda row.
        raw = {r.name: r.seq for r in read_fasta(FIXTURE)}
        vh = chains.process("lam_VH", raw["mabA_VH"])
        vl = chains.process("lam_VL", raw["mabC_VL"])
        self.assertEqual(vl.subclass, "lambda")
        self.store.antibodies.upsert(
            self.store.antibody_row("lam", vh, vl, FIXTURE)
        )

        d = self.client.post(
            "/api/plan", json={"ids": ["GL_VRC01", "mabA", "lam"]}
        ).get_json()
        by_id = {c["id"]: c for c in d["constructs"]}
        self.assertEqual(by_id["GL_VRC01"]["cl_part"], "IGKC")
        self.assertEqual(by_id["mabA"]["cl_part"], "IGKC")
        self.assertEqual(by_id["lam"]["cl_part"], "IGLC2")
        self.assertEqual(by_id["lam"]["subclass"], "lambda")
        self.assertEqual(d["buildable"], 3)


class ConstructLibraryTests(unittest.TestCase):
    """Card 6: built constructs are kept and can be written out again."""

    def setUp(self):
        self.dir = pathlib.Path(tempfile.mkdtemp())
        self.out = self.dir / "out.fasta"
        app = create_app(self.dir)
        app.testing = True
        self.client = app.test_client()
        self.store = Store(self.dir)
        self.client.post("/api/parts/apply", json={"path": str(PARTS)})
        self.client.post("/api/check", json={"path": str(VRC01)})
        self.client.post(
            "/api/approve",
            json={"path": str(VRC01), "ids": ["GL_VRC01"]},
        )

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def build(self):
        return self.client.post(
            "/api/build", json={"ids": ["GL_VRC01"], "output": str(self.out)}
        ).get_json()

    def test_building_stores_the_construct(self):
        d = self.build()
        self.assertEqual(len(d["library"]), 1)
        row = d["library"][0]
        self.assertEqual(row["id"], "GL_VRC01_BCR")
        self.assertEqual(row["antibody_id"], "GL_VRC01")
        self.assertEqual(row["cl_part"], "IGKC")
        self.assertEqual(row["preset"], "single_chain_bcr")
        self.assertEqual(row["has_stop"], "yes")
        self.assertEqual(int(row["residue_count"]), 806)

    def test_stored_sequence_matches_the_written_fasta(self):
        self.build()
        written = read_fasta(self.out)[0]
        row = self.store.constructs.read_map()["GL_VRC01_BCR"]
        self.assertEqual(row["seq"], written.seq)

    def test_rebuilding_updates_rather_than_duplicates(self):
        self.build()
        d = self.build()
        self.assertEqual(len(d["library"]), 1)

    def test_segment_map_round_trips_through_the_csv(self):
        import json

        self.build()
        row = self.store.constructs.read_map()["GL_VRC01_BCR"]
        segments = json.loads(row["segments_json"])
        self.assertEqual(
            [s["label"] for s in segments],
            ["SignalHC", "VL:GL_VRC01", "IGKC", "LinkerSC",
             "VH:GL_VRC01", "CHandCTSC", "Stop"],
        )
        # Coordinates must still address the stored sequence correctly.
        seq = row["seq"]
        for s in segments:
            self.assertEqual(len(seq[s["start"] - 1 : s["end"]]), s["length"])
        self.assertEqual(segments[-1]["end"], len(seq))

    def test_constructs_endpoint_lists_them(self):
        self.build()
        d = self.client.get("/api/constructs").get_json()
        self.assertEqual([c["id"] for c in d["constructs"]], ["GL_VRC01_BCR"])

    def test_delete_removes_a_construct(self):
        self.build()
        res = self.client.post(
            "/api/construct/delete", json={"id": "GL_VRC01_BCR"}
        )
        self.assertEqual(res.status_code, 200)
        self.assertEqual(res.get_json()["library"], [])

    def test_export_writes_selected_constructs(self):
        self.build()
        target = self.dir / "export.fasta"
        res = self.client.post(
            "/api/constructs/export",
            json={"ids": ["GL_VRC01_BCR"], "output": str(target)},
        )
        self.assertEqual(res.status_code, 200)
        self.assertEqual(res.get_json()["written"], 1)
        exported = read_fasta(target)
        self.assertEqual(exported[0].name, "GL_VRC01_BCR")
        self.assertEqual(exported[0].seq, read_fasta(self.out)[0].seq)

    def test_export_refuses_an_unknown_id(self):
        self.build()
        res = self.client.post(
            "/api/constructs/export",
            json={"ids": ["nope"], "output": str(self.dir / "x.fasta")},
        )
        self.assertEqual(res.status_code, 400)
        self.assertIn("Not in the library", res.get_json()["error"])

    def test_export_rejects_an_empty_selection(self):
        res = self.client.post("/api/constructs/export", json={"ids": []})
        self.assertEqual(res.status_code, 400)


if __name__ == "__main__":
    unittest.main()

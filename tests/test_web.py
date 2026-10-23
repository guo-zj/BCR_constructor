"""Tests for the CSV store and the web endpoints.

    python3 -m unittest discover -s tests -t . -v
"""

import csv
import pathlib
import shutil
import tempfile
import unittest

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


class StoreTests(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.store = Store(self.dir)

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_tables_start_empty_without_files(self):
        self.assertEqual(self.store.antibodies.read(), [])
        self.assertEqual(self.store.parts.read(), [])
        self.assertEqual(self.store.constructs.read(), [])

    def test_upsert_adds_then_updates_in_place(self):
        table = self.store.antibodies
        self.assertEqual(table.upsert({"id": "x", "vh_len": "1"}), "added")
        self.assertEqual(table.upsert({"id": "x", "vh_len": "2"}), "updated")
        rows = table.read()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["vh_len"], "2")

    def test_header_is_written_and_stable(self):
        self.store.antibodies.upsert({"id": "x"})
        with open(self.store.antibodies.path, newline="", encoding="utf-8") as fh:
            header = next(csv.reader(fh))
        self.assertEqual(header[0], "id")
        self.assertIn("vh_seq", header)
        self.assertIn("light_subclass", header)

    def test_delete_removes_only_the_named_row(self):
        self.store.antibodies.upsert({"id": "a"})
        self.store.antibodies.upsert({"id": "b"})
        self.assertTrue(self.store.antibodies.delete("a"))
        self.assertEqual([r["id"] for r in self.store.antibodies.read()], ["b"])
        self.assertFalse(self.store.antibodies.delete("missing"))

    def test_part_row_records_role_and_strips_any_stop(self):
        row = self.store.part_row("IGKC", "RTVAAPSV", "parts_sc.fasta")
        self.assertEqual(row["role"], "light_constant")
        self.assertEqual(row["subclass"], "kappa")

        # A stop codon is never stored on a part, whatever the source held.
        row = self.store.part_row("CHandCTSC", "ASTKGPSV*", "parts_sc.fasta")
        self.assertEqual(row["seq"], "ASTKGPSV")
        self.assertEqual(row["length"], 8)
        self.assertNotIn("has_stop", row)

    def test_parts_round_trip_into_assemble_shape(self):
        from bcr import assemble

        parts, problems = assemble.load_parts(PARTS)
        self.assertEqual(problems, [])
        self.store.import_parts(parts, "parts_sc.fasta")
        self.assertEqual(self.store.parts_as_dict(), parts)

    def test_no_stop_is_stored_on_any_part(self):
        from bcr import assemble

        parts, _ = assemble.load_parts(PARTS)
        self.store.import_parts(parts, "parts_sc.fasta")
        for name, seq in self.store.parts_as_dict().items():
            with self.subTest(part=name):
                self.assertNotIn("*", seq)

    def test_construct_row_keeps_segment_coordinates(self):
        import json

        from bcr import assemble, read_fasta

        parts, _ = assemble.load_parts(PARTS)
        clean = clean_vrc01()
        c = assemble.build_from_layout(
            "GL_VRC01", clean["GL_VRC01_VH"], clean["GL_VRC01_VL"], "kappa", parts
        )
        row = self.store.construct_row(c, "single_chain_bcr", "out.fasta")
        self.assertEqual(row["id"], "GL_VRC01_BCR")
        self.assertEqual(row["cl_part"], "IGKC")
        self.assertEqual(row["has_stop"], "yes")
        segments = json.loads(row["segments_json"])
        self.assertEqual(len(segments), 7)
        self.assertEqual(segments[0]["label"], "SignalHC")
        self.assertEqual(segments[-1]["end"], len(c.seq))


class CheckAndApproveTests(unittest.TestCase):
    """Card 1 is two-phase: /api/check proposes, /api/approve commits."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        app = create_app(self.dir)
        app.testing = True
        self.client = app.test_client()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def check(self, path):
        return self.client.post("/api/check", json={"path": str(path)})

    def approve(self, path, ids):
        return self.client.post(
            "/api/approve", json={"path": str(path), "ids": ids}
        )

    def library(self):
        return self.client.get("/api/library").get_json()["antibodies"]

    def test_index_renders_with_a_closed_utility_panel(self):
        res = self.client.get("/")
        self.assertEqual(res.status_code, 200)
        page = res.data
        self.assertIn(b"Constructs library", page)

        # The three input cards live in the utility panel, which starts closed
        # (the panel carries `hidden`) behind a banner along the bottom.
        self.assertIn(b'id="utilBar"', page)
        self.assertIn(b'<div id="utilPanel" hidden', page)
        for pane in (b'id="util-check"', b'id="util-parts"', b'id="util-presets"'):
            self.assertIn(pane, page)
        for title in (b"Sequence sanity check", b"Construct parts", b"Assembly presets"):
            self.assertIn(title, page)

        # Their controls must still be present, since the JS binds to them.
        for control in (b'id="seqPath"', b'id="partsBtn"', b'id="draftChain"',
                        b'id="palette"', b'id="processBtn"'):
            self.assertIn(control, page)

    def test_hidden_elements_are_not_defeated_by_a_display_rule(self):
        """An author `display` rule beats the browser's `[hidden]` default.

        Any element relying on the `hidden` attribute whose selector also sets
        `display` needs an explicit `[hidden] { display: none }` guard, or it
        renders open regardless of the attribute. This caught the utility
        panel showing permanently.
        """
        import re

        page = self.client.get("/").data.decode()
        css = page[page.index("<style>") : page.index("</style>")]
        html = page[page.index("</style>") :]

        for match in re.finditer(r"<(\w+)([^>]*\bhidden\b[^>]*)>", html):
            attrs = match.group(2)
            ident = re.search(r'id="([^"]+)"', attrs)
            classes = re.search(r'class="([^"]+)"', attrs)
            selectors = ([f"#{ident.group(1)}"] if ident else []) + [
                f".{c}" for c in (classes.group(1).split() if classes else [])
            ]
            for selector in selectors:
                sets_display = any(
                    "display:" in block.group(1)
                    for block in re.finditer(
                        re.escape(selector) + r"\s*(?:,[^{]*)?\{([^}]*)\}", css
                    )
                )
                if not sets_display:
                    continue
                with self.subTest(selector=selector):
                    self.assertIn(
                        f"{selector}[hidden]",
                        css,
                        f"{selector} sets display and is used with [hidden], "
                        f"so it needs a {selector}[hidden] rule",
                    )

    def test_input_cards_are_not_in_the_main_flow(self):
        """They moved into the panel, so the page body no longer holds them."""
        page = self.client.get("/").data.decode()
        body = page[page.index('<div class="wrap">'):page.index('id="utilBar"')]
        self.assertIn("Antibody library", body)
        self.assertNotIn('id="seqPath"', body)
        self.assertNotIn('id="draftChain"', body)

    def test_check_writes_nothing(self):
        status, data = self.check(VRC01).status_code, None
        self.assertEqual(status, 200)
        self.assertEqual(self.library(), [])

    def test_check_reports_chain_type_and_where_to_trim(self):
        data = self.check(VRC01).get_json()
        self.assertEqual(data["summary"], {"passed": 2, "failed": 0, "ready": 1})
        by_name = {r["name"]: r for r in data["records"]}

        vh = by_name["GL_VRC01_VH"]
        self.assertEqual(vh["detected"], "heavy")
        self.assertEqual(vh["trim"]["v_start"], 1)
        self.assertEqual(vh["trim"]["v_end"], 121)
        self.assertEqual(vh["trim"]["n_count"], 0)
        self.assertEqual(vh["trim"]["kept"], vh["trim"]["kept"])

        vl = by_name["GL_VRC01_VL"]
        self.assertEqual(vl["detected"], "light")
        self.assertEqual(vl["subclass"], "kappa")

    def test_trim_view_partitions_the_original_sequence(self):
        """n_seq + kept + c_seq must reconstruct the input exactly."""
        from bcr import read_fasta

        data = self.check(FIXTURE).get_json()
        raw = {r.name: r.seq for r in read_fasta(FIXTURE)}
        for rec in data["records"]:
            if rec["status"] != "ok":
                continue
            t = rec["trim"]
            with self.subTest(name=rec["name"]):
                self.assertEqual(
                    t["n_seq"] + t["kept"] + t["c_seq"], raw[rec["name"]]
                )
                self.assertEqual(len(t["n_seq"]), t["n_count"])
                self.assertEqual(len(t["c_seq"]), t["c_count"])

    def test_trim_view_shows_the_leader_and_constant_region(self):
        data = self.check(FIXTURE).get_json()
        mabB = next(r for r in data["records"] if r["name"] == "mabB_VH")
        self.assertEqual(mabB["trim"]["n_count"], 19)
        self.assertEqual(mabB["trim"]["n_seq"], "MGWSCIILFLVATATGVHS")
        self.assertTrue(mabB["trim"]["c_seq"].startswith("ASTKGPSVFPLAP"))

    def test_proposal_marks_a_new_antibody_as_add(self):
        data = self.check(VRC01).get_json()
        p = next(p for p in data["proposals"] if p["id"] == "GL_VRC01")
        self.assertTrue(p["ready"])
        self.assertEqual(p["action"], "add")

    def test_proposal_marks_an_existing_antibody_as_replace(self):
        self.approve(VRC01, ["GL_VRC01"])
        data = self.check(VRC01).get_json()
        p = next(p for p in data["proposals"] if p["id"] == "GL_VRC01")
        self.assertEqual(p["action"], "replace")
        self.assertFalse(p["changed"])

    def test_unpaired_and_failed_chains_cannot_be_added(self):
        data = self.check(FIXTURE).get_json()
        self.assertEqual(data["summary"], {"passed": 4, "failed": 2, "ready": 1})
        by_id = {p["id"]: p for p in data["proposals"]}
        self.assertTrue(by_id["mabA"]["ready"])
        self.assertFalse(by_id["mabB"]["ready"])
        self.assertIn("light chain", by_id["mabB"]["detail"])
        self.assertFalse(by_id["mabC"]["ready"])
        self.assertIn("heavy chain", by_id["mabC"]["detail"])
        self.assertFalse(by_id["mabD"]["ready"])
        self.assertIn("failed", by_id["mabD"]["detail"])

    def test_approve_commits_only_the_approved_ids(self):
        res = self.approve(FIXTURE, ["mabA"])
        self.assertEqual(res.status_code, 200)
        self.assertEqual([s["id"] for s in res.get_json()["saved"]], ["mabA"])
        self.assertEqual([r["id"] for r in self.library()], ["mabA"])

    def test_approve_refuses_an_incomplete_pair(self):
        data = self.approve(FIXTURE, ["mabB"]).get_json()
        self.assertEqual(data["saved"], [])
        self.assertEqual(data["skipped"][0]["id"], "mabB")
        self.assertEqual(self.library(), [])

    def test_approve_reports_an_id_not_in_the_file(self):
        data = self.approve(VRC01, ["nonesuch"]).get_json()
        self.assertEqual(data["saved"], [])
        self.assertIn("not found in this file", data["skipped"][0]["detail"])

    def test_approve_rejects_an_empty_selection(self):
        self.assertEqual(self.approve(VRC01, []).status_code, 400)

    def test_approved_sequences_are_the_trimmed_ones(self):
        from bcr import chains, read_fasta

        self.approve(VRC01, ["GL_VRC01"])
        row = self.library()[0]
        raw = {r.name: r.seq for r in read_fasta(VRC01)}
        self.assertEqual(
            row["vh_seq"], chains.process("x_VH", raw["GL_VRC01_VH"]).clean_seq
        )
        self.assertEqual(
            row["vl_seq"], chains.process("x_VL", raw["GL_VRC01_VL"]).clean_seq
        )
        self.assertEqual(row["light_subclass"], "kappa")

    def test_reapproving_updates_rather_than_duplicates(self):
        self.approve(VRC01, ["GL_VRC01"])
        data = self.approve(VRC01, ["GL_VRC01"]).get_json()
        self.assertEqual(data["saved"][0]["outcome"], "updated")
        self.assertEqual(len(self.library()), 1)

    def test_typing_warnings_are_carried_into_the_row(self):
        self.approve(VRC01, ["GL_VRC01"])
        notes = self.library()[0]["notes"]
        self.assertIn("VL:", notes)
        self.assertIn("103 aa", notes)

    def test_bad_paths_are_rejected(self):
        for path in ("", "/nope/missing.fasta", str(ROOT)):
            with self.subTest(path=path):
                res = self.client.post("/api/check", json={"path": path})
                self.assertEqual(res.status_code, 400)
                self.assertTrue(res.get_json()["error"])

    def test_delete_endpoint(self):
        self.approve(VRC01, ["GL_VRC01"])
        res = self.client.post("/api/antibody/delete", json={"id": "GL_VRC01"})
        self.assertEqual(res.status_code, 200)
        self.assertEqual(res.get_json()["library"], [])


class PartsImportTests(unittest.TestCase):
    """Card 3 is two-phase: preview compares, apply honours a decision."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        app = create_app(self.dir)
        app.testing = True
        self.client = app.test_client()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def preview(self, path):
        return self.client.post("/api/parts/preview", json={"path": str(path)})

    def apply(self, path, decisions=None):
        return self.client.post(
            "/api/parts/apply",
            json={"path": str(path), "decisions": decisions or {}},
        )

    def altered(self, **replacements):
        """A copy of misc.fasta with substrings swapped."""
        text = (PARTS).read_text(encoding="utf-8")
        for old, new in replacements.items():
            text = text.replace(old, new)
        path = pathlib.Path(self.dir) / "alt.fasta"
        path.write_text(text, encoding="utf-8")
        return path

    def test_preview_writes_nothing(self):
        self.assertEqual(self.preview(PARTS).status_code, 200)
        self.assertEqual(
            self.client.get("/api/library").get_json()["parts"], []
        )

    def test_everything_is_new_on_an_empty_table(self):
        data = self.preview(PARTS).get_json()
        self.assertEqual(data["summary"], {"new": 5, "identical": 0, "differs": 0})
        self.assertTrue(all(p["status"] == "new" for p in data["parts"]))

    def test_reimporting_the_same_file_is_all_identical(self):
        self.apply(PARTS)
        data = self.preview(PARTS).get_json()
        self.assertEqual(data["summary"], {"new": 0, "identical": 5, "differs": 0})

    def test_a_changed_part_is_detected_with_its_first_difference(self):
        self.apply(PARTS)
        alt = self.altered(**{"MGWSCIILFLVATATGVHS": "MGWSCIILFLVATATGVHQ"})
        data = self.preview(alt).get_json()
        self.assertEqual(data["summary"]["differs"], 1)
        signal = next(p for p in data["parts"] if p["name"] == "SignalHC")
        self.assertEqual(signal["status"], "differs")
        self.assertEqual(signal["first_diff"], 19)
        self.assertEqual(signal["existing_len"], 19)

    def test_a_length_change_is_described(self):
        self.apply(PARTS)
        alt = self.altered(**{"MGWSCIILFLVATATGVHS": "MGWSCIILFLVATATGVHSAA"})
        signal = next(
            p for p in self.preview(alt).get_json()["parts"] if p["name"] == "SignalHC"
        )
        self.assertIn("2 residue(s) longer", signal["note"])

    def test_default_decision_keeps_the_stored_sequence(self):
        self.apply(PARTS)
        alt = self.altered(**{"MGWSCIILFLVATATGVHS": "MGWSCIILFLVATATGVHQ"})
        data = self.apply(alt).get_json()  # no decisions given
        outcomes = {a["name"]: a["outcome"] for a in data["applied"]}
        self.assertEqual(outcomes["SignalHC"], "kept original")
        stored = {p["name"]: p["seq"] for p in data["parts"]}
        self.assertTrue(stored["SignalHC"].endswith("VHS"))

    def test_override_decision_replaces_only_that_part(self):
        self.apply(PARTS)
        alt = self.altered(
            **{
                "MGWSCIILFLVATATGVHS": "MGWSCIILFLVATATGVHQ",
                "GGSSGSGSGSNWSHPQFEK": "GGSSGSGSGSNWSHPQFEQ",
            }
        )
        data = self.apply(alt, {"LinkerSC": "override"}).get_json()
        outcomes = {a["name"]: a["outcome"] for a in data["applied"]}
        self.assertEqual(outcomes["LinkerSC"], "overridden")
        self.assertEqual(outcomes["SignalHC"], "kept original")
        self.assertEqual(data["changed"], 1)

        stored = {p["name"]: p["seq"] for p in data["parts"]}
        self.assertIn("QFEQ", stored["LinkerSC"])
        self.assertTrue(stored["SignalHC"].endswith("VHS"))

    def test_a_new_part_is_added_without_a_decision(self):
        self.apply(PARTS)
        self.client.post("/api/construct/delete", json={"id": "x"})
        # Remove one part, then re-import: it should come back as "added".
        store = Store(self.dir)
        store.parts.delete("LinkerSC")
        data = self.apply(PARTS).get_json()
        outcomes = {a["name"]: a["outcome"] for a in data["applied"]}
        self.assertEqual(outcomes["LinkerSC"], "added")

    def test_an_invalid_parts_file_is_refused_at_preview(self):
        bad = pathlib.Path(self.dir) / "bad.fasta"
        bad.write_text(">SignalHC\nMGWSCIILFLVATATGVHS\n", encoding="utf-8")
        res = self.preview(bad)
        self.assertEqual(res.status_code, 400)
        self.assertTrue(res.get_json()["problems"])

    def test_a_part_not_starting_at_its_first_residue_is_refused(self):
        """Section 6: a constant region starts at its own first residue."""
        alt = self.altered(**{"RTVAAPSVFIFPPS": "TVAAPSVFIFPPS"})
        res = self.preview(alt)
        self.assertEqual(res.status_code, 400)
        self.assertIn("IGKC", res.get_json()["error"])


if __name__ == "__main__":
    unittest.main()

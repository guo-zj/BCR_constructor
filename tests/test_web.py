"""Tests for the CSV store and the web endpoints.

    python3 -m unittest discover -s tests -t . -v
"""

import csv
import io
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
RESERVED = pathlib.Path(__file__).parent / "fixtures" / "parts_reserved.fasta"
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
        # The layout is part of the name: one antibody can be built several
        # ways into one file.
        self.assertEqual(row["id"], "GL_VRC01_BCR_single_chain_bcr")
        self.assertEqual(row["cl_part"], "IGKC")
        self.assertEqual(row["has_stop"], "yes")
        segments = json.loads(row["segments_json"])
        self.assertEqual(len(segments), 5)
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

    def test_every_table_has_an_editor_slot_outside_the_table(self):
        """A table re-render must not be able to wipe a form in progress.

        Each editor lives in its own element next to the table rather than
        inside it, so re-rendering rows leaves a half-typed record alone.
        """
        page = self.client.get("/").data.decode()
        for slot, table in (('id="abEdit"', 'id="abOut"'),
                            ('id="partEdit"', 'id="partsOut"'),
                            ('id="conEdit"', 'id="conOut"')):
            with self.subTest(slot=slot):
                self.assertIn(slot, page)
                self.assertLess(page.index(slot), page.index(table))
        for control in (b'id="abAdd"', b'id="partAdd"'):
            self.assertIn(control, self.client.get("/").data)

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

    def test_the_scripts_block_comments_are_balanced(self):
        """A stray `*/` silently kills the whole script, and nothing else here
        would notice: there is no JS engine on this machine, so the page is
        never parsed. Markers must strictly alternate `/* … */`, which catches
        both an early close left behind by an edit and an unterminated comment.
        Nothing in the script holds either marker inside a string, and this
        test is the reason to keep it that way.
        """
        import re

        page = self.client.get("/").data.decode()
        js = page[page.rindex("<script>") + len("<script>"):page.rindex("</script>")]
        depth = 0
        for match in re.finditer(r"/\*|\*/", js):
            line = js.count("\n", 0, match.start()) + 1
            if match.group() == "/*":
                self.assertEqual(depth, 0, f"nested /* in the script, line {line}")
                depth = 1
            else:
                self.assertEqual(depth, 1, f"unmatched */ in the script, line {line}")
                depth = 0
        self.assertEqual(depth, 0, "the script ends inside a block comment")

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
        self.assertEqual(data["summary"],
                         {"passed": 2, "failed": 0, "as_given": 0, "ready": 1, "awaiting": 0})
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
        self.assertEqual(data["summary"],
                         {"passed": 4, "failed": 2, "as_given": 0, "ready": 1, "awaiting": 0})
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
        # The intrinsic parts are always listed; nothing of the user's is.
        rows = self.client.get("/api/library").get_json()["parts"]
        self.assertEqual([r["name"] for r in rows if not r["builtin"]], [])

    def test_everything_is_new_on_an_empty_table(self):
        data = self.preview(PARTS).get_json()
        self.assertEqual(
            data["summary"],
            {"new": 8, "identical": 0, "differs": 0, "reserved": 0},
        )
        self.assertTrue(all(p["status"] == "new" for p in data["parts"]))

    def test_reimporting_the_same_file_is_all_identical(self):
        self.apply(PARTS)
        data = self.preview(PARTS).get_json()
        self.assertEqual(
            data["summary"],
            {"new": 0, "identical": 8, "differs": 0, "reserved": 0},
        )

    def test_a_changed_part_is_detected_with_its_first_difference(self):
        self.apply(PARTS)
        alt = self.altered(**{"MEWSWVFLFFLSVTTGVHS": "MEWSWVFLFFLSVTTGVHQ"})
        data = self.preview(alt).get_json()
        self.assertEqual(data["summary"]["differs"], 1)
        signal = next(p for p in data["parts"] if p["name"] == "SignalHC")
        self.assertEqual(signal["status"], "differs")
        self.assertEqual(signal["first_diff"], 19)
        self.assertEqual(signal["existing_len"], 19)

    def test_a_length_change_is_described(self):
        self.apply(PARTS)
        alt = self.altered(**{"MEWSWVFLFFLSVTTGVHS": "MEWSWVFLFFLSVTTGVHSAA"})
        signal = next(
            p for p in self.preview(alt).get_json()["parts"] if p["name"] == "SignalHC"
        )
        self.assertIn("2 residue(s) longer", signal["note"])

    def test_default_decision_keeps_the_stored_sequence(self):
        self.apply(PARTS)
        alt = self.altered(**{"MEWSWVFLFFLSVTTGVHS": "MEWSWVFLFFLSVTTGVHQ"})
        data = self.apply(alt).get_json()  # no decisions given
        outcomes = {a["name"]: a["outcome"] for a in data["applied"]}
        self.assertEqual(outcomes["SignalHC"], "kept original")
        stored = {p["name"]: p["seq"] for p in data["parts"]}
        self.assertTrue(stored["SignalHC"].endswith("TGVHS"))

    def test_override_decision_replaces_only_that_part(self):
        self.apply(PARTS)
        alt = self.altered(
            **{
                "MEWSWVFLFFLSVTTGVHS": "MEWSWVFLFFLSVTTGVHQ",
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
        self.assertTrue(stored["SignalHC"].endswith("TGVHS"))

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
        bad.write_text(">MyPart\nGGGS123OO!\n", encoding="utf-8")
        res = self.preview(bad)
        self.assertEqual(res.status_code, 400)
        self.assertTrue(res.get_json()["problems"])

    def test_a_file_asking_for_nothing_in_particular_is_fine(self):
        """Nothing is demanded of a parts file: the fixed names are intrinsic."""
        thin = pathlib.Path(self.dir) / "thin.fasta"
        thin.write_text(">MyLinker\nGGGGSGGGGS\n", encoding="utf-8")
        res = self.preview(thin)
        self.assertEqual(res.status_code, 200)
        self.assertEqual(res.get_json()["summary"]["new"], 1)

    def test_a_reserved_name_is_offered_a_rename_rather_than_imported(self):
        data = self.preview(RESERVED).get_json()
        # The ordinary records come through; the reserved ones are set apart.
        self.assertEqual([p["name"] for p in data["parts"]], ["SignalHC"])
        by_name = {r["name"]: r for r in data["reserved"]}
        self.assertEqual(sorted(by_name), ["IGKC", "IGLC2"])
        # IGLC2 matches the intrinsic sequence, so it needs no rename.
        self.assertTrue(by_name["IGLC2"]["identical"])
        self.assertEqual(by_name["IGLC2"]["suggestion"], "")
        # IGKC differs, so it is the one that has to be renamed.
        self.assertFalse(by_name["IGKC"]["identical"])
        self.assertEqual(by_name["IGKC"]["suggestion"], "IGKC_mine")
        self.assertEqual(data["summary"]["reserved"], 1)

    def test_applying_without_a_rename_skips_the_reserved_record(self):
        data = self.apply(RESERVED).get_json()
        outcomes = {a["name"]: a["outcome"] for a in data["applied"]}
        self.assertEqual(outcomes["IGKC"], "skipped")
        self.assertEqual(outcomes["IGLC2"], "already intrinsic")
        self.assertEqual(outcomes["SignalHC"], "added")
        stored = {p["name"]: p for p in data["parts"]}
        # The intrinsic IGKC is untouched by the attempt.
        self.assertEqual(stored["IGKC"]["builtin"], "yes")
        self.assertEqual(stored["IGKC"]["seq"], assemble.INTRINSIC_SEQS["IGKC"])

    def test_a_rename_imports_the_sequence_under_the_new_name(self):
        res = self.client.post(
            "/api/parts/apply",
            json={"path": str(RESERVED), "renames": {"IGKC": "IGKC_Km3"}},
        )
        self.assertEqual(res.status_code, 200)
        data = res.get_json()
        outcomes = {a["name"]: a["outcome"] for a in data["applied"]}
        self.assertEqual(outcomes["IGKC"], "renamed")
        stored = {p["name"]: p for p in data["parts"]}
        self.assertIn("IGKC_Km3", stored)
        self.assertEqual(stored["IGKC_Km3"]["builtin"], "")
        self.assertEqual(stored["IGKC"]["seq"], assemble.INTRINSIC_SEQS["IGKC"])

    def test_a_rename_onto_another_reserved_name_is_refused(self):
        data = self.client.post(
            "/api/parts/apply",
            json={"path": str(RESERVED), "renames": {"IGKC": "IGLC2"}},
        ).get_json()
        entry = next(a for a in data["applied"] if a["name"] == "IGKC")
        self.assertEqual(entry["outcome"], "skipped")
        self.assertIn("intrinsic", entry["detail"])

    def test_a_rename_to_an_unusable_name_is_refused(self):
        data = self.client.post(
            "/api/parts/apply",
            json={"path": str(RESERVED), "renames": {"IGKC": "two words"}},
        ).get_json()
        entry = next(a for a in data["applied"] if a["name"] == "IGKC")
        self.assertEqual(entry["outcome"], "skipped")


class TrimDecisionTests(unittest.TestCase):
    """A trim is a modification, so approving a record does not imply it.

    Every record the rules would shorten needs an explicit trim/keep decision;
    without one the pair is skipped and reported as waiting, never stored
    under a guess either way.
    """

    LEADER = "MGWSCIILFLVATATGVHS"

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        app = create_app(self.dir)
        app.testing = True
        self.client = app.test_client()
        self.clean = clean_vrc01()

        # The VRC01 pair with a signal peptide on the heavy chain, so exactly
        # one of the two records carries a proposed trim.
        self.path = pathlib.Path(self.dir) / "leader.fasta"
        self.path.write_text(
            f">GL_VRC01_VH\n{self.LEADER}{self.clean['GL_VRC01_VH']}\n"
            f">GL_VRC01_VL\n{self.clean['GL_VRC01_VL']}\n",
            encoding="utf-8",
        )

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def check(self):
        return self.client.post("/api/check", json={"path": str(self.path)}).get_json()

    def approve(self, decisions=None):
        body = {"path": str(self.path), "ids": ["GL_VRC01"]}
        if decisions is not None:
            body["decisions"] = decisions
        return self.client.post("/api/approve", json=body)

    def library(self):
        return self.client.get("/api/library").get_json()["antibodies"]

    def test_check_flags_which_records_would_be_modified(self):
        data = self.check()
        by_name = {r["name"]: r for r in data["records"]}
        self.assertTrue(by_name["GL_VRC01_VH"]["proposes_trim"])
        self.assertFalse(by_name["GL_VRC01_VL"]["proposes_trim"])
        self.assertEqual(data["summary"]["awaiting"], 1)

        (decision,) = data["proposals"][0]["decisions"]
        self.assertEqual(decision["record"], "GL_VRC01_VH")
        self.assertEqual(decision["chain"], "VH")
        self.assertEqual(decision["n"], len(self.LEADER))
        self.assertEqual(
            decision["as_given"] - decision["trimmed_to"], len(self.LEADER)
        )

    def test_approving_without_a_decision_stores_nothing(self):
        res = self.approve()
        self.assertEqual(res.status_code, 200)
        body = res.get_json()
        self.assertEqual(body["saved"], [])
        (skip,) = body["skipped"]
        self.assertTrue(skip["awaiting"])
        self.assertIn("GL_VRC01_VH", skip["detail"])
        self.assertEqual(self.library(), [])

    def test_choosing_trim_stores_the_v_domain(self):
        body = self.approve({"GL_VRC01_VH": "trim"}).get_json()
        self.assertEqual(body["saved"][0]["kept"], [])
        row = self.library()[0]
        self.assertEqual(row["vh_seq"], self.clean["GL_VRC01_VH"])
        self.assertNotIn("kept as given", row["vh_trim"])

    def test_choosing_keep_stores_the_sequence_untouched(self):
        body = self.approve({"GL_VRC01_VH": "keep"}).get_json()
        self.assertEqual(body["saved"][0]["kept"], ["VH"])
        row = self.library()[0]
        self.assertEqual(row["vh_seq"], self.LEADER + self.clean["GL_VRC01_VH"])
        self.assertIn("kept as given", row["vh_trim"])
        self.assertIn("not a bare V domain", row["notes"])
        # The chain that needed no decision is unaffected either way.
        self.assertEqual(row["vl_seq"], self.clean["GL_VRC01_VL"])

    def test_a_record_needing_no_trim_needs_no_decision(self):
        res = self.client.post(
            "/api/approve", json={"path": str(VRC01), "ids": ["GL_VRC01"]}
        )
        self.assertEqual([s["id"] for s in res.get_json()["saved"]], ["GL_VRC01"])

    def test_a_decision_must_be_trim_or_keep(self):
        res = self.approve({"GL_VRC01_VH": "whatever"})
        self.assertEqual(res.status_code, 400)
        self.assertEqual(self.library(), [])


class EditAntibodyTests(unittest.TestCase):
    """Hand-entered antibodies go through the same rules as imported ones."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        app = create_app(self.dir)
        app.testing = True
        self.client = app.test_client()
        self.clean = clean_vrc01()
        self.vh = self.clean["GL_VRC01_VH"]
        self.vl = self.clean["GL_VRC01_VL"]

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def save(self, **fields):
        return self.client.post("/api/antibody/save", json=fields)

    def rows(self):
        return self.client.get("/api/library").get_json()["antibodies"]

    def test_an_antibody_can_be_added_by_hand(self):
        res = self.save(id="HandMade", vh_seq=self.vh, vl_seq=self.vl,
                        user_notes="typed in from a paper")
        self.assertEqual(res.status_code, 200)
        body = res.get_json()
        self.assertEqual(body["outcome"], "added")

        (row,) = self.rows()
        self.assertEqual(row["id"], "HandMade")
        self.assertEqual(row["vh_seq"], self.vh)
        self.assertEqual(row["light_subclass"], "kappa")
        self.assertEqual(row["source_file"], "manual entry")
        # The user's note sits beside the derived warnings, not on top of them.
        self.assertEqual(row["user_notes"], "typed in from a paper")

    def test_a_whole_fasta_record_can_be_pasted_into_a_field(self):
        pasted = f">GL_VRC01_VH some description\n{self.vh[:40]}\n{self.vh[40:]}\n"
        res = self.save(id="Pasted", vh_seq=pasted, vl_seq=self.vl)
        self.assertEqual(res.status_code, 200)
        self.assertEqual(self.rows()[0]["vh_seq"], self.vh)

    def test_a_trim_is_offered_rather_than_applied(self):
        """A pasted sequence is never shortened on the strength of the paste."""
        raw = "MGWSCIILFLVATATGVHS" + self.vh
        res = self.save(id="Raw", vh_seq=raw, vl_seq=self.vl)
        self.assertEqual(res.status_code, 200)
        body = res.get_json()
        self.assertFalse(body["stored"])
        self.assertEqual(self.rows(), [])          # nothing written

        (pending,) = body["needs_decision"]
        self.assertEqual(pending["chain"], "VH")
        self.assertEqual(pending["n"], 19)
        self.assertEqual(pending["n_seq"], "MGWSCIILFLVATATGVHS")

        res = self.save(id="Raw", vh_seq=raw, vl_seq=self.vl, vh_decision="trim")
        self.assertEqual(res.status_code, 200)
        self.assertEqual(self.rows()[0]["vh_seq"], self.vh)
        self.assertTrue(any("trimmed" in n for n in res.get_json()["notes"]))

    def test_a_declined_trim_stores_the_sequence_as_given(self):
        raw = "MGWSCIILFLVATATGVHS" + self.vh
        res = self.save(id="Raw", vh_seq=raw, vl_seq=self.vl, vh_decision="keep")
        self.assertEqual(res.status_code, 200)
        row = self.rows()[0]
        self.assertEqual(row["vh_seq"], raw)
        self.assertEqual(row["vh_len"], str(len(raw)))
        self.assertIn("kept as given", row["vh_trim"])
        # FR4 still comes from the detected V end, not the stored tail.
        self.assertEqual(row["vh_fr4"], self.vh[-11:])
        self.assertIn("not a bare V domain", row["notes"])

    def test_a_sequence_that_does_not_type_is_refused(self):
        res = self.save(id="Nope", vh_seq="MKKKWWWQQQ", vl_seq=self.vl)
        self.assertEqual(res.status_code, 400)
        self.assertTrue(res.get_json()["problems"])
        self.assertEqual(self.rows(), [])

    def test_the_id_must_not_carry_a_chain_suffix(self):
        res = self.save(id="HandMade_VH", vh_seq=self.vh, vl_seq=self.vl)
        self.assertEqual(res.status_code, 400)
        self.assertIn("_VH", res.get_json()["error"])

    def test_a_name_that_would_break_a_fasta_header_is_refused(self):
        for bad in ("two words", "has,comma", "", "-leading"):
            with self.subTest(name=bad):
                res = self.save(id=bad, vh_seq=self.vh, vl_seq=self.vl)
                self.assertEqual(res.status_code, 400)

    def test_both_sequences_are_required(self):
        res = self.save(id="Half", vh_seq=self.vh, vl_seq="")
        self.assertEqual(res.status_code, 400)
        self.assertIn("VL", res.get_json()["error"])

    def test_renaming_keeps_one_row_and_repoints_its_constructs(self):
        self.save(id="Before", vh_seq=self.vh, vl_seq=self.vl)
        self.client.post("/api/plan", json={"ids": ["Before"]})
        self.client.post("/api/parts/import", json={"path": str(PARTS)})
        build = self.client.post(
            "/api/build",
            json={"ids": ["Before"], "output": str(pathlib.Path(self.dir) / "out.fasta")},
        )
        self.assertEqual(build.status_code, 200)

        res = self.save(old_id="Before", id="After", vh_seq=self.vh, vl_seq=self.vl)
        self.assertEqual(res.status_code, 200)
        body = res.get_json()
        self.assertEqual(body["outcome"], "renamed")
        self.assertEqual([r["id"] for r in self.rows()], ["After"])
        # A construct records what it was built from, so the reference moves too.
        self.assertEqual(
            [c["antibody_id"] for c in body["constructs"]], ["After"]
        )
        self.assertTrue(any("repointed" in n for n in body["notes"]))

    def test_renaming_onto_an_existing_id_is_refused(self):
        self.save(id="One", vh_seq=self.vh, vl_seq=self.vl)
        self.save(id="Two", vh_seq=self.vh, vl_seq=self.vl)
        res = self.save(old_id="One", id="Two", vh_seq=self.vh, vl_seq=self.vl)
        self.assertEqual(res.status_code, 409)
        self.assertEqual(len(self.rows()), 2)

    def test_editing_a_row_that_is_gone_is_refused(self):
        res = self.save(old_id="Ghost", id="Ghost", vh_seq=self.vh, vl_seq=self.vl)
        self.assertEqual(res.status_code, 404)


class EditPartTests(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        app = create_app(self.dir)
        app.testing = True
        self.client = app.test_client()
        self.client.post("/api/parts/import", json={"path": str(PARTS)})

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def save(self, **fields):
        return self.client.post("/api/part/save", json=fields)

    def parts(self):
        return {r["name"]: r for r in self.client.get("/api/library").get_json()["parts"]}

    def test_a_part_can_be_added_by_hand_with_its_role_inferred(self):
        res = self.save(name="FlagTag", seq="GGGSGGGSDYKDDDDK*",
                        user_notes="C-terminal tag")
        self.assertEqual(res.status_code, 200)
        body = res.get_json()
        self.assertEqual(body["outcome"], "added")

        row = self.parts()["FlagTag"]
        # The trailing stop is dropped: parts are residues only.
        self.assertEqual(row["seq"], "GGGSGGGSDYKDDDDK")
        self.assertEqual(row["length"], "16")
        self.assertEqual(row["role"], "other")
        self.assertEqual(row["user_notes"], "C-terminal tag")

    def test_a_hand_typed_signal_peptide_is_recognised(self):
        self.save(name="SignalX", seq="MGWSCIILFLVATATGVHS")
        self.assertEqual(self.parts()["SignalX"]["role"], "signal_peptide")

    def test_a_part_with_an_internal_stop_is_refused(self):
        res = self.save(name="Broken", seq="GGGS*GGGS")
        self.assertEqual(res.status_code, 400)
        self.assertIn("truncated", res.get_json()["error"])
        self.assertNotIn("Broken", self.parts())

    def test_a_part_with_non_amino_acid_characters_is_refused(self):
        res = self.save(name="Junk", seq="GGGS123OO!")
        self.assertEqual(res.status_code, 400)
        self.assertIn("non-amino-acid", res.get_json()["error"])

    def test_an_intrinsic_part_cannot_be_saved_over(self):
        """They ship with the package, so there is nothing here to edit."""
        for name in sorted(assemble.RESERVED_PART_NAMES):
            with self.subTest(name=name):
                res = self.save(name=name, seq="GQPKAAPSVTLFPPSSEELQ")
                self.assertEqual(res.status_code, 400)
                self.assertIn("intrinsic", res.get_json()["error"])
                self.assertEqual(
                    res.get_json()["reserved"],
                    sorted(assemble.RESERVED_PART_NAMES),
                )

    def test_section_six_still_governs_the_intrinsic_constants(self):
        """Their own starting residues are checked where they are loaded."""
        _, problems = assemble.validate_part("IGLC2", "QPKAAPSVTLFPPSSEELQ")
        self.assertTrue(any("GQPKAAPSV" in p for p in problems), problems)
        _, problems = assemble.validate_part("CHTMCTIgM", "SASAPTLFPLVSC")
        self.assertTrue(any("GSASAPTLFPLVSC" in p for p in problems), problems)

    def test_renaming_a_part_follows_through_into_saved_layouts(self):
        slots = [
            {"kind": "part", "part": "SignalHC"},
            {"kind": "vl"}, {"kind": "light_constant"},
            {"kind": "vh"}, {"kind": "stop"},
        ]
        saved = self.client.post(
            "/api/presets/save",
            json={"label": "Mine", "description": "", "slots": slots},
        )
        self.assertEqual(saved.status_code, 200)

        res = self.save(old_name="SignalHC", name="SignalHC2",
                        seq=self.parts()["SignalHC"]["seq"])
        self.assertEqual(res.status_code, 200)
        self.assertEqual(res.get_json()["outcome"], "renamed")
        self.assertNotIn("SignalHC", self.parts())

        presets = self.client.get("/api/presets").get_json()["presets"]
        mine = next(p for p in presets.values() if p["label"] == "Mine")
        self.assertEqual(
            [s.get("part") for s in mine["slots"] if s["kind"] == "part"],
            ["SignalHC2"],
        )
        # The built-in layout still asks for the old name, and says so.
        self.assertTrue(any("built-in" in n for n in res.get_json()["notes"]))

    def test_the_intrinsic_parts_cannot_be_renamed_or_deleted(self):
        """Layouts resolve to them by name, so the names are fixed."""
        for name in sorted(assemble.RESERVED_PART_NAMES):
            with self.subTest(name=name):
                res = self.save(old_name=name, name=f"{name}_v2",
                                seq=self.parts()[name]["seq"])
                self.assertEqual(res.status_code, 400)
                self.assertIn(name, self.parts())

                res = self.client.post("/api/part/delete", json={"name": name})
                self.assertEqual(res.status_code, 400)
                self.assertIn(name, self.parts())

    def test_a_user_part_can_be_saved_over_by_editing_it(self):
        row = self.parts()["LinkerSC"]
        # Adding it as if it were new collides, which is the right answer:
        # it says to edit the stored one rather than silently replacing it.
        self.assertEqual(self.save(name="LinkerSC", seq=row["seq"]).status_code, 409)

        res = self.save(old_name="LinkerSC", name="LinkerSC", seq=row["seq"],
                        user_notes="checked against the plasmid map")
        self.assertEqual(res.status_code, 200)
        self.assertEqual(res.get_json()["outcome"], "updated")
        self.assertEqual(
            self.parts()["LinkerSC"]["user_notes"], "checked against the plasmid map"
        )

    def test_deleting_a_part_a_layout_needs_says_what_will_break(self):
        res = self.client.post("/api/part/delete", json={"name": "SignalHC"})
        self.assertEqual(res.status_code, 200)
        body = res.get_json()
        self.assertTrue(body["ok"])
        self.assertNotIn("SignalHC", {r["name"] for r in body["parts"]})
        self.assertTrue(any("missing part" in n for n in body["notes"]))

    def test_deleting_an_unused_part_is_quiet(self):
        self.save(name="Spare", seq="GGGSGGGS")
        res = self.client.post("/api/part/delete", json={"name": "Spare"})
        self.assertEqual(res.get_json()["notes"], [])


class EditConstructTests(unittest.TestCase):
    """A construct is derived, so only its id and notes are editable."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        app = create_app(self.dir)
        app.testing = True
        self.client = app.test_client()
        clean = clean_vrc01()
        self.client.post("/api/parts/import", json={"path": str(PARTS)})
        self.client.post(
            "/api/antibody/save",
            json={"id": "GL_VRC01", "vh_seq": clean["GL_VRC01_VH"],
                  "vl_seq": clean["GL_VRC01_VL"]},
        )
        self.client.post(
            "/api/build",
            json={"ids": ["GL_VRC01"],
                  "output": str(pathlib.Path(self.dir) / "out.fasta")},
        )

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def rows(self):
        return {r["id"]: r for r in self.client.get("/api/constructs").get_json()["constructs"]}

    def test_a_construct_can_be_renamed_and_annotated(self):
        built_id = Store.construct_id("GL_VRC01", assemble.DEFAULT_PRESET)
        before = self.rows()[built_id]
        res = self.client.post(
            "/api/construct/save",
            json={"old_id": built_id, "id": "VRC01_scBCR_v1",
                  "user_notes": "ordered 2026-10-23"},
        )
        self.assertEqual(res.status_code, 200)
        rows = self.rows()
        self.assertNotIn(built_id, rows)
        row = rows["VRC01_scBCR_v1"]
        self.assertEqual(row["user_notes"], "ordered 2026-10-23")
        # Everything derived from the build survives the rename untouched.
        for field in ("seq", "layout", "segments_json", "residue_count", "preset"):
            self.assertEqual(row[field], before[field])

    def test_a_construct_cannot_be_typed_in_from_scratch(self):
        res = self.client.post(
            "/api/construct/save", json={"id": "Invented", "user_notes": ""}
        )
        self.assertEqual(res.status_code, 400)
        self.assertIn("Assembly layout", res.get_json()["error"])

    def test_renaming_onto_an_existing_construct_is_refused(self):
        self.client.post(
            "/api/construct/save",
            json={"old_id": Store.construct_id("GL_VRC01", assemble.DEFAULT_PRESET),
                  "id": "Taken", "user_notes": ""},
        )
        self.client.post(
            "/api/antibody/save",
            json={"id": "Other", "vh_seq": clean_vrc01()["GL_VRC01_VH"],
                  "vl_seq": clean_vrc01()["GL_VRC01_VL"]},
        )
        self.client.post(
            "/api/build",
            json={"ids": ["Other"],
                  "output": str(pathlib.Path(self.dir) / "out2.fasta")},
        )
        res = self.client.post(
            "/api/construct/save",
            json={"old_id": Store.construct_id("Other", assemble.DEFAULT_PRESET),
                  "id": "Taken", "user_notes": ""},
        )
        self.assertEqual(res.status_code, 409)


class DerivedChainTests(unittest.TestCase):
    """Express mode's three columns: whole chains from intrinsic parts."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        app = create_app(self.dir)
        app.testing = True
        self.client = app.test_client()
        self.clean = clean_vrc01()
        self.client.post(
            "/api/antibody/save",
            json={"id": "GL_VRC01", "vh_seq": self.clean["GL_VRC01_VH"],
                  "vl_seq": self.clean["GL_VRC01_VL"]},
        )

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def row(self):
        return self.client.get("/api/library").get_json()["antibodies"][0]

    def test_the_three_columns_are_filled_on_every_write(self):
        row = self.row()
        vh, vl = self.clean["GL_VRC01_VH"], self.clean["GL_VRC01_VL"]
        self.assertEqual(row["light_chain"], vl + assemble.INTRINSIC_SEQS["IGKC"])
        self.assertEqual(
            row["heavy_chain_igg1"], vh + assemble.INTRINSIC_SEQS["CHTMCTIgG1"]
        )
        self.assertEqual(
            row["heavy_chain_igm"], vh + assemble.INTRINSIC_SEQS["CHTMCTIgM"]
        )

    def test_a_lambda_antibody_takes_iglc2(self):
        """The light chain follows the antibody's own call, per record."""
        from bcr import chains

        vl = chains.process("x_VL", self.clean["GL_VRC01_VL"]).clean_seq
        lam = vl[:-1] + "L"
        self.client.post(
            "/api/antibody/save",
            json={"id": "Lam", "vh_seq": self.clean["GL_VRC01_VH"], "vl_seq": lam},
        )
        rows = {r["id"]: r for r in self.client.get("/api/library").get_json()["antibodies"]}
        self.assertEqual(rows["Lam"]["light_subclass"], "lambda")
        self.assertTrue(
            rows["Lam"]["light_chain"].endswith(assemble.INTRINSIC_SEQS["IGLC2"])
        )
        self.assertTrue(
            rows["GL_VRC01"]["light_chain"].endswith(assemble.INTRINSIC_SEQS["IGKC"])
        )

    def test_a_fresh_row_is_not_stale(self):
        self.assertEqual(self.client.get("/api/library").get_json()["stale_derived"], [])

    def test_stale_columns_are_reported_rather_than_rewritten_on_read(self):
        """A read must not write: the recompute is offered instead."""
        store = Store(self.dir)
        rows = store.antibodies.read()
        rows[0]["light_chain"] = "WRONG"
        store.antibodies.write_all(rows)

        data = self.client.get("/api/library").get_json()
        self.assertEqual(data["stale_derived"], ["GL_VRC01"])
        # Still wrong on disk: reading reported it, nothing more.
        self.assertEqual(Store(self.dir).antibodies.read()[0]["light_chain"], "WRONG")

        res = self.client.post("/api/library/derive", json={})
        self.assertEqual(res.get_json()["changed"], ["GL_VRC01"])
        self.assertEqual(res.get_json()["stale_derived"], [])
        self.assertEqual(
            self.row()["light_chain"],
            self.clean["GL_VRC01_VL"] + assemble.INTRINSIC_SEQS["IGKC"],
        )

    def test_recomputing_when_nothing_is_stale_writes_nothing(self):
        before = (Store(self.dir).antibodies.path).stat().st_mtime_ns
        res = self.client.post("/api/library/derive", json={})
        self.assertEqual(res.get_json()["changed"], [])
        self.assertEqual(
            (Store(self.dir).antibodies.path).stat().st_mtime_ns, before
        )

    def test_no_subclass_means_no_light_chain_rather_than_a_guess(self):
        chains = assemble.derived_chains("QVQLV", "DIQMT", None)
        self.assertEqual(chains["light_chain"], "")
        self.assertTrue(chains["heavy_chain_igg1"].startswith("QVQLV"))


class ExpressLayoutTests(unittest.TestCase):
    """The four express layouts, end to end through the web endpoints."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        app = create_app(self.dir)
        app.testing = True
        self.client = app.test_client()
        self.clean = clean_vrc01()
        self.client.post("/api/parts/import", json={"path": str(PARTS)})
        self.client.post(
            "/api/antibody/save",
            json={"id": "GL_VRC01", "vh_seq": self.clean["GL_VRC01_VH"],
                  "vl_seq": self.clean["GL_VRC01_VL"]},
        )

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def express_keys(self):
        presets = self.client.get("/api/presets").get_json()["presets"]
        return sorted(k for k, v in presets.items() if v["express"])

    def test_all_four_express_layouts_plan_without_missing_parts(self):
        keys = self.express_keys()
        self.assertEqual(len(keys), 4)
        for key in keys:
            with self.subTest(preset=key):
                res = self.client.post(
                    "/api/plan", json={"ids": ["GL_VRC01"], "preset": key}
                )
                self.assertEqual(res.status_code, 200, res.get_json())
                self.assertEqual(res.get_json()["buildable"], 1)

    def test_the_isotype_reaches_the_plan(self):
        for key, label in (("single_chain_bcr_igg1", "IgG1"),
                           ("single_chain_bcr_igm", "IgM")):
            with self.subTest(preset=key):
                data = self.client.post(
                    "/api/plan", json={"ids": ["GL_VRC01"], "preset": key}
                ).get_json()
                construct = data["constructs"][0]
                self.assertEqual(construct["isotype"], label)
                self.assertEqual(construct["cl_part"], "IGKC")

    def test_an_express_construct_is_the_stored_columns_spliced(self):
        """What is built must be what express mode showed."""
        row = self.client.get("/api/library").get_json()["antibodies"][0]
        out = pathlib.Path(self.dir) / "x.fasta"
        res = self.client.post(
            "/api/build",
            json={"ids": ["GL_VRC01"], "preset": "2alinked_igm", "output": str(out)},
        )
        self.assertEqual(res.status_code, 200)
        seq = res.get_json()["constructs"][0]["preview"] or ""
        built = {c["id"]: c for c in self.client.get("/api/constructs").get_json()["constructs"]}
        whole = built[Store.construct_id("GL_VRC01", "2alinked_igm")]["seq"]
        self.assertIn(row["light_chain"], whole)
        self.assertIn(row["heavy_chain_igm"], whole)
        self.assertNotIn(row["heavy_chain_igg1"], whole)
        del seq

    def test_the_igm_and_igg1_forms_differ_only_in_the_heavy_constant(self):
        out = pathlib.Path(self.dir) / "x.fasta"
        seqs = {}
        for key in ("single_chain_bcr_igg1", "single_chain_bcr_igm"):
            self.client.post(
                "/api/build",
                json={"ids": ["GL_VRC01"], "preset": key, "output": str(out)},
            )
            rows = {c["id"]: c for c in self.client.get("/api/constructs").get_json()["constructs"]}
            seqs[key] = rows[Store.construct_id("GL_VRC01", key)]["seq"]
        igg1, igm = seqs["single_chain_bcr_igg1"], seqs["single_chain_bcr_igm"]
        self.assertNotEqual(igg1, igm)
        shared = len(assemble.INTRINSIC_SEQS["CHTMCTIgG1"])
        self.assertEqual(igg1[: -shared - 1], igm[: -len(assemble.INTRINSIC_SEQS["CHTMCTIgM"]) - 1])

    def test_a_stored_construct_records_the_light_constant(self):
        out = pathlib.Path(self.dir) / "x.fasta"
        self.client.post(
            "/api/build",
            json={"ids": ["GL_VRC01"], "preset": "single_chain_bcr_igm",
                  "output": str(out)},
        )
        row = self.client.get("/api/constructs").get_json()["constructs"][0]
        self.assertEqual(row["cl_part"], "IGKC")
        self.assertEqual(row["preset"], "single_chain_bcr_igm")

    def test_a_layout_cannot_place_a_v_domain_and_a_whole_chain_at_once(self):
        slots = [
            {"kind": "vl"}, {"kind": "light_chain"}, {"kind": "stop"},
        ]
        res = self.client.post(
            "/api/presets/save", json={"label": "Doubled", "slots": slots}
        )
        self.assertEqual(res.status_code, 400)
        self.assertIn("not both", res.get_json()["error"])

    def test_a_light_chain_cannot_be_combined_with_a_light_constant(self):
        slots = [
            {"kind": "light_chain"}, {"kind": "light_constant"}, {"kind": "stop"},
        ]
        res = self.client.post(
            "/api/presets/save", json={"label": "Doubled CL", "slots": slots}
        )
        self.assertEqual(res.status_code, 400)
        self.assertIn("already includes its constant region", res.get_json()["error"])

    def test_a_saved_layout_can_use_the_whole_chain_slots(self):
        slots = [
            {"kind": "part", "part": "SignalHC"},
            {"kind": "heavy_chain", "isotype": "igm"},
            {"kind": "stop"},
        ]
        res = self.client.post(
            "/api/presets/save", json={"label": "Heavy IgM only", "slots": slots}
        )
        self.assertEqual(res.status_code, 200, res.get_json())
        key = res.get_json()["key"]
        plan = self.client.post(
            "/api/plan", json={"ids": ["GL_VRC01"], "preset": key}
        )
        self.assertEqual(plan.status_code, 200, plan.get_json())
        self.assertEqual(plan.get_json()["constructs"][0]["isotype"], "IgM")


class ShadowedPartTests(unittest.TestCase):
    """Rows stored before their names became intrinsic."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        app = create_app(self.dir)
        app.testing = True
        self.client = app.test_client()
        self.store = Store(self.dir)
        # Write directly: the endpoints would refuse these names, which is
        # the point -- only a table that predates them can hold one.
        self.store.parts.write_all([
            self.store.part_row("IGKC", assemble.INTRINSIC_SEQS["IGKC"], "old"),
            self.store.part_row("IGLC2", "GQPKAAPSVTLFPPSSEELQ", "old"),
            self.store.part_row("SignalHC", "MEWSWVFLFFLSVTTGVHS", "old"),
        ])

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_a_shadowed_row_is_reported_and_left_out_of_the_table(self):
        data = self.client.get("/api/library").get_json()
        listed = {r["name"]: r for r in data["parts"]}
        self.assertEqual(listed["IGKC"]["builtin"], "yes")
        self.assertEqual(listed["IGKC"]["seq"], assemble.INTRINSIC_SEQS["IGKC"])
        shadowed = {r["name"]: r for r in data["shadowed"]}
        self.assertEqual(sorted(shadowed), ["IGKC", "IGLC2"])
        self.assertTrue(shadowed["IGKC"]["same"])
        self.assertFalse(shadowed["IGLC2"]["same"])

    def test_a_shadowed_row_cannot_win_at_assembly(self):
        """A hand-edited table must not be able to override an intrinsic part."""
        self.assertEqual(
            self.store.parts_as_dict()["IGLC2"], assemble.INTRINSIC_SEQS["IGLC2"]
        )

    def test_shadowed_rows_are_removed_only_when_asked(self):
        res = self.client.post(
            "/api/parts/unshadow", json={"names": ["IGKC", "IGLC2"]}
        )
        self.assertEqual(res.status_code, 200)
        self.assertEqual(sorted(res.get_json()["removed"]), ["IGKC", "IGLC2"])
        self.assertEqual(res.get_json()["shadowed"], [])
        # The user's own row is untouched, and the intrinsic parts remain.
        names = {r["name"] for r in res.get_json()["parts"]}
        self.assertIn("SignalHC", names)
        self.assertTrue(assemble.RESERVED_PART_NAMES <= names)

    def test_unshadow_refuses_a_row_that_is_not_shadowed(self):
        res = self.client.post("/api/parts/unshadow", json={"names": ["SignalHC"]})
        self.assertEqual(res.status_code, 400)
        self.assertIn("SignalHC", {r["name"] for r in self.store.parts.read()})


class IntrinsicPartTests(unittest.TestCase):
    """The shipped parts themselves."""

    def test_the_four_names_are_reserved(self):
        self.assertEqual(
            sorted(assemble.RESERVED_PART_NAMES),
            ["CHTMCTIgG1", "CHTMCTIgM", "IGKC", "IGLC2"],
        )

    def test_the_light_constants_match_the_knowledge_base(self):
        """Section 5 is the source; a drift here reaches every plasmid."""
        from bcr import rules

        self.assertTrue(
            assemble.INTRINSIC_SEQS["IGKC"].startswith(rules.CL_KAPPA_START)
        )
        self.assertEqual(len(assemble.INTRINSIC_SEQS["IGKC"]), 107)
        self.assertTrue(
            assemble.INTRINSIC_SEQS["IGLC2"].startswith(rules.CL_LAMBDA_START)
        )
        self.assertEqual(len(assemble.INTRINSIC_SEQS["IGLC2"]), 106)

    def test_every_intrinsic_part_has_a_role_read_from_its_sequence(self):
        for name, value in assemble.BUILTIN_PARTS.items():
            with self.subTest(name=name):
                self.assertEqual(value["role"], assemble.infer_role(value["seq"])[0])
                self.assertNotEqual(value["role"], "other")
                self.assertTrue(value["description"])
                self.assertNotIn("*", value["seq"])

    def test_a_damaged_intrinsic_file_is_a_broken_installation(self):
        """Carrying on with none of them would look like a pile of typos."""
        import tempfile as tf

        with tf.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            fh.write('{"parts": {}}')
            path = fh.name
        with self.assertRaises(RuntimeError):
            assemble.load_builtin_parts(path)

        with tf.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            fh.write('{"parts": {"IGKC": {"seq": "RTVAAPSVFIFPPS"}}}')
            path = fh.name
        with self.assertRaises(RuntimeError) as caught:
            assemble.load_builtin_parts(path)
        self.assertIn("missing the intrinsic part", str(caught.exception))


class UploadTests(unittest.TestCase):
    """Dropping a FASTA puts it on disk and hands back its path.

    An upload is deliberately not a third way to submit sequences: it only
    produces a path, which then goes through the ordinary two-phase flows.
    """

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        app = create_app(self.dir)
        app.testing = True
        self.client = app.test_client()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def send(self, filename, data):
        if isinstance(data, str):
            data = data.encode()
        return self.client.post(
            "/api/upload",
            data={"file": (io.BytesIO(data), filename)},
            content_type="multipart/form-data",
        )

    def test_a_dropped_fasta_lands_in_the_uploads_directory(self):
        res = self.send("vrc01.fasta", VRC01.read_bytes())
        self.assertEqual(res.status_code, 200)
        body = res.get_json()
        saved = pathlib.Path(body["path"])
        self.assertEqual(saved.parent, pathlib.Path(self.dir) / "uploads")
        self.assertEqual(saved.read_bytes(), VRC01.read_bytes())
        self.assertEqual(body["bytes"], len(VRC01.read_bytes()))

    def test_the_returned_path_drives_the_ordinary_check_flow(self):
        path = self.send("vrc01.fa", VRC01.read_bytes()).get_json()["path"]
        res = self.client.post("/api/check", json={"path": path})
        self.assertEqual(res.status_code, 200)
        self.assertEqual(res.get_json()["summary"]["passed"], 2)

    def test_two_uploads_of_one_name_do_not_overwrite_each_other(self):
        first = self.send("parts.fasta", ">SignalHC\nMGWSCIILFLVATATGVHS\n")
        second = self.send("parts.fasta", ">SignalHC\nMGWSCIILFLVATATGVHX\n")
        a, b = first.get_json()["path"], second.get_json()["path"]
        self.assertNotEqual(a, b)
        self.assertTrue(pathlib.Path(a).exists())
        self.assertTrue(pathlib.Path(b).exists())

    def test_only_fa_and_fasta_are_accepted(self):
        for name in ("notes.txt", "sequences", "archive.fasta.gz", "x.fa.exe"):
            with self.subTest(name=name):
                res = self.send(name, ">A\nQVQ\n")
                self.assertEqual(res.status_code, 400)
        self.assertFalse(list((pathlib.Path(self.dir) / "uploads").glob("*"))
                         if (pathlib.Path(self.dir) / "uploads").exists() else [])

    def test_a_filename_cannot_escape_the_uploads_directory(self):
        res = self.send("../../../etc/evil.fasta", ">A\nQVQLVQ\n")
        self.assertEqual(res.status_code, 200)
        saved = pathlib.Path(res.get_json()["path"])
        self.assertEqual(saved.parent, pathlib.Path(self.dir) / "uploads")
        self.assertNotIn("..", saved.name)
        self.assertEqual(res.get_json()["name"], "evil.fasta")

    def test_something_that_is_not_a_fasta_is_refused(self):
        res = self.send("mystery.fasta", "QVQLVQSGAEVKKPG\n")
        self.assertEqual(res.status_code, 400)
        self.assertIn(">", res.get_json()["error"])

    def test_an_empty_file_is_refused(self):
        res = self.send("blank.fasta", "   \n")
        self.assertEqual(res.status_code, 400)
        self.assertIn("empty", res.get_json()["error"].lower())

    def test_a_missing_file_field_is_refused(self):
        res = self.client.post("/api/upload", data={}, content_type="multipart/form-data")
        self.assertEqual(res.status_code, 400)

    def test_the_drop_zones_and_naming_example_are_on_the_page(self):
        page = self.client.get("/").data
        for control in (b'id="seqDrop"', b'id="seqFile"', b'id="partsDrop"',
                        b'id="partsFile"', b'accept=".fa,.fasta"'):
            self.assertIn(control, page)
        # The naming-rule cartoon, with a text alternative rather than an
        # image nobody using a screen reader can read.
        self.assertIn(b"CR9114", page)
        self.assertIn(b'<figure class="naming">', page)
        self.assertIn(b"aria-label=", page[page.index(b'<figure class="naming">'):])


if __name__ == "__main__":
    unittest.main()


class AddAnywayTests(unittest.TestCase):
    """A sanity-check failure is a judgement, and it can be overruled.

    Overruling only ever stores the sequence exactly as it came in: nothing is
    derived from it, and the row says as much. The two failures that are not
    judgements -- a sequence that is not amino acids, an empty one -- cannot be
    overruled at all.
    """

    ANYWAY = pathlib.Path(__file__).parent / "fixtures" / "add_anyway.fasta"

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        app = create_app(self.dir)
        app.testing = True
        self.client = app.test_client()
        self.clean = clean_vrc01()
        self.records = {
            r.name: r.seq for r in read_fasta(self.ANYWAY)
        }

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_overruling_re_renders_without_losing_the_reading_position(self):
        """Pressing "add anyway" must leave the panel where it stands.

        The handler re-runs the check rather than patching the view, so the
        whole block is replaced. Emptying it first collapses the panel body and
        its scroll offset clamps to zero, which throws the user back to the top
        of a long file. So the re-run is quiet -- no "Processing…" placeholder
        -- and the offset is captured and restored around it.
        """
        page = self.client.get("/").data.decode()
        handler = page[page.index('querySelectorAll("[data-anyway]")'):]
        handler = handler[:handler.index("[data-unanyway]")]
        self.assertIn('holdScroll($("checkOut"))', page)
        self.assertIn("runCheck({ quiet: true })", page)
        self.assertNotIn("Processing", handler)

        run = page[page.index("async function runCheck"):]
        run = run[:run.index("\n}")]
        self.assertIn("if (!opts.quiet)", run)
        self.assertIn("Processing", run)

        self.assertIn("function holdScroll(el)", page)

    def check(self, **extra):
        return self.client.post(
            "/api/check", json={"path": str(self.ANYWAY), **extra}
        ).get_json()

    def rows(self):
        return self.client.get("/api/library").get_json()["antibodies"]

    # --- the offer --------------------------------------------------------

    def test_a_failed_record_is_offered_as_given(self):
        body = self.check()
        failed = next(r for r in body["records"] if r["status"] == "fail")
        self.assertEqual(failed["name"], "mabX_VH")
        self.assertTrue(failed["add_anyway"])
        self.assertEqual(failed["add_anyway_chain"], "heavy")
        self.assertEqual(failed["add_anyway_blocked"], "")
        # and the pair is not addable while the failure stands
        (proposal,) = body["proposals"]
        self.assertFalse(proposal["ready"])

    def test_a_passing_record_is_not_offered_anything(self):
        passed = next(r for r in self.check()["records"] if r["status"] == "ok")
        self.assertFalse(passed["add_anyway"])
        self.assertFalse(passed["accepted_as_given"])

    def test_non_amino_acid_characters_cannot_be_overruled(self):
        path = pathlib.Path(self.dir) / "junk.fasta"
        path.write_text(">junk_VH\nQVQLZZ123X@@\n", encoding="utf-8")
        body = self.client.post("/api/check", json={"path": str(path)}).get_json()
        (record,) = body["records"]
        self.assertFalse(record["add_anyway"])
        self.assertIn("not amino acids", record["add_anyway_blocked"])

        res = self.client.post(
            "/api/check", json={"path": str(path), "overrides": ["junk_VH"]}
        )
        self.assertEqual(res.status_code, 400)
        self.assertIn("cannot be stored as given", res.get_json()["error"])

    def test_an_unlabelled_untypable_record_cannot_be_overruled(self):
        """With no _VH/_VL suffix and no detected class, nothing says which
        chain it is, so there is no slot to store it in."""
        path = pathlib.Path(self.dir) / "nameless.fasta"
        path.write_text(">mystery\nMVSKGEELFTGVVPILVELDGDVNGHKFSVSG\n",
                        encoding="utf-8")
        (record,) = self.client.post(
            "/api/check", json={"path": str(path)}
        ).get_json()["records"]
        self.assertFalse(record["add_anyway"])
        self.assertIn("_VH/_VL", record["add_anyway_blocked"])

    # --- overruling it ----------------------------------------------------

    def test_overruling_makes_the_pair_addable(self):
        body = self.check(overrides=["mabX_VH"])
        record = next(r for r in body["records"] if r["name"] == "mabX_VH")
        self.assertEqual(record["status"], "ok")
        self.assertTrue(record["accepted_as_given"])
        self.assertEqual(record["confidence"], "none")
        self.assertEqual(record["length"], len(self.records["mabX_VH"]))
        self.assertEqual(body["summary"]["ready"], 1)

        (proposal,) = body["proposals"]
        self.assertTrue(proposal["ready"])
        self.assertEqual(proposal["accepted_as_given"], ["VH"])
        # overruling the check is not a trim decision: nothing was trimmed
        self.assertEqual(proposal["decisions"], [])

    def test_check_still_writes_nothing(self):
        self.check(overrides=["mabX_VH"])
        self.assertEqual(self.rows(), [])

    def test_approve_stores_the_sequence_exactly_as_given(self):
        self.check(overrides=["mabX_VH"])
        body = self.client.post("/api/approve", json={
            "path": str(self.ANYWAY),
            "ids": ["mabX"],
            "overrides": ["mabX_VH"],
        }).get_json()
        (saved,) = body["saved"]
        self.assertEqual(saved["as_given"], ["VH"])

        (row,) = self.rows()
        self.assertEqual(row["vh_seq"], self.records["mabX_VH"])
        self.assertEqual(row["vh_confidence"], "none")
        self.assertIn("failed the sanity check", row["vh_trim"])
        # nothing was derived, so the derived columns say nothing
        self.assertEqual(row["vh_fr4"], "")
        self.assertIn("stored as given", row["notes"])
        # the half that typed normally is untouched by any of this
        self.assertEqual(row["vl_seq"], self.records["mabX_VL"])
        self.assertEqual(row["vl_confidence"], "high")

    def test_approve_without_the_override_refuses_the_pair(self):
        """The override travels with the request. Approving without it leaves
        the file typing exactly as the rules say, so there is no pair."""
        body = self.client.post("/api/approve", json={
            "path": str(self.ANYWAY), "ids": ["mabX"],
        }).get_json()
        self.assertEqual(body["saved"], [])
        (skipped,) = body["skipped"]
        self.assertIn("failed the sanity check", skipped["detail"])
        self.assertEqual(self.rows(), [])

    def test_an_override_naming_a_record_not_in_the_file_is_an_error(self):
        res = self.client.post(
            "/api/check", json={"path": str(self.ANYWAY), "overrides": ["nope_VH"]}
        )
        self.assertEqual(res.status_code, 400)
        self.assertIn("no such record", res.get_json()["error"])

    def test_overrides_must_be_a_list(self):
        for url in ("/api/check", "/api/approve"):
            res = self.client.post(url, json={
                "path": str(self.ANYWAY), "ids": ["mabX"], "overrides": "mabX_VH",
            })
            self.assertEqual(res.status_code, 400)
            self.assertIn("must be a list", res.get_json()["error"])

    def test_an_overruled_record_still_passes_typing_for_the_record_view(self):
        """Overruling does not rewrite history: the errors that stood are
        still reported, as warnings, on the record and in the stored row."""
        record = next(
            r for r in self.check(overrides=["mabX_VH"])["records"]
            if r["name"] == "mabX_VH"
        )
        self.assertEqual(record["errors"], [])
        self.assertTrue(any("framework 4" in w for w in record["warnings"]))
        self.assertTrue(any("verbatim" in w for w in record["warnings"]))

    # --- the hand editor --------------------------------------------------

    def test_a_pasted_chain_that_fails_is_offered_as_given(self):
        res = self.client.post("/api/antibody/save", json={
            "id": "Pasted",
            "vh_seq": self.records["mabX_VH"],
            "vl_seq": self.clean["GL_VRC01_VL"],
        })
        self.assertEqual(res.status_code, 400)
        body = res.get_json()
        (offer,) = body["add_anyway"]
        self.assertEqual(offer["chain"], "VH")
        self.assertEqual(offer["field"], "vh_as_given")
        self.assertTrue(offer["can"])
        self.assertEqual(offer["as_given"], len(self.records["mabX_VH"]))
        self.assertEqual(self.rows(), [])

    def test_saying_so_stores_the_pasted_chain_as_given(self):
        res = self.client.post("/api/antibody/save", json={
            "id": "Pasted",
            "vh_seq": self.records["mabX_VH"],
            "vl_seq": self.clean["GL_VRC01_VL"],
            "vh_as_given": True,
        })
        self.assertEqual(res.status_code, 200)
        self.assertTrue(any("as given" in n for n in res.get_json()["notes"]))

        (row,) = self.rows()
        self.assertEqual(row["vh_seq"], self.records["mabX_VH"])
        self.assertEqual(row["vh_confidence"], "none")

    def test_the_flag_is_per_chain(self):
        """Insisting on the heavy chain says nothing about the light one."""
        res = self.client.post("/api/antibody/save", json={
            "id": "Pasted",
            "vh_seq": self.records["mabX_VH"],
            "vl_seq": self.records["mabX_VH"],   # not a light chain either
            "vh_as_given": True,
        })
        self.assertEqual(res.status_code, 400)
        (offer,) = res.get_json()["add_anyway"]
        self.assertEqual(offer["chain"], "VL")
        self.assertEqual(self.rows(), [])

    def test_the_summary_does_not_count_an_overruled_record_as_passing(self):
        body = self.check(overrides=["mabX_VH"])
        self.assertEqual(body["summary"]["passed"], 1)
        self.assertEqual(body["summary"]["failed"], 0)
        self.assertEqual(body["summary"]["as_given"], 1)

    # --- a light chain needs kappa or lambda stated ------------------------

    def light_pair(self):
        """A clean VH beside a light record the rules reject. The light half
        is the one that cannot be stored without a kappa/lambda call."""
        path = pathlib.Path(self.dir) / "light.fasta"
        path.write_text(
            f">mabY_VH\n{self.clean['GL_VRC01_VH']}\n"
            f">mabY_VL\n{self.records['mabX_VH']}\n",
            encoding="utf-8",
        )
        return path

    def test_a_failed_light_chain_is_offered_both_subclasses(self):
        path = self.light_pair()
        body = self.client.post("/api/check", json={"path": str(path)}).get_json()
        failed = next(r for r in body["records"] if r["status"] == "fail")
        self.assertEqual(failed["name"], "mabY_VL")
        self.assertTrue(failed["add_anyway"])
        self.assertEqual(failed["add_anyway_chain"], "light")
        self.assertEqual(failed["add_anyway_subclasses"], ["kappa", "lambda"])

    def test_a_failed_heavy_chain_is_offered_no_subclass(self):
        failed = next(r for r in self.check()["records"] if r["status"] == "fail")
        self.assertEqual(failed["add_anyway_chain"], "heavy")
        self.assertEqual(failed["add_anyway_subclasses"], [])

    def test_overruling_a_light_chain_without_saying_which_is_refused(self):
        path = self.light_pair()
        for url, extra in (
            ("/api/check", {}),
            ("/api/approve", {"ids": ["mabY"]}),
        ):
            with self.subTest(url=url):
                res = self.client.post(url, json={
                    "path": str(path), "overrides": ["mabY_VL"], **extra,
                })
                self.assertEqual(res.status_code, 400)
                self.assertIn("kappa or lambda", res.get_json()["error"])
        self.assertEqual(self.rows(), [])

    def test_a_stated_subclass_chooses_the_constant_region(self):
        path = self.light_pair()
        body = self.client.post("/api/check", json={
            "path": str(path),
            "overrides": ["mabY_VL"],
            "subclasses": {"mabY_VL": "lambda"},
        }).get_json()
        record = next(r for r in body["records"] if r["name"] == "mabY_VL")
        self.assertEqual(record["subclass"], "lambda")
        (proposal,) = body["proposals"]
        self.assertTrue(proposal["ready"])
        self.assertEqual(proposal["subclass"], "lambda")

        self.client.post("/api/approve", json={
            "path": str(path),
            "ids": ["mabY"],
            "overrides": ["mabY_VL"],
            "subclasses": {"mabY_VL": "lambda"},
        })
        (row,) = self.rows()
        self.assertEqual(row["light_subclass"], "lambda")
        # the whole light chain is derivable again, with the CL the user named
        self.assertEqual(
            row["light_chain"],
            self.records["mabX_VH"] + assemble.INTRINSIC_SEQS["IGLC2"],
        )
        self.assertTrue(any("stated by you" in n for n in row["notes"].split("; ")),
                        row["notes"])

    def test_kappa_takes_igkc(self):
        path = self.light_pair()
        self.client.post("/api/approve", json={
            "path": str(path),
            "ids": ["mabY"],
            "overrides": ["mabY_VL"],
            "subclasses": {"mabY_VL": "kappa"},
        })
        (row,) = self.rows()
        self.assertEqual(row["light_subclass"], "kappa")
        self.assertTrue(
            row["light_chain"].endswith(assemble.INTRINSIC_SEQS["IGKC"])
        )

    def test_a_subclass_that_is_neither_is_refused(self):
        path = self.light_pair()
        res = self.client.post("/api/check", json={
            "path": str(path),
            "overrides": ["mabY_VL"],
            "subclasses": {"mabY_VL": "IGKC"},
        })
        self.assertEqual(res.status_code, 400)
        self.assertIn("kappa or lambda, not 'igkc'", res.get_json()["error"])

    def test_subclasses_must_be_an_object(self):
        res = self.client.post("/api/check", json={
            "path": str(self.ANYWAY), "subclasses": ["kappa"],
        })
        self.assertEqual(res.status_code, 400)
        self.assertIn("must be an object", res.get_json()["error"])

    def test_the_hand_editor_asks_which_subclass_before_storing(self):
        res = self.client.post("/api/antibody/save", json={
            "id": "Pasted",
            "vh_seq": self.clean["GL_VRC01_VH"],
            "vl_seq": self.records["mabX_VH"],
            "vl_as_given": True,          # yes, but not which
        })
        self.assertEqual(res.status_code, 400)
        (offer,) = res.get_json()["add_anyway"]
        self.assertEqual(offer["chain"], "VL")
        self.assertTrue(offer["can"])
        self.assertEqual(offer["subclasses"], ["kappa", "lambda"])
        self.assertEqual(self.rows(), [])

    def test_the_hand_editor_takes_the_subclass_in_the_flag(self):
        res = self.client.post("/api/antibody/save", json={
            "id": "Pasted",
            "vh_seq": self.clean["GL_VRC01_VH"],
            "vl_seq": self.records["mabX_VH"],
            "vl_as_given": "kappa",
        })
        self.assertEqual(res.status_code, 200)
        (row,) = self.rows()
        self.assertEqual(row["light_subclass"], "kappa")
        self.assertEqual(row["vl_confidence"], "none")


class MultiLayoutBuildTests(unittest.TestCase):
    """One run, several layouts.

    The build mode is part of every construct name, so a run that builds the
    same antibody four ways writes four records to one file and keeps four
    rows, none of them shadowing another.
    """

    EXPRESS = ("single_chain_bcr_igg1", "single_chain_bcr_igm",
               "2alinked_igg1", "2alinked_igm")

    def setUp(self):
        self.dir = pathlib.Path(tempfile.mkdtemp())
        app = create_app(self.dir)
        app.testing = True
        self.client = app.test_client()
        self.out = self.dir / "built.fasta"
        self.client.post("/api/parts/apply", json={"path": str(PARTS)})
        self.client.post("/api/check", json={"path": str(VRC01)})
        self.client.post("/api/approve",
                         json={"path": str(VRC01), "ids": ["GL_VRC01"]})

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def plan(self, **fields):
        return self.client.post("/api/plan", json={"ids": ["GL_VRC01"], **fields})

    def build(self, **fields):
        return self.client.post(
            "/api/build",
            json={"ids": ["GL_VRC01"], "output": str(self.out), **fields},
        )

    def test_a_plan_covers_every_layout_asked_for(self):
        body = self.plan(presets=list(self.EXPRESS)).get_json()
        self.assertEqual(body["buildable"], 4)
        self.assertEqual(len(body["constructs"]), 4)
        self.assertEqual([p["key"] for p in body["presets"]], list(self.EXPRESS))
        self.assertEqual([p["buildable"] for p in body["presets"]], [1, 1, 1, 1])
        # every planned construct says which layout it is and what it will
        # be called
        for entry in body["constructs"]:
            self.assertIn(entry["preset"], self.EXPRESS)
            self.assertEqual(
                entry["construct_id"],
                Store.construct_id("GL_VRC01", entry["preset"]),
            )

    def test_the_matrix_is_antibodies_times_layouts(self):
        clean = clean_vrc01()
        self.client.post("/api/antibody/save", json={
            "id": "Second", "vh_seq": clean["GL_VRC01_VH"],
            "vl_seq": clean["GL_VRC01_VL"],
        })
        body = self.client.post("/api/plan", json={
            "ids": ["GL_VRC01", "Second"], "presets": list(self.EXPRESS[:2]),
        }).get_json()
        self.assertEqual(len(body["constructs"]), 4)
        self.assertEqual(body["buildable"], 4)

    def test_one_run_writes_one_record_per_layout(self):
        body = self.build(presets=list(self.EXPRESS)).get_json()
        self.assertEqual(body["written"], 4)
        self.assertEqual([p["written"] for p in body["by_preset"]], [1, 1, 1, 1])

        written = read_fasta(self.out)
        self.assertEqual(len(written), 4)
        names = [r.name for r in written]
        self.assertEqual(
            names, [Store.construct_id("GL_VRC01", k) for k in self.EXPRESS]
        )
        # the whole point: unique names in one file
        self.assertEqual(len(set(names)), 4)
        # and four distinct sequences, each ending in exactly one stop
        self.assertEqual(len({r.seq for r in written}), 4)
        for record in written:
            self.assertTrue(record.seq.endswith("*"))
            self.assertEqual(record.seq.count("*"), 1)

    def test_the_build_mode_is_in_the_fasta_description_too(self):
        self.build(presets=["single_chain_bcr_igm"])
        (record,) = read_fasta(self.out)
        self.assertIn("Single chain BCR_IGM", record.description)

    def test_each_layout_keeps_its_own_row(self):
        self.build(presets=list(self.EXPRESS))
        rows = {
            c["id"]: c
            for c in self.client.get("/api/constructs").get_json()["constructs"]
        }
        self.assertEqual(len(rows), 4)
        for key in self.EXPRESS:
            row = rows[Store.construct_id("GL_VRC01", key)]
            self.assertEqual(row["preset"], key)
            self.assertEqual(row["antibody_id"], "GL_VRC01")

    def test_rebuilding_the_same_layouts_updates_rather_than_duplicates(self):
        self.build(presets=list(self.EXPRESS))
        body = self.build(presets=list(self.EXPRESS)).get_json()
        self.assertEqual(len(body["library"]), 4)
        self.assertEqual(body["replaced"], 4)
        self.assertEqual(len(read_fasta(self.out)), 4)

    def test_a_single_layout_run_is_unchanged(self):
        """The one-layout case still answers in the flat shape it always did."""
        body = self.plan(preset="single_chain_bcr_igm").get_json()
        self.assertEqual(body["preset"], "single_chain_bcr_igm")
        self.assertEqual(body["preset_label"], "Single chain BCR_IGM")
        self.assertEqual(len(body["constructs"]), 1)

    def test_a_run_of_several_names_no_single_preset(self):
        body = self.plan(presets=list(self.EXPRESS[:2])).get_json()
        self.assertEqual(body["preset"], "")
        self.assertEqual(body["preset_label"], "")

    def test_the_same_layout_twice_is_built_once(self):
        body = self.build(presets=["2alinked_igg1", "2alinked_igg1"]).get_json()
        self.assertEqual(body["written"], 1)

    def test_an_unknown_layout_in_the_list_stops_the_run(self):
        for url in ("/api/plan", "/api/build"):
            with self.subTest(url=url):
                res = self.client.post(url, json={
                    "ids": ["GL_VRC01"],
                    "presets": ["single_chain_bcr_igg1", "nonesuch"],
                    "output": str(self.out),
                })
                self.assertEqual(res.status_code, 400)
                self.assertIn("nonesuch", res.get_json()["error"])
        self.assertFalse(self.out.exists())

    def test_an_empty_layout_list_is_refused(self):
        res = self.plan(presets=[])
        self.assertEqual(res.status_code, 400)
        self.assertIn("No layout selected", res.get_json()["error"])

    def test_presets_must_be_a_list(self):
        res = self.plan(presets="single_chain_bcr_igg1")
        self.assertEqual(res.status_code, 400)
        self.assertIn("must be a list", res.get_json()["error"])

    def test_a_layout_missing_parts_stops_the_whole_run(self):
        """Not a silent omission: the user picked that layout on purpose."""
        slots = [{"kind": "part", "part": "Nowhere"},
                 {"kind": "light_chain"}, {"kind": "stop"}]
        saved = self.client.post("/api/presets/save", json={
            "label": "Needs a missing part", "slots": slots,
        })
        key = saved.get_json()["key"]
        res = self.plan(presets=["single_chain_bcr_igg1", key])
        self.assertEqual(res.status_code, 400)
        body = res.get_json()
        self.assertIn("Nowhere", body["error"])
        self.assertIn(key, body["missing_by_preset"])
        self.assertNotIn("single_chain_bcr_igg1", body["missing_by_preset"])

    def test_skipped_entries_say_which_layout_they_were_for(self):
        body = self.client.post("/api/build", json={
            "ids": ["GL_VRC01", "nope"],
            "presets": ["single_chain_bcr_igg1", "single_chain_bcr_igm"],
            "output": str(self.out),
        }).get_json()
        self.assertEqual(body["written"], 2)
        self.assertEqual(
            sorted((s["id"], s["preset"]) for s in body["skipped"]),
            [("nope", "single_chain_bcr_igg1"), ("nope", "single_chain_bcr_igm")],
        )

    def test_exporting_several_rows_keeps_their_names_apart(self):
        self.build(presets=list(self.EXPRESS))
        target = self.dir / "export.fasta"
        ids = [Store.construct_id("GL_VRC01", k) for k in self.EXPRESS]
        res = self.client.post("/api/constructs/export",
                               json={"ids": ids, "output": str(target)})
        self.assertEqual(res.status_code, 200)
        exported = read_fasta(target)
        self.assertEqual([r.name for r in exported], ids)
        self.assertEqual(len({r.name for r in exported}), 4)


class LegacyConstructNameTests(unittest.TestCase):
    """Rows keyed before the layout was part of the name.

    `<antibody>_BCR` cannot tell two layouts of one antibody apart. The row
    records its layout, so the name it should have is computable -- but a
    rename is a modification, so it is reported and done only when asked.
    """

    def setUp(self):
        self.dir = pathlib.Path(tempfile.mkdtemp())
        app = create_app(self.dir)
        app.testing = True
        self.client = app.test_client()
        self.store = Store(self.dir)
        self.client.post("/api/parts/apply", json={"path": str(PARTS)})
        self.client.post("/api/check", json={"path": str(VRC01)})
        self.client.post("/api/approve",
                         json={"path": str(VRC01), "ids": ["GL_VRC01"]})
        self.client.post("/api/build", json={
            "ids": ["GL_VRC01"], "preset": "single_chain_bcr_igm",
            "output": str(self.dir / "out.fasta"),
        })
        # Rewrite the one row back to the old key, as a table from before
        # this change would hold it.
        rows = self.store.constructs.read()
        self.new_id = rows[0]["id"]
        rows[0]["id"] = "GL_VRC01_BCR"
        self.store.constructs.write_all(rows)

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def legacy(self):
        return self.client.get("/api/constructs").get_json()["legacy"]

    def test_a_legacy_row_is_reported_with_the_name_it_should_have(self):
        (entry,) = self.legacy()
        self.assertEqual(entry["id"], "GL_VRC01_BCR")
        self.assertEqual(entry["new_id"], self.new_id)
        self.assertEqual(entry["preset"], "single_chain_bcr_igm")

    def test_reporting_changes_nothing(self):
        self.legacy()
        self.assertIn("GL_VRC01_BCR", self.store.constructs.read_map())

    def test_renaming_moves_the_row_and_leaves_the_build_alone(self):
        before = self.store.constructs.read_map()["GL_VRC01_BCR"]
        body = self.client.post("/api/constructs/rename-legacy",
                                json={"ids": []}).get_json()
        self.assertEqual([m["new_id"] for m in body["moved"]], [self.new_id])
        self.assertEqual(body["legacy"], [])

        rows = self.store.constructs.read_map()
        self.assertNotIn("GL_VRC01_BCR", rows)
        row = rows[self.new_id]
        for field in ("seq", "layout", "segments_json", "preset", "antibody_id"):
            self.assertEqual(row[field], before[field])

    def test_a_rename_that_would_collide_is_skipped(self):
        """The built row wins; the legacy one is left for the user to look at."""
        self.client.post("/api/build", json={
            "ids": ["GL_VRC01"], "preset": "single_chain_bcr_igm",
            "output": str(self.dir / "out2.fasta"),
        })
        self.assertEqual(len(self.store.constructs.read()), 2)
        body = self.client.post("/api/constructs/rename-legacy",
                                json={"ids": []}).get_json()
        self.assertEqual(body["moved"], [])
        self.assertIn("GL_VRC01_BCR", self.store.constructs.read_map())

    def test_a_name_the_user_chose_is_not_touched(self):
        self.client.post("/api/construct/save", json={
            "old_id": "GL_VRC01_BCR", "id": "VRC01_v1", "user_notes": "",
        })
        self.assertEqual(self.legacy(), [])

    def test_a_freshly_built_row_is_never_legacy(self):
        self.client.post("/api/build", json={
            "ids": ["GL_VRC01"], "presets": ["2alinked_igg1"],
            "output": str(self.dir / "out3.fasta"),
        })
        self.assertEqual(
            [e["id"] for e in self.legacy()], ["GL_VRC01_BCR"]
        )

    def test_renaming_with_nothing_to_rename_is_refused(self):
        self.client.post("/api/constructs/rename-legacy", json={"ids": []})
        res = self.client.post("/api/constructs/rename-legacy", json={"ids": []})
        self.assertEqual(res.status_code, 400)

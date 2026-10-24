"""Tests for layout validation, preset storage, and the preset builder API.

    python3 -m unittest discover -s tests -t . -v
"""

import json
import pathlib
import shutil
import tempfile
import unittest

from bcr import assemble, read_fasta
from bcr.presets import PresetStore, slugify
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

HEAVY_ONLY = [
    {"kind": "part", "part": "SignalHC"},
    {"kind": "vh"},
    {"kind": "part", "part": "CHandCTSC"},
    {"kind": "stop"},
]


class SlugTests(unittest.TestCase):
    def test_label_becomes_a_key(self):
        self.assertEqual(slugify("Heavy chain only"), "heavy_chain_only")
        self.assertEqual(slugify("scFv-Fc"), "scfv_fc")
        self.assertEqual(slugify("  Single  chain  BCR "), "single_chain_bcr")

    def test_unusable_label_still_yields_a_key(self):
        self.assertEqual(slugify("!!!"), "preset")
        self.assertEqual(slugify(""), "preset")


class BuiltinPresetFileTests(unittest.TestCase):
    """The shipped layouts are data in builtin_presets.json, not code."""

    def test_the_file_exists_and_is_what_is_loaded(self):
        self.assertTrue(assemble.BUILTIN_PRESET_FILE.exists())
        presets, default = assemble.load_builtin_presets()
        self.assertEqual(presets, assemble.PRESETS)
        self.assertEqual(default, assemble.DEFAULT_PRESET)

    def test_the_default_layout_comes_from_the_file(self):
        data = json.loads(assemble.BUILTIN_PRESET_FILE.read_text(encoding="utf-8"))
        key = data["default"]
        self.assertEqual(key, assemble.DEFAULT_PRESET)
        self.assertEqual(
            data["presets"][key]["slots"], assemble.PRESETS[key]["slots"]
        )
        self.assertEqual(
            data["presets"][key]["label"], assemble.PRESETS[key]["label"]
        )

    def test_no_preset_label_is_hardcoded_in_python(self):
        """Editing the JSON must be enough to rename or relayout a preset."""
        for module in ("assemble.py", "presets.py", "web.py", "cli.py", "store.py"):
            source = (pathlib.Path(assemble.__file__).parent / module).read_text(
                encoding="utf-8"
            )
            with self.subTest(module=module):
                self.assertNotIn("Single chain BCR", source)
                self.assertNotIn("single_chain_bcr", source)

    def _load(self, payload):
        path = pathlib.Path(self.tmp) / "p.json"
        path.write_text(
            payload if isinstance(payload, str) else json.dumps(payload),
            encoding="utf-8",
        )
        return assemble.load_builtin_presets(path)

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_a_valid_custom_file_loads(self):
        presets, default = self._load(
            {"default": "x", "presets": {"x": {"slots": [{"kind": "vh"}]}}}
        )
        self.assertEqual(default, "x")
        # label and description default to something usable.
        self.assertEqual(presets["x"]["label"], "x")
        self.assertEqual(presets["x"]["description"], "")

    def test_a_broken_file_raises_rather_than_silently_emptying(self):
        cases = [
            ("{ not json", "not valid JSON"),
            ({"presets": {}}, "defines no presets"),
            ({"default": "x", "presets": {"y": {"slots": []}}}, "is not one of"),
            ({"default": "x", "presets": {"x": {"label": "x"}}}, "no slot list"),
        ]
        for payload, phrase in cases:
            with self.subTest(payload=payload):
                with self.assertRaises(RuntimeError) as caught:
                    self._load(payload)
                self.assertIn(phrase, str(caught.exception))

    def test_a_missing_file_raises(self):
        with self.assertRaises(RuntimeError) as caught:
            assemble.load_builtin_presets(pathlib.Path(self.tmp) / "nope.json")
        self.assertIn("missing", str(caught.exception))

    def test_the_shipped_layouts_are_structurally_valid(self):
        for key, value in assemble.PRESETS.items():
            with self.subTest(preset=key):
                self.assertEqual(assemble.validate_slots(value["slots"]), [])


class ValidateSlotsTests(unittest.TestCase):
    def test_the_builtin_layout_validates(self):
        parts, _ = assemble.load_parts(PARTS)
        slots = assemble.PRESETS[assemble.DEFAULT_PRESET]["slots"]
        self.assertEqual(assemble.validate_slots(slots), [])
        self.assertEqual(assemble.missing_parts(slots, parts), [])

    def test_an_empty_layout_is_refused(self):
        self.assertTrue(assemble.validate_slots([]))
        self.assertTrue(assemble.validate_slots("nonsense"))

    def test_unknown_kind_is_refused(self):
        problems = assemble.validate_slots([{"kind": "nope"}])
        self.assertIn("unknown kind", problems[0])

    def test_part_slot_needs_a_name(self):
        self.assertIn("no part named", assemble.validate_slots([{"kind": "part"}])[0])

    def test_structure_check_ignores_whether_a_part_is_loaded(self):
        """Availability is `missing_parts`'s job, not the structure check's.

        A layout may be designed before its parts file exists, so naming a
        part that is not loaded is not a structural fault.
        """
        slots = [{"kind": "part", "part": "Nope"}]
        self.assertEqual(assemble.validate_slots(slots), [])
        self.assertEqual(assemble.missing_parts(slots, {"SignalHC": "M"}), ["Nope"])

    def test_stop_must_be_last_and_unique(self):
        self.assertIn(
            "must be the last component",
            assemble.validate_slots([{"kind": "stop"}, {"kind": "vh"}])[0],
        )
        self.assertIn(
            "at most one stop codon",
            assemble.validate_slots([{"kind": "stop"}, {"kind": "stop"}])[0],
        )

    def test_v_domains_and_light_constant_are_single(self):
        for slot, phrase in (
            ({"kind": "vh"}, "at most one VH"),
            ({"kind": "vl"}, "at most one VL"),
            ({"kind": "light_constant"}, "at most one light constant"),
        ):
            with self.subTest(kind=slot["kind"]):
                problems = assemble.validate_slots([slot, slot])
                self.assertTrue(any(phrase in p for p in problems), problems)


class PresetStoreTests(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.store = PresetStore(self.dir)

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_builtins_are_present_without_a_file(self):
        self.assertFalse(self.store.path.exists())
        presets = self.store.all()
        self.assertIn(assemble.DEFAULT_PRESET, presets)
        self.assertTrue(presets[assemble.DEFAULT_PRESET]["builtin"])

    def test_save_then_read_back(self):
        self.assertEqual(
            self.store.save("heavy_only", "Heavy only", "d", HEAVY_ONLY, "now"),
            "added",
        )
        self.assertEqual(
            self.store.save("heavy_only", "Heavy only", "d", HEAVY_ONLY, "now"),
            "updated",
        )
        entry = self.store.all()["heavy_only"]
        self.assertFalse(entry["builtin"])
        self.assertEqual(entry["slots"], HEAVY_ONLY)
        self.assertEqual(self.store.slots_for("heavy_only"), HEAVY_ONLY)

    def test_json_on_disk_is_readable(self):
        self.store.save("heavy_only", "Heavy only", "", HEAVY_ONLY, "now")
        data = json.loads(self.store.path.read_text(encoding="utf-8"))
        self.assertEqual(data["heavy_only"]["slots"], HEAVY_ONLY)

    def test_a_builtin_key_cannot_be_overwritten(self):
        with self.assertRaises(ValueError):
            self.store.save(
                assemble.DEFAULT_PRESET, "Hijack", "", HEAVY_ONLY, "now"
            )
        with self.assertRaises(ValueError):
            self.store.delete(assemble.DEFAULT_PRESET)

    def test_a_builtin_key_in_the_file_is_ignored(self):
        """A hand-edited file must not be able to shadow a built-in layout."""
        self.store.path.parent.mkdir(parents=True, exist_ok=True)
        self.store.path.write_text(
            json.dumps({assemble.DEFAULT_PRESET: {"label": "x", "slots": []}}),
            encoding="utf-8",
        )
        entry = self.store.all()[assemble.DEFAULT_PRESET]
        self.assertTrue(entry["builtin"])
        self.assertEqual(
            entry["label"], assemble.PRESETS[assemble.DEFAULT_PRESET]["label"]
        )
        self.assertNotEqual(entry["label"], "x")

    def test_delete(self):
        self.store.save("heavy_only", "Heavy only", "", HEAVY_ONLY, "now")
        self.assertTrue(self.store.delete("heavy_only"))
        self.assertFalse(self.store.delete("heavy_only"))
        self.assertNotIn("heavy_only", self.store.all())

    def test_a_damaged_file_does_not_break_the_builtins(self):
        self.store.path.parent.mkdir(parents=True, exist_ok=True)
        self.store.path.write_text("{ not json", encoding="utf-8")
        self.assertEqual(self.store.read_user(), {})
        self.assertIn(assemble.DEFAULT_PRESET, self.store.all())

    def test_slots_for_an_unknown_key_is_none(self):
        self.assertIsNone(self.store.slots_for("nonesuch"))


class PresetApiTests(unittest.TestCase):
    def setUp(self):
        self.dir = pathlib.Path(tempfile.mkdtemp())
        app = create_app(self.dir)
        app.testing = True
        self.client = app.test_client()
        self.client.post("/api/parts/apply", json={"path": str(PARTS)})
        self.client.post("/api/check", json={"path": str(VRC01)})
        self.client.post(
            "/api/approve",
            json={"path": str(VRC01), "ids": ["GL_VRC01"]},
        )

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def palette(self):
        return self.client.get("/api/presets").get_json()["palette"]

    def test_palette_has_one_tag_per_part_plus_the_domains(self):
        labels = [i["label"] for i in self.palette()]
        self.assertEqual(
            labels,
            [
                # the antibody's own domains, and the light constant as one tag
                "VL", "VH", "IGKC / IGLC2",
                # the whole chains express layouts place
                "Light chain", "Heavy chain IgG1", "Heavy chain IgM",
                # the intrinsic heavy constants come first, then the user's
                "CHTMCTIgG1", "CHTMCTIgM",
                "SignalHC", "LinkerSC", "SignalLC", "Furin", "F2A", "T2A",
                "EGFP", "CHandCTSC",
                "Stop",
            ],
        )

    def test_the_whole_chain_tags_carry_their_isotype(self):
        by_label = {i["label"]: i for i in self.palette()}
        self.assertEqual(by_label["Heavy chain IgG1"]["isotype"], "igg1")
        self.assertEqual(by_label["Heavy chain IgM"]["isotype"], "igm")
        self.assertEqual(by_label["Light chain"]["kind"], "light_chain")

    def test_intrinsic_parts_are_grouped_apart_from_the_users(self):
        groups = {i["label"]: i["group"] for i in self.palette()}
        self.assertEqual(groups["CHTMCTIgG1"], "Intrinsic")
        self.assertEqual(groups["SignalHC"], "Parts")

    def test_igkc_and_iglc2_share_a_single_tag(self):
        """Not one tag each: a layout chooses between them per antibody."""
        items = self.palette()
        self.assertNotIn("IGKC", [i["label"] for i in items])
        self.assertNotIn("IGLC2", [i["label"] for i in items])
        light = [i for i in items if i["kind"] == "light_constant"]
        self.assertEqual(len(light), 1)
        self.assertEqual(light[0]["label"], "IGKC / IGLC2")

    def test_palette_part_tags_carry_their_length(self):
        signal = next(i for i in self.palette() if i["label"] == "SignalHC")
        self.assertEqual(signal["kind"], "part")
        self.assertEqual(signal["part"], "SignalHC")
        self.assertEqual(signal["detail"], "19 aa")

    def test_describe_renders_a_draft_without_saving(self):
        d = self.client.post(
            "/api/presets/describe", json={"slots": HEAVY_ONLY}
        ).get_json()
        self.assertEqual(d["problems"], [])
        self.assertEqual(
            [s["label"] for s in d["layout"]],
            ["SignalHC", "VH", "CHandCTSC", "Stop"],
        )
        self.assertNotIn("heavy_only", self.client.get("/api/presets").get_json()["presets"])

    def test_describe_reports_problems(self):
        d = self.client.post(
            "/api/presets/describe",
            json={"slots": [{"kind": "stop"}, {"kind": "vh"}]},
        ).get_json()
        self.assertTrue(d["problems"])

    def test_save_and_use_a_custom_preset(self):
        res = self.client.post(
            "/api/presets/save",
            json={"label": "Heavy only", "description": "no light chain", "slots": HEAVY_ONLY},
        )
        self.assertEqual(res.status_code, 200)
        self.assertEqual(res.get_json()["key"], "heavy_only")

        presets = self.client.get("/api/presets").get_json()["presets"]
        self.assertIn("heavy_only", presets)
        self.assertFalse(presets["heavy_only"]["builtin"])

        plan = self.client.post(
            "/api/plan", json={"ids": ["GL_VRC01"], "preset": "heavy_only"}
        ).get_json()
        self.assertEqual(plan["buildable"], 1)
        construct = plan["constructs"][0]
        self.assertEqual(
            [s["label"] for s in construct["segments"]],
            ["SignalHC", "VH:GL_VRC01", "CHandCTSC", "Stop"],
        )
        # No light constant slot, so no kappa/lambda resolution is reported.
        self.assertEqual(construct["cl_part"], "")

    def test_building_with_a_custom_preset_writes_that_layout(self):
        self.client.post(
            "/api/presets/save", json={"label": "Heavy only", "slots": HEAVY_ONLY}
        )
        out = self.dir / "heavy.fasta"
        res = self.client.post(
            "/api/build",
            json={"ids": ["GL_VRC01"], "preset": "heavy_only", "output": str(out)},
        )
        self.assertEqual(res.status_code, 200)
        written = read_fasta(out)[0]
        parts = {r.name: r.seq.rstrip("*") for r in read_fasta(PARTS)}
        clean = clean_vrc01()
        self.assertEqual(
            written.seq,
            parts["SignalHC"] + clean["GL_VRC01_VH"] + parts["CHandCTSC"] + "*",
        )
        self.assertEqual(res.get_json()["library"][0]["preset"], "heavy_only")

    def test_save_refuses_an_unnamed_preset(self):
        res = self.client.post("/api/presets/save", json={"slots": HEAVY_ONLY})
        self.assertEqual(res.status_code, 400)
        self.assertIn("name", res.get_json()["error"])

    def test_save_refuses_a_structurally_invalid_layout(self):
        res = self.client.post(
            "/api/presets/save",
            json={"label": "Bad", "slots": [{"kind": "stop"}, {"kind": "vh"}]},
        )
        self.assertEqual(res.status_code, 400)
        self.assertTrue(res.get_json()["problems"])

    def test_save_allows_a_layout_whose_parts_are_not_loaded_yet(self):
        """Design first, import the parts later -- but it cannot be built."""
        res = self.client.post(
            "/api/presets/save",
            json={"label": "Future", "slots": [{"kind": "part", "part": "NotYet"},
                                                {"kind": "vh"}]},
        )
        self.assertEqual(res.status_code, 200)
        presets = self.client.get("/api/presets").get_json()["presets"]
        self.assertEqual(
            [m["part"] for m in presets["future"]["missing"]], ["NotYet"]
        )
        blocked = self.client.post(
            "/api/plan", json={"ids": ["GL_VRC01"], "preset": "future"}
        )
        self.assertEqual(blocked.status_code, 400)

    def test_save_refuses_a_builtin_key(self):
        res = self.client.post(
            "/api/presets/save",
            json={"label": "x", "key": assemble.DEFAULT_PRESET, "slots": HEAVY_ONLY},
        )
        self.assertEqual(res.status_code, 400)
        self.assertIn("built-in", res.get_json()["error"])

    def test_delete_a_saved_preset(self):
        self.client.post(
            "/api/presets/save", json={"label": "Heavy only", "slots": HEAVY_ONLY}
        )
        res = self.client.post("/api/presets/delete", json={"key": "heavy_only"})
        self.assertEqual(res.status_code, 200)
        self.assertNotIn(
            "heavy_only", self.client.get("/api/presets").get_json()["presets"]
        )

    def test_delete_refuses_a_builtin(self):
        res = self.client.post(
            "/api/presets/delete", json={"key": assemble.DEFAULT_PRESET}
        )
        self.assertEqual(res.status_code, 400)

    def test_plan_refuses_an_unknown_preset(self):
        res = self.client.post(
            "/api/plan", json={"ids": ["GL_VRC01"], "preset": "nonesuch"}
        )
        self.assertEqual(res.status_code, 400)

    def test_builder_card_is_served(self):
        page = self.client.get("/").data
        self.assertIn(b"Assembly presets", page)
        self.assertIn(b'id="draftChain"', page)
        self.assertIn(b'id="palette"', page)
        self.assertIn(b'draggable="true"', page)


class MissingPartAlertTests(unittest.TestCase):
    """A preset naming parts the table lacks must alert, not build silently.

    OLD_PARTS is the original parts file, whose headers are Signal / Linker /
    CHandCT. The built-in preset now names SignalHC / LinkerSC / CHandCTSC, so
    importing it is exactly the mismatch the alert exists for.
    """

    OLD_PARTS = ROOT / "misc.fasta"

    def setUp(self):
        self.dir = pathlib.Path(tempfile.mkdtemp())
        app = create_app(self.dir)
        app.testing = True
        self.client = app.test_client()
        self.client.post("/api/parts/apply", json={"path": str(self.OLD_PARTS)})
        self.client.post("/api/check", json={"path": str(VRC01)})
        self.client.post(
            "/api/approve",
            json={"path": str(VRC01), "ids": ["GL_VRC01"]},
        )

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_missing_parts_are_listed_on_the_preset(self):
        presets = self.client.get("/api/presets").get_json()["presets"]
        missing = presets[assemble.DEFAULT_PRESET]["missing"]
        # CHTMCTIgG1 is intrinsic, so it is never missing; the two parts the
        # layout names by hand are what the old file does not carry.
        self.assertEqual([m["part"] for m in missing], ["SignalHC", "LinkerSC"])

    def test_a_near_miss_name_gets_a_suggestion(self):
        presets = self.client.get("/api/presets").get_json()["presets"]
        hints = {
            m["part"]: m["suggestion"]
            for m in presets[assemble.DEFAULT_PRESET]["missing"]
        }
        self.assertEqual(hints["SignalHC"], "Signal")
        self.assertEqual(hints["LinkerSC"], "Linker")

    def test_plan_refuses_and_names_what_is_missing(self):
        res = self.client.post("/api/plan", json={"ids": ["GL_VRC01"]})
        self.assertEqual(res.status_code, 400)
        data = res.get_json()
        self.assertIn("not in the parts table", data["error"])
        self.assertIn("SignalHC", data["error"])
        self.assertIn("did you mean 'Signal'", data["error"])
        self.assertEqual(len(data["missing"]), 2)

    def test_build_refuses_rather_than_assembling_a_gap(self):
        out = self.dir / "x.fasta"
        res = self.client.post(
            "/api/build", json={"ids": ["GL_VRC01"], "output": str(out)}
        )
        self.assertEqual(res.status_code, 400)
        self.assertFalse(out.exists())

    def test_describe_flags_missing_parts_while_editing(self):
        d = self.client.post(
            "/api/presets/describe",
            json={"slots": [{"kind": "part", "part": "SignalHC"}, {"kind": "vh"}]},
        ).get_json()
        self.assertEqual([m["part"] for m in d["missing"]], ["SignalHC"])
        # The layout itself is structurally valid; the part is simply absent.
        self.assertEqual(d["problems"], [])

    def test_renaming_the_parts_resolves_the_alert(self):
        self.client.post("/api/parts/apply", json={"path": str(PARTS)})
        presets = self.client.get("/api/presets").get_json()["presets"]
        self.assertEqual(presets[assemble.DEFAULT_PRESET]["missing"], [])
        res = self.client.post("/api/plan", json={"ids": ["GL_VRC01"]})
        self.assertEqual(res.status_code, 200)
        self.assertEqual(res.get_json()["buildable"], 1)

    def test_a_light_constant_slot_requires_both_igkc_and_iglc2(self):
        """Either could be chosen, so both must be present."""
        kappa_only = self.dir / "kappa_only.fasta"
        kappa_only.write_text(
            ">IGKC\nRTVAAPSVFIFPPSDEQLK\n>IGLC2\nGQPKAAPSVTLFPPSSEELQ\n",
            encoding="utf-8",
        )
        parts = {"IGKC": "RTVAAPSVFIFPPSDEQLK"}
        self.assertEqual(
            assemble.missing_parts([{"kind": "light_constant"}], parts), ["IGLC2"]
        )


if __name__ == "__main__":
    unittest.main()

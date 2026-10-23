"""Local web UI for the BCR constructor.

    python3 -m bcr.web

Binds to 127.0.0.1 only. The app reads FASTA files by path from the machine it
runs on, so it is a local tool and must not be exposed to a network.

The sanity-check card drives bcr.chains.process -- the same function the
`clean` CLI calls -- so the web UI and the command line cannot drift apart.
"""

import argparse
import traceback
from pathlib import Path

from flask import Flask, jsonify, render_template, request

from . import assemble
from .chains import parse_declared, process
from .fasta import Record, read_fasta, write_fasta
from .presets import PresetStore, slugify
from .store import Store, utc_now


def _trim_note(result) -> str:
    bits = []
    if result.n_trimmed:
        bits.append(f"-{result.n_trimmed} N")
    if result.c_trimmed:
        bits.append(f"-{result.c_trimmed} C")
    return ", ".join(bits) if bits else "already clean"


def _resolve(raw_path: str) -> Path:
    """Turn user input into a readable file path, or raise ValueError."""
    if not raw_path or not raw_path.strip():
        raise ValueError("No file path given.")
    path = Path(raw_path.strip()).expanduser()
    if not path.is_absolute():
        path = Path.cwd() / path
    if not path.exists():
        raise ValueError(f"No such file: {path}")
    if path.is_dir():
        raise ValueError(f"That is a directory, not a FASTA file: {path}")
    if not path.is_file():
        raise ValueError(f"Not a regular file: {path}")
    return path


def create_app(data_dir=None) -> Flask:
    app = Flask(__name__)
    store = Store(data_dir)
    preset_store = PresetStore(store.data_dir)

    @app.get("/")
    def index():
        return render_template(
            "index.html",
            default_sequences=str(Path.cwd() / "sequences.fasta"),
            default_parts=str(Path.cwd() / "misc.fasta"),
            default_output=str(Path.cwd() / "bcr_constructs.fasta"),
            data_dir=str(store.data_dir),
        )

    @app.get("/api/library")
    def library():
        return jsonify(
            antibodies=store.antibodies.read(),
            parts=store.parts.read(),
            constructs=store.constructs.read(),
        )

    def _read_and_type(path):
        """Read a FASTA and type every record. Raises ValueError on bad input."""
        try:
            records = read_fasta(path)
        except OSError as exc:
            raise ValueError(f"Cannot read {path}: {exc}") from exc
        if not records:
            raise ValueError(f"No sequences found in {path}")
        # Keep each record beside its result: the approval view needs the
        # untrimmed sequence to show what would be cut away.
        return [(rec, process(rec.name, rec.seq)) for rec in records]

    def _trim_detail(record, result) -> dict:
        """Exactly which residues would be removed, for the approval view."""
        if not result.ok:
            return {}
        # process() ignores a trailing stop codon, so index against the same
        # string it indexed against.
        seq = record.seq.rstrip("*")
        return {
            "v_start": result.v_start + 1,  # 1-based for display
            "v_end": result.v_end + 1,
            "total": len(seq),
            "n_count": result.n_trimmed,
            "c_count": result.c_trimmed,
            "n_seq": seq[: result.v_start],
            "kept": result.clean_seq,
            "c_seq": seq[result.v_end + 1 :],
        }

    def _pair_groups(typed):
        """Group results by antibody, deciding which pairs are complete.

        Every record is grouped, passed or failed, so a pair that cannot be
        saved can say which half is missing or which half failed.
        """
        grouped: dict[str, list] = {}
        for _, result in typed:
            base, _suffix = parse_declared(result.name)
            grouped.setdefault(base, []).append(result)

        groups = []
        for base, group in sorted(grouped.items()):
            vh = next((r for r in group if r.ok and r.detected == "heavy"), None)
            vl = next((r for r in group if r.ok and r.detected == "light"), None)
            if vh and vl:
                groups.append({"id": base, "ready": True, "vh": vh, "vl": vl})
                continue
            failed = [r.name for r in group if not r.ok]
            if failed:
                reason = f"{', '.join(failed)} failed the sanity check"
            elif vh:
                reason = "no light chain (VL) in this file"
            elif vl:
                reason = "no heavy chain (VH) in this file"
            else:
                reason = "no usable heavy or light chain"
            groups.append({"id": base, "ready": False, "reason": reason})
        return groups

    @app.post("/api/check")
    def check():
        """Card 1, phase 1: type and trim, and report. Writes nothing.

        The result is a proposal the user approves in phase 2; nothing reaches
        the library until /api/approve is called with the ids to keep.
        """
        payload = request.get_json(silent=True) or {}
        try:
            path = _resolve(payload.get("path", ""))
            typed = _read_and_type(path)
        except ValueError as exc:
            return jsonify(ok=False, error=str(exc)), 400

        records = [
            {
                "name": result.name,
                "status": "ok" if result.ok else "fail",
                "declared": result.declared or "",
                "detected": result.detected or "",
                "subclass": result.subclass or "",
                "length": len(result.clean_seq),
                "trim_note": _trim_note(result) if result.ok else "",
                "trim": _trim_detail(record, result),
                "confidence": result.confidence if result.ok else "",
                "evidence": result.evidence,
                "warnings": result.warnings,
                "errors": result.errors,
            }
            for record, result in typed
        ]

        existing = store.antibodies.read_map()
        proposals = []
        for group in _pair_groups(typed):
            if not group["ready"]:
                proposals.append(
                    {
                        "id": group["id"],
                        "ready": False,
                        "action": "cannot add",
                        "detail": group["reason"],
                    }
                )
                continue
            vh, vl = group["vh"], group["vl"]
            prior = existing.get(group["id"])
            changed = bool(prior) and (
                prior.get("vh_seq") != vh.clean_seq
                or prior.get("vl_seq") != vl.clean_seq
            )
            proposals.append(
                {
                    "id": group["id"],
                    "ready": True,
                    "action": "replace" if prior else "add",
                    "changed": changed,
                    "detail": (
                        f"VH {len(vh.clean_seq)} aa, VL {len(vl.clean_seq)} aa "
                        f"({vl.subclass or 'subclass unknown'})"
                    ),
                    "subclass": vl.subclass or "",
                    "vh_len": len(vh.clean_seq),
                    "vl_len": len(vl.clean_seq),
                    "confidence": f"{vh.confidence} / {vl.confidence}",
                }
            )

        return jsonify(
            ok=True,
            file=str(path),
            records=records,
            proposals=proposals,
            summary={
                "passed": sum(1 for _, r in typed if r.ok),
                "failed": sum(1 for _, r in typed if not r.ok),
                "ready": sum(1 for p in proposals if p["ready"]),
            },
        )

    @app.post("/api/approve")
    def approve():
        """Card 1, phase 2: commit the approved antibodies to the library.

        The file is re-read and re-typed rather than trusting sequences sent
        back by the browser, so what is stored is always what the rules derive
        from the file on disk.
        """
        payload = request.get_json(silent=True) or {}
        approved = {str(i) for i in (payload.get("ids") or [])}
        if not approved:
            return jsonify(ok=False, error="Nothing approved."), 400

        try:
            path = _resolve(payload.get("path", ""))
            typed = _read_and_type(path)
        except ValueError as exc:
            return jsonify(ok=False, error=str(exc)), 400

        saved, skipped = [], []
        for group in _pair_groups(typed):
            if group["id"] not in approved:
                continue
            if not group["ready"]:
                skipped.append({"id": group["id"], "detail": group["reason"]})
                continue
            row = store.antibody_row(group["id"], group["vh"], group["vl"], path)
            outcome = store.antibodies.upsert(row)
            saved.append({"id": group["id"], "outcome": outcome})

        missing = approved - {s["id"] for s in saved} - {s["id"] for s in skipped}
        for antibody_id in sorted(missing):
            skipped.append({"id": antibody_id, "detail": "not found in this file"})

        return jsonify(
            ok=True,
            saved=saved,
            skipped=skipped,
            library=store.antibodies.read(),
        )

    @app.post("/api/parts/import")
    def import_parts():
        """Card 3: validate a parts FASTA and load it into the parts table."""
        payload = request.get_json(silent=True) or {}
        try:
            path = _resolve(payload.get("path", ""))
        except ValueError as exc:
            return jsonify(ok=False, error=str(exc)), 400

        parts, problems = assemble.load_parts(path)
        if problems:
            return jsonify(ok=False, error="; ".join(problems), problems=problems), 400

        outcome = store.import_parts(parts, path)
        return jsonify(
            ok=True,
            file=str(path),
            outcome=outcome,
            parts=store.parts.read(),
        )

    def _seq_diff(old: str, new: str) -> dict:
        """Compare a stored part against an incoming one."""
        if old == new:
            return {
                "status": "identical",
                "first_diff": None,
                "note": "identical to the stored sequence",
            }
        limit = min(len(old), len(new))
        first = next((i for i in range(limit) if old[i] != new[i]), limit)
        bits = []
        if len(old) != len(new):
            delta = len(new) - len(old)
            bits.append(
                f"incoming is {abs(delta)} residue(s) "
                f"{'longer' if delta > 0 else 'shorter'}"
            )
        bits.append(f"first differs at residue {first + 1}")
        return {"status": "differs", "first_diff": first + 1, "note": "; ".join(bits)}

    @app.post("/api/parts/preview")
    def preview_parts():
        """Card 3, phase 1: summarise what a parts file would import.

        Nothing is written. A part already in the table is compared against the
        incoming one so the user can decide per part whether to override it.
        """
        payload = request.get_json(silent=True) or {}
        try:
            path = _resolve(payload.get("path", ""))
        except ValueError as exc:
            return jsonify(ok=False, error=str(exc)), 400

        parts, problems = assemble.load_parts(path)
        if problems:
            return jsonify(ok=False, error="; ".join(problems), problems=problems), 400

        stored = store.parts.read_map()
        incoming = []
        for name, seq in parts.items():
            row = store.part_row(name, seq, path)
            entry = {
                "name": name,
                "role": row["role"],
                "subclass": row["subclass"],
                "incoming_len": row["length"],
                "incoming_seq": row["seq"],
            }
            prior = stored.get(name)
            if prior is None:
                entry.update(
                    status="new",
                    note="not in the parts table yet",
                    existing_len=None,
                    existing_seq="",
                    first_diff=None,
                )
            else:
                entry.update(
                    existing_len=int(prior.get("length") or 0),
                    existing_seq=prior.get("seq", ""),
                    **_seq_diff(prior.get("seq", ""), seq),
                )
            incoming.append(entry)

        incoming.sort(key=lambda e: (e["status"] != "differs", e["name"]))
        return jsonify(
            ok=True,
            file=str(path),
            parts=incoming,
            summary={
                "new": sum(1 for e in incoming if e["status"] == "new"),
                "identical": sum(1 for e in incoming if e["status"] == "identical"),
                "differs": sum(1 for e in incoming if e["status"] == "differs"),
            },
        )

    @app.post("/api/parts/apply")
    def apply_parts():
        """Card 3, phase 2: import, honouring a per-part override decision.

        `decisions` maps a part name to "override" or "keep". A part that is
        new is always added; one that is unchanged is a no-op; one that differs
        is only replaced when its decision says so, so the default of leaving a
        decision out is to keep what is already stored.
        """
        payload = request.get_json(silent=True) or {}
        decisions = {str(k): str(v) for k, v in (payload.get("decisions") or {}).items()}
        try:
            path = _resolve(payload.get("path", ""))
        except ValueError as exc:
            return jsonify(ok=False, error=str(exc)), 400

        parts, problems = assemble.load_parts(path)
        if problems:
            return jsonify(ok=False, error="; ".join(problems), problems=problems), 400

        rows = store.parts.read_map()
        applied = []
        for name, seq in parts.items():
            prior = rows.get(name)
            if prior is None:
                rows[name] = store.part_row(name, seq, path)
                applied.append({"name": name, "outcome": "added"})
            elif prior.get("seq", "") == seq:
                applied.append({"name": name, "outcome": "unchanged"})
            elif decisions.get(name) == "override":
                rows[name] = store.part_row(name, seq, path)
                applied.append({"name": name, "outcome": "overridden"})
            else:
                applied.append({"name": name, "outcome": "kept original"})

        store.parts.write_all(list(rows.values()))
        return jsonify(
            ok=True,
            file=str(path),
            applied=applied,
            changed=sum(
                1 for a in applied if a["outcome"] in ("added", "overridden")
            ),
            parts=store.parts.read(),
        )

    def _describe(construct, row=None) -> dict:
        """One construct rendered for the build-plan UI."""
        cl_segment = next(
            (s for s in construct.segments if s.kind == "light_constant"), None
        )
        seq = construct.seq
        return {
            "id": construct.name,
            "ok": construct.ok,
            "subclass": construct.subclass or "",
            "cl_part": cl_segment.label if cl_segment else "",
            "size_note": construct.size_note if construct.ok else "",
            "residue_count": construct.residue_count if construct.ok else 0,
            "has_stop": construct.has_stop,
            "layout": construct.layout,
            "segments": [
                {
                    "label": s.label,
                    "kind": s.kind,
                    "source": s.source,
                    "start": s.start,
                    "end": s.end,
                    "length": len(s.seq),
                }
                for s in construct.segments
            ],
            "preview": (
                f"{seq[:48]}…{seq[-24:]}" if len(seq) > 80 else seq
            ),
            "warnings": construct.warnings,
            "errors": construct.errors,
            "vh_confidence": (row or {}).get("vh_confidence", ""),
            "vl_confidence": (row or {}).get("vl_confidence", ""),
        }

    def _plan(ids, preset):
        """Resolve the selected library rows into constructs, without writing."""
        rows = store.antibodies.read_map()
        parts = store.parts_as_dict()
        slots = preset_store.slots_for(preset)
        planned, constructs = [], []
        for antibody_id in ids:
            row = rows.get(antibody_id)
            if row is None:
                planned.append(
                    {
                        "id": antibody_id,
                        "ok": False,
                        "errors": ["not in the antibody library"],
                        "warnings": [],
                        "segments": [],
                        "subclass": "",
                        "cl_part": "",
                        "size_note": "",
                        "preview": "",
                    }
                )
                continue
            construct = assemble.build_from_layout(
                antibody_id,
                row.get("vh_seq", ""),
                row.get("vl_seq", ""),
                row.get("light_subclass") or None,
                parts,
                preset,
                slots,
            )
            constructs.append(construct)
            planned.append(_describe(construct, row))
        return planned, constructs

    def _missing(slots, parts):
        """Parts a layout needs but the table lacks, with a near-match hint."""
        return [
            {"part": name, "suggestion": assemble.suggest_part(name, parts)}
            for name in assemble.missing_parts(slots, parts)
        ]

    def _palette(parts):
        """The components a layout can be built from, as draggable tags.

        Every part in the table gets its own tag except IGKC and IGLC2: those
        two share a single light-constant tag, because a layout picks between
        them per antibody rather than naming one outright.
        """
        light_constants = set(assemble.CL_FOR_SUBCLASS.values())
        items = [
            {
                "kind": "vl",
                "label": "VL",
                "detail": "light V domain, per antibody",
                "group": "Antibody",
            },
            {
                "kind": "vh",
                "label": "VH",
                "detail": "heavy V domain, per antibody",
                "group": "Antibody",
            },
            {
                "kind": "light_constant",
                "label": " / ".join(sorted(light_constants)),
                "detail": "resolved from the kappa/lambda call",
                "group": "Antibody",
            },
        ]
        for name, seq in parts.items():
            if name in light_constants:
                continue
            items.append(
                {
                    "kind": "part",
                    "part": name,
                    "label": name,
                    "detail": f"{len(seq.rstrip('*'))} aa",
                    "group": "Parts",
                }
            )
        items.append(
            {
                "kind": "stop",
                "label": "Stop",
                "detail": "stop codon, must come last",
                "group": "Terminator",
            }
        )
        return items

    @app.get("/api/presets")
    def presets():
        parts = store.parts_as_dict()
        return jsonify(
            default=assemble.DEFAULT_PRESET,
            presets={
                key: {
                    "label": value["label"],
                    "description": value["description"],
                    "builtin": value["builtin"],
                    "slots": value["slots"],
                    "layout": assemble.describe_layout(
                        parts=parts, slots=value["slots"]
                    ),
                    "missing": _missing(value["slots"], parts),
                }
                for key, value in preset_store.all().items()
            },
            palette=_palette(parts),
            parts_loaded=bool(parts),
        )

    @app.post("/api/presets/describe")
    def describe_draft():
        """Render a layout being edited, without saving it."""
        payload = request.get_json(silent=True) or {}
        slots = payload.get("slots") or []
        parts = store.parts_as_dict()
        problems = assemble.validate_slots(slots) if slots else []
        return jsonify(
            ok=True,
            problems=problems,
            missing=_missing(slots, parts),
            layout=assemble.describe_layout(parts=parts, slots=slots),
        )

    @app.post("/api/presets/save")
    def save_preset():
        payload = request.get_json(silent=True) or {}
        label = (payload.get("label") or "").strip()
        description = (payload.get("description") or "").strip()
        slots = payload.get("slots") or []
        key = (payload.get("key") or "").strip() or slugify(label)

        if not label:
            return jsonify(ok=False, error="Give the preset a name."), 400

        problems = assemble.validate_slots(slots)
        if problems:
            return (
                jsonify(
                    ok=False,
                    error="; ".join(problems),
                    problems=problems,
                ),
                400,
            )

        try:
            outcome = preset_store.save(key, label, description, slots, utc_now())
        except ValueError as exc:
            return jsonify(ok=False, error=str(exc)), 400

        return jsonify(ok=True, key=key, outcome=outcome)

    @app.post("/api/presets/delete")
    def delete_preset():
        payload = request.get_json(silent=True) or {}
        key = (payload.get("key") or "").strip()
        if not key:
            return jsonify(ok=False, error="No preset named."), 400
        try:
            removed = preset_store.delete(key)
        except ValueError as exc:
            return jsonify(ok=False, error=str(exc)), 400
        if not removed:
            return jsonify(ok=False, error=f"No saved preset '{key}'."), 400
        return jsonify(ok=True, key=key)

    @app.post("/api/plan")
    def plan():
        """Card 4 -> 5: show what would be assembled, resolving IGKC/IGLC2."""
        payload = request.get_json(silent=True) or {}
        ids = [str(i) for i in (payload.get("ids") or [])]
        preset = payload.get("preset") or assemble.DEFAULT_PRESET

        if not ids:
            return jsonify(ok=False, error="No antibodies selected."), 400
        if preset_store.slots_for(preset) is None:
            return jsonify(ok=False, error=f"Unknown preset {preset!r}."), 400
        if not store.parts_as_dict():
            return (
                jsonify(
                    ok=False,
                    error="The parts table is empty — open Utilities \u2192 Construct parts first.",
                ),
                400,
            )

        missing = _missing(preset_store.slots_for(preset), store.parts_as_dict())
        if missing:
            detail = "; ".join(
                m["part"] + (f" (did you mean '{m['suggestion']}'?)" if m["suggestion"] else "")
                for m in missing
            )
            return (
                jsonify(
                    ok=False,
                    error=(
                        f"This preset needs parts that are not in the parts "
                        f"table: {detail}. Open Utilities to import a parts "
                        f"FASTA containing them, or edit the preset there."
                    ),
                    missing=missing,
                ),
                400,
            )

        planned, _ = _plan(ids, preset)
        return jsonify(
            ok=True,
            preset=preset,
            preset_label=preset_store.all()[preset]["label"],
            constructs=planned,
            buildable=sum(1 for p in planned if p["ok"]),
        )

    @app.post("/api/build")
    def build_constructs():
        """Card 5: write the confirmed constructs to a FASTA file."""
        payload = request.get_json(silent=True) or {}
        ids = [str(i) for i in (payload.get("ids") or [])]
        preset = payload.get("preset") or assemble.DEFAULT_PRESET
        raw_output = (payload.get("output") or "bcr_constructs.fasta").strip()

        if not ids:
            return jsonify(ok=False, error="No antibodies selected."), 400
        if preset_store.slots_for(preset) is None:
            return jsonify(ok=False, error=f"Unknown preset {preset!r}."), 400

        output = Path(raw_output).expanduser()
        if not output.is_absolute():
            output = Path.cwd() / output
        if output.is_dir():
            return jsonify(ok=False, error=f"{output} is a directory."), 400

        planned, constructs = _plan(ids, preset)
        built = [c for c in constructs if c.ok]
        if not built:
            return (
                jsonify(
                    ok=False,
                    error="Nothing can be built from the current selection.",
                    constructs=planned,
                ),
                400,
            )

        # Say what is being displaced: the file is replaced, not appended to.
        replaced = 0
        if output.exists():
            try:
                replaced = len(read_fasta(output))
            except OSError:
                replaced = 0

        try:
            write_fasta(
                output,
                [
                    Record(
                        name=f"{c.name}_BCR",
                        seq=c.seq,
                        description=(
                            f"{c.subclass} BCR, {c.size_note} | {c.layout}"
                        ),
                    )
                    for c in built
                ],
            )
        except OSError as exc:
            return jsonify(ok=False, error=f"Cannot write {output}: {exc}"), 400

        # Keep what was built, so the constructs library reflects the file.
        stored = store.constructs.read_map()
        for construct in built:
            row = store.construct_row(construct, preset, output)
            stored[row["id"]] = row
        store.constructs.write_all(list(stored.values()))

        return jsonify(
            ok=True,
            output=str(output),
            written=len(built),
            replaced=replaced,
            skipped=[p["id"] for p in planned if not p["ok"]],
            constructs=planned,
            library=store.constructs.read(),
        )

    @app.get("/api/constructs")
    def constructs():
        return jsonify(constructs=store.constructs.read())

    @app.post("/api/construct/delete")
    def delete_construct():
        payload = request.get_json(silent=True) or {}
        construct_id = (payload.get("id") or "").strip()
        if not construct_id:
            return jsonify(ok=False, error="No id given."), 400
        removed = store.constructs.delete(construct_id)
        return jsonify(ok=removed, library=store.constructs.read())

    @app.post("/api/constructs/export")
    def export_constructs():
        """Write selected stored constructs back out to a FASTA."""
        payload = request.get_json(silent=True) or {}
        ids = [str(i) for i in (payload.get("ids") or [])]
        raw_output = (payload.get("output") or "bcr_constructs.fasta").strip()
        if not ids:
            return jsonify(ok=False, error="No constructs selected."), 400

        output = Path(raw_output).expanduser()
        if not output.is_absolute():
            output = Path.cwd() / output
        if output.is_dir():
            return jsonify(ok=False, error=f"{output} is a directory."), 400

        rows = store.constructs.read_map()
        missing = [i for i in ids if i not in rows]
        if missing:
            return (
                jsonify(ok=False, error=f"Not in the library: {', '.join(missing)}"),
                400,
            )

        replaced = 0
        if output.exists():
            try:
                replaced = len(read_fasta(output))
            except OSError:
                replaced = 0

        try:
            write_fasta(
                output,
                [
                    Record(
                        name=rows[i]["id"],
                        seq=rows[i]["seq"],
                        description=(
                            f"{rows[i]['subclass']} BCR, "
                            f"{rows[i]['residue_count']} aa"
                            + (" + stop" if rows[i]["has_stop"] == "yes" else "")
                            + f" | {rows[i]['layout']}"
                        ),
                    )
                    for i in ids
                ],
            )
        except OSError as exc:
            return jsonify(ok=False, error=f"Cannot write {output}: {exc}"), 400

        return jsonify(
            ok=True, output=str(output), written=len(ids), replaced=replaced
        )

    @app.post("/api/antibody/delete")
    def delete_antibody():
        payload = request.get_json(silent=True) or {}
        antibody_id = (payload.get("id") or "").strip()
        if not antibody_id:
            return jsonify(ok=False, error="No id given."), 400
        removed = store.antibodies.delete(antibody_id)
        return jsonify(ok=removed, library=store.antibodies.read())

    @app.errorhandler(500)
    def internal_error(exc):  # pragma: no cover - diagnostic path
        return jsonify(ok=False, error=f"Server error: {exc}"), 500

    return app


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="bcr.web", description="Local web UI for the BCR constructor."
    )
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument(
        "--host",
        default="127.0.0.1",
        help="bind address; leave as localhost unless you know what you are doing",
    )
    parser.add_argument("--data-dir", default=None)
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args(argv)

    app = create_app(args.data_dir)
    print(f"BCR constructor UI: http://{args.host}:{args.port}")
    app.run(host=args.host, port=args.port, debug=args.debug)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

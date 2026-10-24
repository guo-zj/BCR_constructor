"""Local web UI for the BCR constructor.

    python3 -m bcr.web

Binds to 127.0.0.1 only. The app reads FASTA files by path from the machine it
runs on, so it is a local tool and must not be exposed to a network.

The sanity-check card drives bcr.chains.process -- the same function the
`clean` CLI calls -- so the web UI and the command line cannot drift apart.
"""

import argparse
import re
import traceback
from datetime import datetime, timezone
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


#: Extensions accepted by the drop zones. A FASTA by any other name is not
#: one we will guess at -- the user can still type a path to it.
UPLOAD_SUFFIXES = (".fa", ".fasta")

#: A FASTA of antibody parts is a few kilobytes. This is a generous ceiling
#: that still keeps a stray drop of something enormous from being written.
MAX_UPLOAD_BYTES = 8 * 1024 * 1024


def _safe_upload_name(raw: str) -> str:
    """Reduce a browser-supplied filename to a safe basename.

    The name arrives from the client, so it is untrusted: any directory part
    is dropped and the stem is reduced to a conservative character set. The
    extension gates the upload rather than merely decorating it.
    """
    stem = Path(raw or "").name
    suffix = Path(stem).suffix.lower()
    if suffix not in UPLOAD_SUFFIXES:
        accepted = " or ".join(UPLOAD_SUFFIXES)
        raise ValueError(f"Only {accepted} files can be dropped here (got {stem!r}).")
    base = re.sub(r"[^A-Za-z0-9._-]+", "_", Path(stem).stem).strip("._-")
    return (base[:60] or "upload") + suffix


#: A name that goes into a FASTA header and a CSV key: no whitespace (a FASTA
#: name is the first token on its line), no comma, nothing exotic.
NAME_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:+-]*$")


def _name_problems(name: str, kind: str) -> list[str]:
    if not name:
        return [f"A {kind} needs a name."]
    if not NAME_PATTERN.match(name):
        return [
            f"'{name}' is not a usable {kind} name. It becomes a FASTA header "
            f"and a CSV key, so it must start with a letter or digit and hold "
            f"only letters, digits and . _ : + -"
        ]
    if len(name) > 80:
        return [f"That {kind} name is longer than 80 characters."]
    return []


def _clean_residues(raw: str) -> str:
    """Turn pasted text into a bare sequence.

    Accepts a plain sequence, a wrapped one, or a whole FASTA record pasted in
    -- header lines, line breaks, spaces and residue numbering all come out.
    Anything left that is not an amino acid is reported by the caller's rules,
    not silently dropped.
    """
    body = "\n".join(
        line for line in (raw or "").splitlines() if not line.lstrip().startswith(">")
    )
    return re.sub(r"[\s\d]+", "", body).upper()


def _subclass_of(part_name: str) -> str:
    """Which light subclass a light-constant part serves."""
    for subclass, name in assemble.CL_FOR_SUBCLASS.items():
        if name == part_name:
            return subclass
    return ""


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
    app.config["MAX_CONTENT_LENGTH"] = MAX_UPLOAD_BYTES
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
            parts=store.all_part_rows(),
            constructs=store.constructs.read(),
            reserved=sorted(assemble.RESERVED_PART_NAMES),
            shadowed=store.shadowed_parts(),
            stale_derived=store.stale_derived(),
            isotypes=[
                {"key": k, "label": assemble.ISOTYPE_LABEL[k], "part": v}
                for k, v in assemble.CH_FOR_ISOTYPE.items()
            ],
        )

    @app.post("/api/library/derive")
    def derive_library():
        """Recompute the derived whole-chain columns for every antibody.

        Offered rather than done on read: it rewrites antibodies.csv, and a
        write is the user's to approve. Nothing is written when nothing is
        stale.
        """
        changed = store.refresh_derived()
        return jsonify(
            ok=True,
            changed=changed,
            library=store.antibodies.read(),
            stale_derived=store.stale_derived(),
        )

    @app.post("/api/parts/unshadow")
    def unshadow_parts():
        """Delete stored part rows whose names have become intrinsic."""
        payload = request.get_json(silent=True) or {}
        names = [str(n) for n in (payload.get("names") or [])]
        if not names:
            return jsonify(ok=False, error="No parts named."), 400
        shadowed = {row["name"] for row in store.shadowed_parts()}
        unknown = [n for n in names if n not in shadowed]
        if unknown:
            return jsonify(
                ok=False,
                error=f"Not a shadowed row: {', '.join(unknown)}",
            ), 400
        removed = [n for n in names if store.parts.delete(n)]
        return jsonify(
            ok=True,
            removed=removed,
            parts=store.all_part_rows(),
            shadowed=store.shadowed_parts(),
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
                # A trim is a modification, so it is offered rather than
                # applied: nothing is stored for this record until the user
                # says trim or keep. See /api/approve.
                "proposes_trim": bool(
                    result.ok and (result.n_trimmed or result.c_trimmed)
                ),
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
            decisions = [
                {
                    "record": result.name,
                    "chain": tag,
                    "n": result.n_trimmed,
                    "c": result.c_trimmed,
                    "note": _trim_note(result),
                    "confidence": result.confidence,
                    "as_given": len(result.input_seq),
                    "trimmed_to": len(result.clean_seq),
                }
                for tag, result in (("VH", vh), ("VL", vl))
                if result.n_trimmed or result.c_trimmed
            ]
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
                    "decisions": decisions,
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
                "awaiting": sum(
                    len(p.get("decisions") or []) for p in proposals if p["ready"]
                ),
            },
        )

    @app.post("/api/approve")
    def approve():
        """Card 1, phase 2: commit the approved antibodies to the library.

        The file is re-read and re-typed rather than trusting sequences sent
        back by the browser, so what is stored is always what the rules derive
        from the file on disk.

        `decisions` maps a record name to "trim" or "keep". Every record the
        rules would shorten needs one: with no decision the pair is skipped
        and said to be waiting, rather than stored under a guess. Records that
        need no trim need no decision.
        """
        payload = request.get_json(silent=True) or {}
        approved = {str(i) for i in (payload.get("ids") or [])}
        raw_decisions = payload.get("decisions") or {}
        if not isinstance(raw_decisions, dict):
            return jsonify(ok=False, error="`decisions` must be an object."), 400
        decisions = {str(k): str(v) for k, v in raw_decisions.items()}
        unknown = sorted({v for v in decisions.values()} - {"trim", "keep"})
        if unknown:
            return jsonify(
                ok=False,
                error=f"A decision must be 'trim' or 'keep', not {unknown[0]!r}.",
            ), 400
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

            keep, undecided = [], []
            for tag, result in (("VH", group["vh"]), ("VL", group["vl"])):
                if not (result.n_trimmed or result.c_trimmed):
                    continue
                choice = decisions.get(result.name)
                if choice is None:
                    undecided.append(
                        f"{result.name} ({_trim_note(result)}, "
                        f"{result.confidence} confidence)"
                    )
                elif choice == "keep":
                    keep.append(tag)
            if undecided:
                skipped.append({
                    "id": group["id"],
                    "detail": "waiting on a trim decision for "
                              + ", ".join(undecided),
                    "awaiting": True,
                })
                continue

            row = store.antibody_row(
                group["id"], group["vh"], group["vl"], path, keep=keep
            )
            outcome = store.antibodies.upsert(row)
            saved.append({"id": group["id"], "outcome": outcome, "kept": keep})

        missing = approved - {s["id"] for s in saved} - {s["id"] for s in skipped}
        for antibody_id in sorted(missing):
            skipped.append({"id": antibody_id, "detail": "not found in this file"})

        return jsonify(
            ok=True,
            saved=saved,
            skipped=skipped,
            library=store.antibodies.read(),
            stale_derived=store.stale_derived(),
        )

    @app.post("/api/parts/import")
    def import_parts():
        """Card 3: validate a parts FASTA and load it into the parts table."""
        payload = request.get_json(silent=True) or {}
        try:
            path = _resolve(payload.get("path", ""))
        except ValueError as exc:
            return jsonify(ok=False, error=str(exc)), 400

        parts, problems, reserved = assemble.read_parts_file(path)
        if problems:
            return jsonify(ok=False, error="; ".join(problems), problems=problems), 400

        outcome = store.import_parts(parts, path)
        return jsonify(
            ok=True,
            file=str(path),
            outcome=outcome,
            parts=store.all_part_rows(),
        )

    def _reserved_entries(reserved: dict) -> list[dict]:
        """The records in a parts file that claim an intrinsic name.

        Each is offered a new name rather than being dropped: the user chose
        that sequence deliberately, they just cannot have that name. A record
        identical to the intrinsic sequence needs no rename, since it would
        override nothing -- it is simply redundant.
        """
        entries = []
        for name, seq in sorted(reserved.items()):
            same = seq == assemble.INTRINSIC_SEQS[name]
            entries.append({
                "name": name,
                "seq": seq,
                "length": len(seq),
                "identical": same,
                "intrinsic_length": len(assemble.INTRINSIC_SEQS[name]),
                "suggestion": "" if same else f"{name}_mine",
                "note": (
                    "identical to the intrinsic part, so there is nothing to "
                    "import"
                    if same
                    else "an intrinsic name cannot be reused; give it another"
                ),
            })
        return entries

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

        parts, problems, reserved = assemble.read_parts_file(path)
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
        reserved_entries = _reserved_entries(reserved)
        return jsonify(
            ok=True,
            file=str(path),
            parts=incoming,
            reserved=reserved_entries,
            summary={
                "new": sum(1 for e in incoming if e["status"] == "new"),
                "identical": sum(1 for e in incoming if e["status"] == "identical"),
                "differs": sum(1 for e in incoming if e["status"] == "differs"),
                "reserved": sum(1 for e in reserved_entries if not e["identical"]),
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
        renames = {str(k): str(v).strip() for k, v in (payload.get("renames") or {}).items()}
        try:
            path = _resolve(payload.get("path", ""))
        except ValueError as exc:
            return jsonify(ok=False, error=str(exc)), 400

        parts, problems, reserved = assemble.read_parts_file(path)
        if problems:
            return jsonify(ok=False, error="; ".join(problems), problems=problems), 400

        rows = store.parts.read_map()
        applied = []

        # A record claiming an intrinsic name is imported only under a new
        # name the user supplies. Without one it is skipped and said to be
        # skipped -- it is never imported under the reserved name, and never
        # quietly dropped either.
        for name, seq in sorted(reserved.items()):
            if seq == assemble.INTRINSIC_SEQS[name]:
                applied.append({
                    "name": name,
                    "outcome": "already intrinsic",
                    "detail": "identical to the intrinsic part",
                })
                continue
            new_name = renames.get(name, "")
            if not new_name:
                applied.append({
                    "name": name,
                    "outcome": "skipped",
                    "detail": f"'{name}' is an intrinsic name; give it another "
                              f"to import this sequence",
                })
                continue
            trouble = _name_problems(new_name, "part")
            if new_name in assemble.RESERVED_PART_NAMES:
                trouble = [f"'{new_name}' is also an intrinsic name."]
            if not trouble:
                _, trouble = assemble.validate_part(new_name, seq)
            if trouble:
                applied.append({
                    "name": name, "outcome": "skipped", "detail": trouble[0]
                })
                continue
            if new_name in rows and rows[new_name].get("seq", "") != seq:
                applied.append({
                    "name": name,
                    "outcome": "skipped",
                    "detail": f"'{new_name}' is already stored with a different "
                              f"sequence; pick another name",
                })
                continue
            rows[new_name] = store.part_row(new_name, seq, path)
            applied.append({
                "name": name,
                "outcome": "renamed",
                "detail": f"imported as '{new_name}'",
            })
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
            parts=store.all_part_rows(),
            shadowed=store.shadowed_parts(),
        )

    def _describe(construct, row=None) -> dict:
        """One construct rendered for the build-plan UI."""
        cl_segment = next(
            (s for s in construct.segments if s.kind == "light_constant"), None
        )
        cl_part = cl_segment.label if cl_segment else ""
        if not cl_part and any(s.kind == "light_chain" for s in construct.segments):
            # An express layout folds the constant region into the light
            # chain, so name it from the subclass instead of the segment.
            cl_part = assemble.CL_FOR_SUBCLASS.get(construct.subclass or "", "")
        isotypes = sorted(
            {
                assemble.ISOTYPE_LABEL[iso]
                for s in construct.segments
                if s.kind == "heavy_chain"
                for iso, part in assemble.CH_FOR_ISOTYPE.items()
                if f"+{part})" in s.label
            }
        )
        seq = construct.seq
        return {
            "id": construct.name,
            "ok": construct.ok,
            "subclass": construct.subclass or "",
            "cl_part": cl_part,
            "isotype": " / ".join(isotypes),
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
            {
                "kind": "light_chain",
                "label": "Light chain",
                "detail": "VL + IGKC / IGLC2, per antibody",
                "group": "Whole chains",
            },
        ]
        for isotype, part in assemble.CH_FOR_ISOTYPE.items():
            items.append(
                {
                    "kind": "heavy_chain",
                    "isotype": isotype,
                    "label": f"Heavy chain {assemble.ISOTYPE_LABEL[isotype]}",
                    "detail": f"VH + {part}, per antibody",
                    "group": "Whole chains",
                }
            )
        for name, seq in parts.items():
            if name in light_constants:
                continue
            items.append(
                {
                    "kind": "part",
                    "part": name,
                    "label": name,
                    "detail": f"{len(seq.rstrip('*'))} aa",
                    "group": (
                        "Intrinsic"
                        if name in assemble.RESERVED_PART_NAMES
                        else "Parts"
                    ),
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
                    "express": bool(value.get("express")),
                    "slots": value["slots"],
                    "layout": assemble.describe_layout(
                        parts=parts, slots=value["slots"]
                    ),
                    "missing": _missing(value["slots"], parts),
                }
                for key, value in preset_store.all().items()
            },
            palette=_palette(parts),
            # The intrinsic parts are always present, so this says whether the
            # user has imported any of their own -- which is what the UI hint
            # is actually about.
            parts_loaded=bool(store.parts.read()),
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
        return jsonify(
            ok=removed,
            library=store.antibodies.read(),
            stale_derived=store.stale_derived(),
        )

    # --- editing the tables by hand -------------------------------------
    #
    # Hand-entered sequences go through exactly the same typing and part rules
    # as imported ones: the editor is a second way in, not a way around. A
    # record that cannot be typed is refused with the reasons, rather than
    # stored as an untyped row the rest of the app cannot reason about.

    @app.post("/api/antibody/save")
    def save_antibody():
        payload = request.get_json(silent=True) or {}
        new_id = (payload.get("id") or "").strip()
        old_id = (payload.get("old_id") or "").strip()
        vh_raw = _clean_residues(payload.get("vh_seq") or "")
        vl_raw = _clean_residues(payload.get("vl_seq") or "")
        user_notes = (payload.get("user_notes") or "").strip()

        problems = _name_problems(new_id, "antibody")
        upper = new_id.upper()
        if upper.endswith("_VH") or upper.endswith("_VL"):
            problems.append(
                "The id is the antibody name on its own; drop the "
                "_VH / _VL suffix, which says which chain a FASTA record is."
            )
        if not vh_raw:
            problems.append("No VH sequence given.")
        if not vl_raw:
            problems.append("No VL sequence given.")
        if problems:
            return jsonify(ok=False, error=problems[0], problems=problems), 400

        rows = store.antibodies.read_map()
        if old_id and old_id not in rows:
            return jsonify(ok=False, error=f"'{old_id}' is not in the library."), 404
        if new_id != old_id and new_id in rows:
            return jsonify(
                ok=False,
                error=f"'{new_id}' is already in the library. "
                f"Rename that one first, or edit it instead.",
            ), 409

        vh = process(f"{new_id}_VH", vh_raw)
        vl = process(f"{new_id}_VL", vl_raw)
        failures = [f"VH: {e}" for e in vh.errors] + [f"VL: {e}" for e in vl.errors]
        if failures:
            return jsonify(
                ok=False,
                error="These sequences do not type as a VH / VL pair.",
                problems=failures,
            ), 400

        # Trimming is a modification, so it is never applied to a pasted
        # sequence on the strength of the paste alone. A chain the rules would
        # shorten comes back as a decision to make, with nothing written.
        keep, pending = [], []
        for tag, result, field in (
            ("VH", vh, "vh_decision"), ("VL", vl, "vl_decision")
        ):
            if not (result.n_trimmed or result.c_trimmed):
                continue
            choice = (payload.get(field) or "").strip()
            if choice not in ("trim", "keep"):
                pending.append({
                    "chain": tag,
                    "field": field,
                    "n": result.n_trimmed,
                    "c": result.c_trimmed,
                    "note": _trim_note(result),
                    "confidence": result.confidence,
                    "as_given": len(result.input_seq),
                    "trimmed_to": len(result.clean_seq),
                    "n_seq": result.input_seq[: result.v_start],
                    "c_seq": result.input_seq[result.v_end + 1 :],
                })
            elif choice == "keep":
                keep.append(tag)
        if pending:
            return jsonify(ok=True, stored=False, needs_decision=pending)

        row = store.antibody_row(new_id, vh, vl, "manual entry", keep=keep)
        row["user_notes"] = user_notes
        renaming = bool(old_id) and old_id != new_id
        if renaming:
            store.antibodies.delete(old_id)
        outcome = store.antibodies.upsert(row)
        if renaming:
            outcome = "renamed"   # upsert sees a fresh key and would say "added"

        # A construct records which antibody it was built from, so a rename
        # has to follow through or the reference dangles.
        renamed_constructs = 0
        if old_id and old_id != new_id:
            constructs = store.constructs.read()
            for construct in constructs:
                if construct.get("antibody_id") == old_id:
                    construct["antibody_id"] = new_id
                    renamed_constructs += 1
            if renamed_constructs:
                store.constructs.write_all(constructs)

        notes = []
        for tag, result in (("VH", vh), ("VL", vl)):
            if result.n_trimmed or result.c_trimmed:
                notes.append(
                    f"{tag} kept as given; the {_trim_note(result)} the rules "
                    f"proposed was declined, so this is not a bare V domain"
                    if tag in keep
                    else f"{tag} trimmed to its V domain ({_trim_note(result)})"
                )
            if result.confidence != "high":
                notes.append(f"{tag} boundary confidence is {result.confidence}")
        if renamed_constructs:
            notes.append(
                f"{renamed_constructs} construct"
                f"{'' if renamed_constructs == 1 else 's'} repointed at '{new_id}'"
            )

        return jsonify(
            ok=True,
            outcome=outcome,
            id=new_id,
            notes=notes,
            warnings=list(vh.warnings) + list(vl.warnings),
            library=store.antibodies.read(),
            constructs=store.constructs.read(),
            stale_derived=store.stale_derived(),
        )

    @app.post("/api/part/save")
    def save_part():
        payload = request.get_json(silent=True) or {}
        new_name = (payload.get("name") or "").strip()
        old_name = (payload.get("old_name") or "").strip()
        seq = _clean_residues(payload.get("seq") or "")
        user_notes = (payload.get("user_notes") or "").strip()

        problems = _name_problems(new_name, "part")
        if not seq:
            problems.append("No sequence given.")
        if problems:
            return jsonify(ok=False, error=problems[0], problems=problems), 400

        # Renaming a light constant away would leave nothing for a
        # light_constant slot to resolve to, for any layout at all.
        for candidate in (new_name, old_name):
            if candidate in assemble.RESERVED_PART_NAMES:
                return jsonify(
                    ok=False,
                    error=f"'{candidate}' is an intrinsic part: it ships with "
                    f"the package and every layout resolves to it, so it "
                    f"cannot be edited, renamed or shadowed. Reserved names "
                    f"are {', '.join(sorted(assemble.RESERVED_PART_NAMES))}.",
                    reserved=sorted(assemble.RESERVED_PART_NAMES),
                ), 400

        residues, trouble = assemble.validate_part(new_name, seq)
        if trouble:
            return jsonify(ok=False, error=trouble[0], problems=trouble), 400

        rows = store.parts.read_map()
        if old_name and old_name not in rows:
            return jsonify(ok=False, error=f"'{old_name}' is not a stored part."), 404
        if new_name != old_name and new_name in rows:
            return jsonify(
                ok=False, error=f"A part called '{new_name}' is already stored."
            ), 409

        row = store.part_row(new_name, residues, "manual entry")
        row["user_notes"] = user_notes

        notes = []
        renaming = bool(old_name) and old_name != new_name
        if renaming:
            store.parts.delete(old_name)
            touched = preset_store.rename_part(old_name, new_name)
            if touched:
                notes.append(
                    f"{len(touched)} saved preset"
                    f"{'' if len(touched) == 1 else 's'} now point at "
                    f"'{new_name}': {', '.join(touched)}"
                )
            builtin = sorted(
                key
                for key, value in assemble.PRESETS.items()
                if any(s.get("part") == old_name for s in value["slots"])
            )
            if builtin:
                notes.append(
                    f"built-in preset{'' if len(builtin) == 1 else 's'} "
                    f"{', '.join(builtin)} still ask for '{old_name}', which no "
                    f"longer exists — those will report a missing part"
                )
        outcome = store.parts.upsert(row)
        if renaming:
            outcome = "renamed"   # upsert sees a fresh key and would say "added"

        return jsonify(
            ok=True, outcome=outcome, name=new_name, role=row["role"],
            length=row["length"], notes=notes, parts=store.all_part_rows(),
            shadowed=store.shadowed_parts(),
        )

    @app.post("/api/part/delete")
    def delete_part():
        payload = request.get_json(silent=True) or {}
        name = (payload.get("name") or "").strip()
        if not name:
            return jsonify(ok=False, error="No name given."), 400
        if name in assemble.RESERVED_PART_NAMES:
            return jsonify(
                ok=False,
                error=f"'{name}' is an intrinsic part and cannot be deleted: "
                f"it ships with the package and layouts resolve to it.",
                reserved=sorted(assemble.RESERVED_PART_NAMES),
            ), 400
        used_by = sorted(
            key
            for key, value in preset_store.all().items()
            if any(s.get("part") == name for s in value["slots"])
        )
        removed = store.parts.delete(name)
        notes = []
        if removed and used_by:
            notes.append(
                f"still named by {len(used_by)} layout"
                f"{'' if len(used_by) == 1 else 's'} ({', '.join(used_by)}), "
                f"which will report a missing part until it is edited"
            )
        return jsonify(
            ok=removed, notes=notes, parts=store.all_part_rows(),
            shadowed=store.shadowed_parts(),
        )

    @app.post("/api/construct/save")
    def save_construct():
        """A construct is derived, so only its id and notes are editable."""
        payload = request.get_json(silent=True) or {}
        new_id = (payload.get("id") or "").strip()
        old_id = (payload.get("old_id") or "").strip()
        user_notes = (payload.get("user_notes") or "").strip()

        problems = _name_problems(new_id, "construct")
        if not old_id:
            problems.append(
                "A construct is assembled, not typed in: build one from the "
                "Assembly layout card rather than adding it here."
            )
        if problems:
            return jsonify(ok=False, error=problems[0], problems=problems), 400

        rows = store.constructs.read_map()
        if old_id not in rows:
            return jsonify(ok=False, error=f"'{old_id}' is not in the library."), 404
        if new_id != old_id and new_id in rows:
            return jsonify(
                ok=False, error=f"A construct called '{new_id}' is already stored."
            ), 409

        row = dict(rows[old_id])
        row["id"] = new_id
        row["user_notes"] = user_notes
        row["updated_utc"] = utc_now()
        if new_id != old_id:
            store.constructs.delete(old_id)
        store.constructs.upsert(row)
        return jsonify(ok=True, id=new_id, library=store.constructs.read())

    @app.post("/api/upload")
    def upload():
        """Save a dropped FASTA to data/uploads/ and hand back its path.

        An upload is deliberately only a way of getting a file onto disk. The
        path it returns is fed into the ordinary path-driven flows, so the
        two-phase invariant holds unchanged: /api/approve still re-reads and
        re-types the file from disk rather than trusting anything the browser
        sent.
        """
        upload = request.files.get("file")
        if upload is None:
            return jsonify(ok=False, error="No file was uploaded."), 400
        try:
            name = _safe_upload_name(upload.filename or "")
        except ValueError as exc:
            return jsonify(ok=False, error=str(exc)), 400

        data = upload.read(MAX_UPLOAD_BYTES + 1)
        if len(data) > MAX_UPLOAD_BYTES:
            mb = MAX_UPLOAD_BYTES // (1024 * 1024)
            return jsonify(ok=False, error=f"File is larger than {mb} MB."), 400
        if not data.strip():
            return jsonify(ok=False, error="That file is empty."), 400
        if not data.lstrip().startswith(b">"):
            return jsonify(
                ok=False,
                error="That does not look like a FASTA file: it has no '>' header line.",
            ), 400

        dest_dir = store.data_dir / "uploads"
        dest_dir.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        dest = dest_dir / f"{stamp}-{name}"
        n = 2
        while dest.exists():
            dest = dest_dir / f"{stamp}-{n}-{name}"
            n += 1
        dest.write_bytes(data)
        return jsonify(ok=True, path=str(dest), name=name, bytes=len(data))

    @app.errorhandler(413)
    def too_large(exc):
        mb = MAX_UPLOAD_BYTES // (1024 * 1024)
        return jsonify(ok=False, error=f"Request body is larger than {mb} MB."), 413

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

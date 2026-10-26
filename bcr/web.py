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
from .chains import (
    LIGHT_SUBCLASSES,
    accept_as_given,
    can_accept_as_given,
    parse_declared,
    process,
)
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


def _add_anyway(result) -> dict:
    """Whether a failed record could be stored as given, for the check view.

    `add_anyway` gates the button; `add_anyway_blocked` is the reason it is
    not offered, which is shown rather than swallowed -- a user who insists
    deserves to know why insisting will not help.
    """
    if result.ok:
        return {"add_anyway": False, "add_anyway_blocked": ""}
    chain, reason = can_accept_as_given(result)
    return {
        "add_anyway": bool(chain),
        "add_anyway_chain": chain or "",
        "add_anyway_blocked": "" if chain else reason,
        # A light chain cannot be stored without a kappa/lambda call, and an
        # overruled record has no boundary for the rules to read one from, so
        # the user states it. The UI offers the two outright.
        "add_anyway_subclasses": (
            list(LIGHT_SUBCLASSES) if chain == "light" else []
        ),
    }


def _overrides_from(payload) -> tuple[list, dict]:
    """The override fields of a check/approve payload, validated.

    Raises ValueError, which both endpoints already turn into a 400.
    """
    overrides = payload.get("overrides") or []
    if not isinstance(overrides, list):
        raise ValueError("`overrides` must be a list.")
    raw = payload.get("subclasses") or {}
    if not isinstance(raw, dict):
        raise ValueError("`subclasses` must be an object.")
    subclasses = {str(k): str(v).strip().lower() for k, v in raw.items()}
    unknown = sorted(set(subclasses.values()) - set(LIGHT_SUBCLASSES))
    if unknown:
        raise ValueError(
            f"A light chain is kappa or lambda, not {unknown[0]!r}."
        )
    return overrides, subclasses


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
            legacy_constructs=store.legacy_constructs(),
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

    def _read_and_type(path, overrides=(), subclasses=None):
        """Read a FASTA and type every record. Raises ValueError on bad input.

        `overrides` names records the user has chosen to store even though the
        rules rejected them. Those are accepted verbatim here, at the one
        place typing happens, so pairing, the approval view and the stored row
        all see the same thing -- there is no second path into the library.

        `subclasses` maps a record name to "kappa" or "lambda", which a light
        chain stored as given cannot do without.
        """
        try:
            records = read_fasta(path)
        except OSError as exc:
            raise ValueError(f"Cannot read {path}: {exc}") from exc
        if not records:
            raise ValueError(f"No sequences found in {path}")
        insisted = {str(name) for name in overrides}
        stated = subclasses or {}
        # Keep each record beside its result: the approval view needs the
        # untrimmed sequence to show what would be cut away.
        typed = []
        for rec in records:
            result = process(rec.name, rec.seq)
            if rec.name in insisted and not result.ok:
                # ValueError if the failure cannot be overruled, or if a light
                # chain arrives without its kappa/lambda call.
                result = accept_as_given(result, stated.get(rec.name))
            typed.append((rec, result))
        unknown = sorted(insisted - {rec.name for rec, _ in typed})
        if unknown:
            raise ValueError(
                f"Asked to store {unknown[0]!r} as given, but there is no such "
                f"record in {path}"
            )
        return typed

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
            overrides, subclasses = _overrides_from(payload)
            path = _resolve(payload.get("path", ""))
            typed = _read_and_type(path, overrides, subclasses)
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
                # A failure is a judgement, and the user may overrule it: the
                # record is then stored exactly as it came in. `add_anyway`
                # says whether this one could be, and `add_anyway_blocked`
                # why not when it could not.
                "accepted_as_given": result.accepted_as_given,
                **_add_anyway(result),
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
                    # Which halves of this pair are only here because the user
                    # overruled the sanity check. The pair is addable, but it
                    # should not look like any other pair while it is.
                    "accepted_as_given": [
                        tag
                        for tag, result in (("VH", vh), ("VL", vl))
                        if result.accepted_as_given
                    ],
                }
            )

        return jsonify(
            ok=True,
            file=str(path),
            records=records,
            proposals=proposals,
            summary={
                # An overruled record is addable but it did not pass, and the
                # count must not say it did.
                "passed": sum(1 for _, r in typed if r.ok and not r.accepted_as_given),
                "failed": sum(1 for _, r in typed if not r.ok),
                "as_given": sum(1 for _, r in typed if r.accepted_as_given),
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

        `overrides` names records the rules rejected that the user has chosen
        to store anyway, and `subclasses` gives the kappa/lambda call for any
        light chain among them. They are re-typed and re-accepted here from
        the file on disk, exactly as in /api/check, so insisting still cannot
        smuggle a sequence past the rules -- it only changes what is done
        about them.
        """
        payload = request.get_json(silent=True) or {}
        approved = {str(i) for i in (payload.get("ids") or [])}
        try:
            overrides, subclasses = _overrides_from(payload)
        except ValueError as exc:
            return jsonify(ok=False, error=str(exc)), 400
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
            typed = _read_and_type(path, overrides, subclasses)
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
            saved.append({
                "id": group["id"],
                "outcome": outcome,
                "kept": keep,
                "as_given": [
                    tag
                    for tag, result in (("VH", group["vh"]), ("VL", group["vl"]))
                    if result.accepted_as_given
                ],
            })

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

    def _describe(construct, row=None, preset="") -> dict:
        """One construct rendered for the build-plan UI.

        `preset` is the layout it was built under. It travels with every
        construct because one run can build the same antibody several ways,
        and the name it will be written under depends on it.
        """
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
            # The name this construct will carry in the FASTA and the table.
            "construct_id": store.construct_id(construct.name, preset)
            if preset else "",
            "preset": preset,
            "preset_label": (
                preset_store.all().get(preset, {}).get("label", preset)
                if preset else ""
            ),
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

    def _plan(ids, presets):
        """Resolve the selected library rows into constructs, without writing.

        Every antibody is planned under every layout asked for, so one run can
        produce a whole matrix. Returns (planned, built) where `built` pairs
        each construct with the preset it came from -- the preset is not
        recoverable from the construct, and the name depends on it.
        """
        rows = store.antibodies.read_map()
        parts = store.parts_as_dict()
        planned, constructs = [], []
        for preset in presets:
            slots = preset_store.slots_for(preset)
            label = preset_store.all().get(preset, {}).get("label", preset)
            for antibody_id in ids:
                row = rows.get(antibody_id)
                if row is None:
                    planned.append(
                        {
                            "id": antibody_id,
                            "construct_id": "",
                            "preset": preset,
                            "preset_label": label,
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
                constructs.append((preset, construct))
                planned.append(_describe(construct, row, preset))
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

    def _presets_from(payload) -> list[str]:
        """The layouts a plan or build asks for, in order, de-duplicated.

        One run can name several, so `presets` is the general form; `preset`
        is the one-layout case, and an empty payload means the shipped
        default. Raises ValueError, which both callers turn into a 400.
        """
        raw = payload.get("presets")
        if raw is None:
            raw = [payload.get("preset") or assemble.DEFAULT_PRESET]
        if not isinstance(raw, list):
            raise ValueError("`presets` must be a list.")
        wanted, seen = [], set()
        for key in (str(k) for k in raw):
            if key not in seen:
                seen.add(key)
                wanted.append(key)
        if not wanted:
            raise ValueError("No layout selected.")
        unknown = [k for k in wanted if preset_store.slots_for(k) is None]
        if unknown:
            raise ValueError(f"Unknown preset {unknown[0]!r}.")
        return wanted

    def _missing_for(presets) -> dict:
        """Which of these layouts cannot be completed from the parts table."""
        parts = store.parts_as_dict()
        return {
            key: _missing(preset_store.slots_for(key), parts)
            for key in presets
            if _missing(preset_store.slots_for(key), parts)
        }

    @app.post("/api/plan")
    def plan():
        """Card 4 -> 5: show what would be assembled, resolving IGKC/IGLC2.

        Every selected antibody is planned under every selected layout, so a
        plan is a matrix: the response is one flat list, each entry naming the
        layout it belongs to and the name it would be written under.
        """
        payload = request.get_json(silent=True) or {}
        ids = [str(i) for i in (payload.get("ids") or [])]
        try:
            presets = _presets_from(payload)
        except ValueError as exc:
            return jsonify(ok=False, error=str(exc)), 400

        if not ids:
            return jsonify(ok=False, error="No antibodies selected."), 400
        # A layout that cannot be completed stops the whole plan rather than
        # quietly contributing nothing: the user picked it on purpose.
        shortfall = _missing_for(presets)
        if shortfall:
            labels = preset_store.all()
            detail = "; ".join(
                f"{labels[key]['label']} needs "
                + ", ".join(
                    m["part"]
                    + (f" (did you mean '{m['suggestion']}'?)" if m["suggestion"] else "")
                    for m in missing
                )
                for key, missing in shortfall.items()
            )
            return (
                jsonify(
                    ok=False,
                    error=(
                        f"Parts are missing from the parts table: {detail}. "
                        f"Open Utilities to import a parts FASTA containing "
                        f"them, or edit the layout there."
                    ),
                    missing=[m for missing in shortfall.values() for m in missing],
                    missing_by_preset=shortfall,
                ),
                400,
            )

        planned, _ = _plan(ids, presets)
        all_presets = preset_store.all()
        return jsonify(
            ok=True,
            presets=[
                {
                    "key": key,
                    "label": all_presets[key]["label"],
                    "buildable": sum(
                        1 for p in planned if p["preset"] == key and p["ok"]
                    ),
                    "planned": sum(1 for p in planned if p["preset"] == key),
                }
                for key in presets
            ],
            # The single-layout view, kept because most runs are one layout.
            preset=presets[0] if len(presets) == 1 else "",
            preset_label=(
                all_presets[presets[0]]["label"] if len(presets) == 1 else ""
            ),
            constructs=planned,
            buildable=sum(1 for p in planned if p["ok"]),
        )

    @app.post("/api/build")
    def build_constructs():
        """Card 5: write the confirmed constructs to a FASTA file."""
        payload = request.get_json(silent=True) or {}
        ids = [str(i) for i in (payload.get("ids") or [])]
        raw_output = (payload.get("output") or "bcr_constructs.fasta").strip()
        try:
            presets = _presets_from(payload)
        except ValueError as exc:
            return jsonify(ok=False, error=str(exc)), 400

        if not ids:
            return jsonify(ok=False, error="No antibodies selected."), 400

        output = Path(raw_output).expanduser()
        if not output.is_absolute():
            output = Path.cwd() / output
        if output.is_dir():
            return jsonify(ok=False, error=f"{output} is a directory."), 400

        planned, constructs = _plan(ids, presets)
        built = [(preset, c) for preset, c in constructs if c.ok]
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

        labels = preset_store.all()
        try:
            write_fasta(
                output,
                [
                    Record(
                        # The build mode is part of the name: several layouts
                        # of one antibody share this file, and a reader with
                        # only the FASTA must still be able to tell them apart.
                        name=store.construct_id(c.name, preset),
                        seq=c.seq,
                        description=(
                            f"{c.subclass} BCR, {c.size_note} | "
                            f"{labels[preset]['label']} | {c.layout}"
                        ),
                    )
                    for preset, c in built
                ],
            )
        except OSError as exc:
            return jsonify(ok=False, error=f"Cannot write {output}: {exc}"), 400

        # Keep what was built, so the constructs library reflects the file.
        stored = store.constructs.read_map()
        for preset, construct in built:
            row = store.construct_row(construct, preset, output)
            stored[row["id"]] = row
        store.constructs.write_all(list(stored.values()))

        return jsonify(
            ok=True,
            output=str(output),
            written=len(built),
            replaced=replaced,
            skipped=[
                {"id": p["id"], "preset": p["preset"], "label": p["preset_label"]}
                for p in planned
                if not p["ok"]
            ],
            # One line per layout, so a run of several says what each produced.
            by_preset=[
                {
                    "key": key,
                    "label": labels[key]["label"],
                    "written": sum(1 for preset, _ in built if preset == key),
                }
                for key in presets
            ],
            names=[store.construct_id(c.name, preset) for preset, c in built],
            constructs=planned,
            library=store.constructs.read(),
            legacy=store.legacy_constructs(),
        )

    @app.get("/api/constructs")
    def constructs():
        return jsonify(
            constructs=store.constructs.read(),
            legacy=store.legacy_constructs(),
        )

    @app.post("/api/constructs/rename-legacy")
    def rename_legacy_constructs():
        """Give pre-layout-naming rows their layout-qualified ids.

        Offered rather than done: it renames rows the user may have referred
        to elsewhere, so it waits to be asked. The sequences are untouched.
        """
        payload = request.get_json(silent=True) or {}
        ids = [str(i) for i in (payload.get("ids") or [])]
        if not ids:
            ids = [entry["id"] for entry in store.legacy_constructs()]
        if not ids:
            return jsonify(ok=False, error="No rows need renaming."), 400
        moved = store.rename_legacy_constructs(ids)
        return jsonify(
            ok=True,
            moved=moved,
            library=store.constructs.read(),
            legacy=store.legacy_constructs(),
        )

    @app.post("/api/construct/delete")
    def delete_construct():
        payload = request.get_json(silent=True) or {}
        construct_id = (payload.get("id") or "").strip()
        if not construct_id:
            return jsonify(ok=False, error="No id given."), 400
        removed = store.constructs.delete(construct_id)
        return jsonify(
            ok=removed,
            library=store.constructs.read(),
            legacy=store.legacy_constructs(),
        )

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
                            + f" | {rows[i].get('preset', '')}"
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

        typed = {
            "VH": process(f"{new_id}_VH", vh_raw),
            "VL": process(f"{new_id}_VL", vl_raw),
        }
        # A failed chain can be stored as it stands, if the user says so for
        # that chain by name. The flag is per chain, so insisting on one half
        # never quietly waves the other half through. For the light chain the
        # flag carries the kappa/lambda call, since nothing can be read off a
        # sequence whose V domain was never found.
        failures, offers = [], []
        for tag, field in (("VH", "vh_as_given"), ("VL", "vl_as_given")):
            result = typed[tag]
            if result.ok:
                continue
            chain, reason = can_accept_as_given(result)
            said = payload.get(field)
            if chain and said:
                stated = said if isinstance(said, str) else ""
                try:
                    typed[tag] = accept_as_given(result, stated)
                    continue
                except ValueError as exc:
                    # Said yes, but not which subclass: ask rather than guess.
                    reason = str(exc)
                    failures.append(f"{tag}: {reason}")
                    offers.append({
                        "chain": tag,
                        "field": field,
                        "can": True,
                        "blocked": "",
                        "as_given": len(result.input_seq),
                        "subclasses": list(LIGHT_SUBCLASSES)
                        if chain == "light" else [],
                    })
                    continue
            failures.extend(f"{tag}: {e}" for e in result.errors)
            offers.append({
                "chain": tag,
                "field": field,
                "can": bool(chain),
                "blocked": "" if chain else reason,
                "as_given": len(result.input_seq),
                # The light chain is offered as two: kappa or lambda.
                "subclasses": list(LIGHT_SUBCLASSES) if chain == "light" else [],
            })
        if failures:
            return jsonify(
                ok=False,
                error="These sequences do not type as a VH / VL pair.",
                problems=failures,
                add_anyway=offers,
            ), 400
        vh, vl = typed["VH"], typed["VL"]

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
            if result.accepted_as_given:
                notes.append(
                    f"{tag} stored as given at your request after failing the "
                    f"sanity check: all {len(result.clean_seq)} residues, with "
                    f"no chain class, boundary or FR4 established"
                )
                continue
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
            legacy_constructs=store.legacy_constructs(),
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
        return jsonify(
            ok=True,
            id=new_id,
            library=store.constructs.read(),
            legacy=store.legacy_constructs(),
        )

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

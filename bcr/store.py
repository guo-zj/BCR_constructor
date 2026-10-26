"""CSV-backed store for clean antibody V regions and construct parts.

Two tables, both plain CSV so they can be opened in a spreadsheet, diffed, or
read by any other tool without this package:

    data/antibodies.csv   one row per antibody, keyed on its name
    data/parts.csv        one row per construct part, keyed on its name

Sequences are stored alongside the rule-derived facts that were established
when they were typed -- chain class, kappa/lambda, boundary confidence, how
much was trimmed -- so a downstream step never has to re-derive them, and a
low-confidence call stays visible instead of being forgotten at import time.

Writes are atomic: the whole table is rendered to a temporary file in the same
directory and moved into place, so an interrupted write cannot truncate an
existing table.
"""

import csv
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from .assemble import derived_chains

DEFAULT_DATA_DIR = Path(__file__).resolve().parents[1] / "data"

ANTIBODY_FIELDS = [
    "id",
    "vh_seq",
    "vl_seq",
    "vh_len",
    "vl_len",
    "light_subclass",
    "vh_confidence",
    "vl_confidence",
    "vh_trim",
    "vl_trim",
    "vh_fr4",
    "vl_fr4",
    # The three whole chains the pair implies, spliced with the intrinsic
    # constant regions by assemble.derived_chains. Materialised here so the
    # CSV is usable on its own -- they are derived, never a source of truth,
    # and are recomputed from vh_seq / vl_seq on every write.
    "light_chain",
    "heavy_chain_igg1",
    "heavy_chain_igm",
    "source_file",
    "updated_utc",
    "notes",
    # `notes` is derived -- it carries the warnings typing produced.
    # `user_notes` is free text the user owns, so an edit in the UI never
    # overwrites a rule-derived warning. Both columns are appended last, so a
    # table written by an older version still reads.
    "user_notes",
]

# Parts hold residues only. A stop codon is contributed by the layout at build
# time (assemble's `stop` slot), so there is no per-part stop column.
PART_FIELDS = [
    "name",
    "seq",
    "length",
    "role",
    "subclass",
    "source_file",
    "updated_utc",
    "user_notes",
]

CONSTRUCT_FIELDS = [
    "id",
    "antibody_id",
    "preset",
    "subclass",
    "cl_part",
    "seq",
    "residue_count",
    "has_stop",
    "layout",
    "segments_json",
    "source_file",
    "updated_utc",
    "notes",
    "user_notes",
]

# A part's role is inferred from its sequence by bcr.assemble.infer_role,
# because part names come from the parts FASTA and are not fixed.


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


class CsvTable:
    """A CSV file treated as a table keyed on one column."""

    def __init__(self, path, fieldnames: list[str], key: str):
        self.path = Path(path)
        self.fieldnames = fieldnames
        self.key = key

    def exists(self) -> bool:
        return self.path.exists()

    def read(self) -> list[dict]:
        if not self.path.exists():
            return []
        with open(self.path, newline="", encoding="utf-8") as fh:
            return [dict(row) for row in csv.DictReader(fh)]

    def read_map(self) -> dict[str, dict]:
        return {row[self.key]: row for row in self.read() if row.get(self.key)}

    def write_all(self, rows: list[dict]) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        fd, tmp = tempfile.mkstemp(
            dir=self.path.parent, prefix=self.path.name, suffix=".tmp"
        )
        try:
            with os.fdopen(fd, "w", newline="", encoding="utf-8") as fh:
                writer = csv.DictWriter(fh, fieldnames=self.fieldnames)
                writer.writeheader()
                for row in rows:
                    writer.writerow({k: row.get(k, "") for k in self.fieldnames})
            os.replace(tmp, self.path)
        except BaseException:
            if os.path.exists(tmp):
                os.unlink(tmp)
            raise

    def upsert(self, row: dict) -> str:
        """Insert or replace one row. Returns "added" or "updated"."""
        rows = self.read()
        key_value = row[self.key]
        for index, existing in enumerate(rows):
            if existing.get(self.key) == key_value:
                rows[index] = row
                self.write_all(rows)
                return "updated"
        rows.append(row)
        self.write_all(rows)
        return "added"

    def delete(self, key_value: str) -> bool:
        rows = self.read()
        kept = [r for r in rows if r.get(self.key) != key_value]
        if len(kept) == len(rows):
            return False
        self.write_all(kept)
        return True


class Store:
    """Both tables, rooted at one data directory."""

    def __init__(self, data_dir=None):
        self.data_dir = Path(data_dir or DEFAULT_DATA_DIR)
        self.antibodies = CsvTable(
            self.data_dir / "antibodies.csv", ANTIBODY_FIELDS, "id"
        )
        self.parts = CsvTable(self.data_dir / "parts.csv", PART_FIELDS, "name")
        self.constructs = CsvTable(
            self.data_dir / "constructs.csv", CONSTRUCT_FIELDS, "id"
        )

    # --- antibodies ------------------------------------------------------

    def antibody_row(self, antibody_id, vh, vl, source_file, keep=()) -> dict:
        """Build a row from two validated ChainResults.

        `vh` and `vl` are bcr.chains.ChainResult objects that have already
        passed typing; the caller is responsible for that check.

        `keep` names the chains ("VH", "VL") whose proposed trim the user
        declined. Those are stored exactly as they came in, and the row says
        so: a sequence is never shortened without someone having said yes.

        A result carrying `accepted_as_given` failed typing and is being
        stored on the user's insistence (chains.accept_as_given). Nothing was
        derived for it, so the columns that would otherwise be rule-derived
        say nothing rather than saying something untrue.
        """
        keep = set(keep)
        notes = []
        for tag, result in (("VH", vh), ("VL", vl)):
            for warning in result.warnings:
                notes.append(f"{tag}: {warning}")
            if tag in keep and (result.n_trimmed or result.c_trimmed):
                notes.append(
                    f"{tag}: kept as given at your request -- the "
                    f"{_trim_note(result)} the rules proposed is still present, "
                    f"so this is not a bare V domain"
                )

        def stored(tag, result):
            return result.input_seq if tag in keep else result.clean_seq

        def trim_note(tag, result):
            if result.accepted_as_given:
                return "stored as given (failed the sanity check)"
            if tag not in keep:
                return _trim_note(result)
            return (
                f"kept as given ({_trim_note(result)} declined)"
                if result.n_trimmed or result.c_trimmed
                else "none"
            )

        vh_seq, vl_seq = stored("VH", vh), stored("VL", vl)
        chains = derived_chains(vh_seq, vl_seq, vl.subclass)
        if not chains["light_chain"]:
            notes.append(
                "no whole light chain could be derived: the kappa/lambda call "
                "is unknown, so neither IGKC nor IGLC2 can be chosen"
            )
        if keep:
            notes.append(
                "the derived whole chains splice a sequence that was kept as "
                "given, so they are not a V domain plus a constant region"
            )
        unverified = [
            tag for tag, r in (("VH", vh), ("VL", vl)) if r.accepted_as_given
        ]
        if unverified:
            notes.append(
                f"{' and '.join(unverified)} stored as given after failing the "
                f"sanity check: no chain class, boundary or FR4 was established, "
                f"so the derived whole chains splice an unverified sequence"
            )

        def fr4(tag, result, back):
            """FR4 as observed -- or nothing, if no boundary was ever found."""
            if result.accepted_as_given:
                return ""
            return result.input_seq[max(0, result.v_end - back) : result.v_end + 1]

        return {
            "id": antibody_id,
            "vh_seq": vh_seq,
            "vl_seq": vl_seq,
            "vh_len": len(vh_seq),
            "vl_len": len(vl_seq),
            "light_subclass": vl.subclass or "",
            "vh_confidence": vh.confidence,
            "vl_confidence": vl.confidence,
            "vh_trim": trim_note("VH", vh),
            "vl_trim": trim_note("VL", vl),
            # FR4 as actually observed: 11 residues for heavy (WGxGTxVTVSS),
            # 10 for light (FGxGTxVxxx). Identifies the J gene family and
            # records that the boundary landed where the rules say it should.
            # Taken from the detected V end rather than the tail of what is
            # stored, so it is still FR4 when a C-terminal trim was declined.
            "vh_fr4": fr4("VH", vh, 10),
            "vl_fr4": fr4("VL", vl, 9),
            "light_chain": chains["light_chain"],
            "heavy_chain_igg1": chains["heavy_chain_igg1"],
            "heavy_chain_igm": chains["heavy_chain_igm"],
            "source_file": str(source_file),
            "updated_utc": utc_now(),
            "notes": "; ".join(notes),
        }

    # --- parts -----------------------------------------------------------

    def part_row(self, name, seq, source_file) -> dict:
        """Build a parts row. A trailing stop codon is never stored."""
        from . import assemble

        residues = seq.rstrip("*")
        role, subclass = assemble.infer_role(residues)
        return {
            "name": name,
            "seq": residues,
            "length": len(residues),
            "role": role,
            "subclass": subclass,
            "source_file": str(source_file),
            "updated_utc": utc_now(),
        }

    def import_parts(self, parts: dict[str, str], source_file) -> dict[str, str]:
        """Upsert every part. Returns {name: "added"|"updated"}."""
        outcome = {}
        rows = self.parts.read_map()
        for name, seq in parts.items():
            row = self.part_row(name, seq, source_file)
            outcome[name] = "updated" if name in rows else "added"
            rows[name] = row
        self.parts.write_all(list(rows.values()))
        return outcome

    def intrinsic_part_rows(self) -> list[dict]:
        """The intrinsic parts as parts-table rows, marked as such.

        They are never written to parts.csv: a layout resolves to them, so
        they have to be beyond editing, and the only way to guarantee that is
        not to keep an editable copy. They are merged in at read time instead,
        the same way built-in presets are layered under the user's.
        """
        from . import assemble

        return [
            {
                "name": name,
                "seq": value["seq"],
                "length": len(value["seq"]),
                "role": value["role"],
                "subclass": value["subclass"],
                "source_file": value["source"],
                "updated_utc": "",
                "user_notes": value["description"],
                "builtin": "yes",
            }
            for name, value in assemble.BUILTIN_PARTS.items()
        ]

    def all_part_rows(self) -> list[dict]:
        """Intrinsic parts first, then the user's own.

        A stored row whose name is now intrinsic is left out: the intrinsic
        sequence is what every layout resolves to, so showing both would
        suggest a choice that does not exist. Those rows are reported by
        `shadowed_parts` instead, rather than being deleted on the quiet.
        """
        from . import assemble

        user = [
            dict(row, builtin="")
            for row in self.parts.read()
            if row["name"] not in assemble.RESERVED_PART_NAMES
        ]
        return self.intrinsic_part_rows() + user

    def shadowed_parts(self) -> list[dict]:
        """Stored rows whose names have become intrinsic.

        These predate the intrinsic parts: they are no longer used by anything
        and can be removed, but that is the user's call, so they are reported
        rather than cleaned up. `same` says whether the stored sequence still
        matches the intrinsic one.
        """
        from . import assemble

        return [
            {
                "name": row["name"],
                "length": row.get("length", ""),
                "same": row["seq"] == assemble.INTRINSIC_SEQS[row["name"]],
            }
            for row in self.parts.read()
            if row["name"] in assemble.RESERVED_PART_NAMES
        ]

    def parts_as_dict(self) -> dict[str, str]:
        """The parts table in the shape bcr.assemble.build expects.

        Intrinsic parts are included and win over a same-named user row, which
        import and editing both refuse to create in any case.
        """
        from . import assemble

        # Intrinsic first, and a reserved name in the CSV is skipped rather
        # than merged: a hand-edited table must not be able to shadow one.
        merged = dict(assemble.INTRINSIC_SEQS)
        for row in self.parts.read():
            if row["name"] not in assemble.RESERVED_PART_NAMES:
                merged[row["name"]] = row["seq"]
        return merged

    def stale_derived(self) -> list[str]:
        """Antibody ids whose derived whole-chain columns are out of date."""
        stale = []
        for row in self.antibodies.read():
            chains = derived_chains(
                row.get("vh_seq", ""),
                row.get("vl_seq", ""),
                row.get("light_subclass", ""),
            )
            if any(row.get(key, "") != value for key, value in chains.items()):
                stale.append(row["id"])
        return stale

    def refresh_derived(self) -> list[str]:
        """Recompute the derived columns for every antibody. Returns the ids
        that changed. Rewrites nothing when nothing is stale."""
        rows = self.antibodies.read()
        changed = []
        for row in rows:
            chains = derived_chains(
                row.get("vh_seq", ""),
                row.get("vl_seq", ""),
                row.get("light_subclass", ""),
            )
            if any(row.get(key, "") != value for key, value in chains.items()):
                row.update(chains)
                changed.append(row["id"])
        if changed:
            self.antibodies.write_all(rows)
        return changed

    # --- constructs ------------------------------------------------------

    @staticmethod
    def construct_id(antibody_id: str, preset: str) -> str:
        """The row key for one antibody built under one layout.

        The same name the FASTA record carries: assemble.construct_name owns
        the rule, so the table and the file cannot disagree.
        """
        from . import assemble

        return assemble.construct_name(antibody_id, preset)

    def legacy_constructs(self) -> list[dict]:
        """Rows keyed before the layout was part of the name.

        Their id is `<antibody>_BCR`, which cannot tell two layouts of one
        antibody apart. Each row already records the layout it was built
        under, so the name it should have is computable -- but renaming a row
        is a modification, so these are reported and renamed only when asked.
        A row whose preset is blank cannot be mended this way and is left out.
        """
        out = []
        for row in self.constructs.read():
            antibody_id = row.get("antibody_id", "")
            preset = row.get("preset", "")
            if not antibody_id or not preset:
                continue
            if row["id"] == self.construct_id(antibody_id, preset):
                continue
            if row["id"] != f"{antibody_id}_BCR":
                continue   # a name the user chose; not ours to tidy
            out.append({
                "id": row["id"],
                "new_id": self.construct_id(antibody_id, preset),
                "antibody_id": antibody_id,
                "preset": preset,
            })
        return out

    def rename_legacy_constructs(self, ids) -> list[dict]:
        """Give the named legacy rows their layout-qualified ids.

        Returns what moved. A rename that would collide with a row already
        there is skipped, not merged: the newer row is the built one.
        """
        wanted = {str(i) for i in ids}
        rows = self.constructs.read()
        taken = {row["id"] for row in rows}
        moved = []
        for entry in self.legacy_constructs():
            if entry["id"] not in wanted or entry["new_id"] in taken:
                continue
            for row in rows:
                if row["id"] == entry["id"]:
                    row["id"] = entry["new_id"]
                    taken.add(entry["new_id"])
                    moved.append(entry)
        if moved:
            self.constructs.write_all(rows)
        return moved

    def construct_row(self, construct, preset, source_file) -> dict:
        """Build a row from an assembled bcr.assemble.Construct."""
        from . import assemble

        cl_segment = next(
            (s for s in construct.segments if s.kind == "light_constant"), None
        )
        cl_part = cl_segment.label if cl_segment else ""
        if not cl_part and any(s.kind == "light_chain" for s in construct.segments):
            # An express layout folds the constant region into the whole light
            # chain, so there is no segment of its own to read it off; the
            # subclass names it just as well.
            cl_part = assemble.CL_FOR_SUBCLASS.get(construct.subclass or "", "")
        return {
            "id": self.construct_id(construct.name, preset),
            "antibody_id": construct.name,
            "preset": preset,
            "subclass": construct.subclass or "",
            "cl_part": cl_part,
            "seq": construct.seq,
            "residue_count": construct.residue_count,
            "has_stop": "yes" if construct.has_stop else "no",
            # `layout` is for reading in a spreadsheet; `segments_json` is the
            # machine-readable form, so coordinates survive the round trip.
            "layout": construct.layout,
            "segments_json": json.dumps(
                [
                    {
                        "label": s.label,
                        "kind": s.kind,
                        "source": s.source,
                        "start": s.start,
                        "end": s.end,
                        "length": len(s.seq),
                    }
                    for s in construct.segments
                ]
            ),
            "source_file": str(source_file),
            "updated_utc": utc_now(),
            "notes": "; ".join(construct.warnings),
        }


def _trim_note(result) -> str:
    bits = []
    if result.n_trimmed:
        bits.append(f"-{result.n_trimmed}N")
    if result.c_trimmed:
        bits.append(f"-{result.c_trimmed}C")
    return "/".join(bits) if bits else "none"

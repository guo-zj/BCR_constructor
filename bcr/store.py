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
    "source_file",
    "updated_utc",
    "notes",
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

    def antibody_row(self, antibody_id, vh, vl, source_file) -> dict:
        """Build a row from two validated ChainResults.

        `vh` and `vl` are bcr.chains.ChainResult objects that have already
        passed typing; the caller is responsible for that check.
        """
        notes = []
        for tag, result in (("VH", vh), ("VL", vl)):
            for warning in result.warnings:
                notes.append(f"{tag}: {warning}")

        return {
            "id": antibody_id,
            "vh_seq": vh.clean_seq,
            "vl_seq": vl.clean_seq,
            "vh_len": len(vh.clean_seq),
            "vl_len": len(vl.clean_seq),
            "light_subclass": vl.subclass or "",
            "vh_confidence": vh.confidence,
            "vl_confidence": vl.confidence,
            "vh_trim": _trim_note(vh),
            "vl_trim": _trim_note(vl),
            # FR4 as actually observed: 11 residues for heavy (WGxGTxVTVSS),
            # 10 for light (FGxGTxVxxx). Identifies the J gene family and
            # records that the boundary landed where the rules say it should.
            "vh_fr4": vh.clean_seq[-11:],
            "vl_fr4": vl.clean_seq[-10:],
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

    def parts_as_dict(self) -> dict[str, str]:
        """The parts table in the shape bcr.assemble.build expects."""
        return {row["name"]: row["seq"] for row in self.parts.read()}

    # --- constructs ------------------------------------------------------

    def construct_row(self, construct, preset, source_file) -> dict:
        """Build a row from an assembled bcr.assemble.Construct."""
        cl_segment = next(
            (s for s in construct.segments if s.kind == "light_constant"), None
        )
        return {
            "id": f"{construct.name}_BCR",
            "antibody_id": construct.name,
            "preset": preset,
            "subclass": construct.subclass or "",
            "cl_part": cl_segment.label if cl_segment else "",
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

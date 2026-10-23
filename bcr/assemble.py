"""Assemble a membrane-bound BCR construct from trimmed V domains.

Layout, N- to C-terminus:

    Signal | VL | IGKC or IGLC2 | Linker | VH | CHandCT

The light constant region is chosen from the VL's own kappa/lambda call, so
the two halves of the light chain always match. Section 6 of the knowledge
base is what sets the first junction: a signal peptide butts directly against
the mature N-terminus with no linker, because signal peptidase reads the local
sequence and padding it would leave GS... on the mature protein.
"""

import json
from dataclasses import dataclass, field
from pathlib import Path

from . import rules
from .chains import ChainResult, parse_declared

# Which light constant region each subclass takes. These two names are fixed:
# a layout cannot name one outright, because the choice is made per antibody
# from its own kappa/lambda call, so the code has to know what to look for.
CL_FOR_SUBCLASS = {"kappa": "IGKC", "lambda": "IGLC2"}

# Every other part name comes from the parts FASTA and is not known in advance,
# so only the light constants are required.
REQUIRED_PARTS = tuple(CL_FOR_SUBCLASS.values())

# The first residues a part must carry, where we can identify it by name.
# The knowledge base (section 6) is explicit that a constant region starts at
# its own first residue; checking that catches a part file that has been
# truncated or pasted mid-domain.
EXPECTED_START = {
    "IGKC": rules.CL_KAPPA_START[:8],  # RTVAAPSV
    "IGLC2": rules.CL_LAMBDA_START[:9],  # GQPKAAPSV
}

# Since part names are free-form, a part's role is read from its sequence
# rather than from what it is called. These openings come from the knowledge
# base (sections 2, 3c and 5).
ROLE_SIGNATURES = (
    (rules.CL_KAPPA_START[:8], "light_constant", "kappa"),
    (rules.CL_LAMBDA_START[:9], "light_constant", "lambda"),
    (rules.CH1_IGG1_START, "heavy_constant", ""),
)

# A signal peptide is short and begins at Met. Nothing else in a construct
# part does both, so this is a safe fallback; anything unrecognised is "other"
# rather than guessed at.
SIGNAL_MAX_LEN = 40


def infer_role(seq: str) -> tuple[str, str]:
    """(role, subclass) for a part, read from its sequence, not its name."""
    for prefix, role, subclass in ROLE_SIGNATURES:
        if seq.startswith(prefix):
            return role, subclass
    if seq.startswith("M") and len(seq) <= SIGNAL_MAX_LEN:
        return "signal_peptide", ""
    return "other", ""


def missing_parts(slots, parts) -> list[str]:
    """Part names a layout needs that the parts table does not hold."""
    names: list[str] = []
    uses_light_constant = False
    for slot in slots:
        if not isinstance(slot, dict):
            continue
        if slot.get("kind") == "part":
            name = slot.get("part")
            if name and name not in parts and name not in names:
                names.append(name)
        elif slot.get("kind") == "light_constant":
            uses_light_constant = True
    if uses_light_constant:
        # Both are needed: which one is used depends on the antibody.
        for name in CL_FOR_SUBCLASS.values():
            if name not in parts and name not in names:
                names.append(name)
    return names


def suggest_part(name: str, parts) -> str | None:
    """The closest part name actually present, to catch a misspelling."""
    import difflib

    matches = difflib.get_close_matches(name, list(parts), n=1, cutoff=0.6)
    return matches[0] if matches else None

# --- layouts --------------------------------------------------------------
#
# A layout is an ordered list of slots, so the order of a construct is data
# rather than control flow: the UI can show it before anything is built, and a
# new format needs a new entry here rather than a new code path.
#
# Slot kinds:
#   part            a fixed sequence from the parts table, named outright
#   vl / vh         that antibody's own trimmed V domain
#   light_constant  resolved per antibody to IGKC or IGLC2 from its subclass
#   stop            the stop codon, contributed by the layout
#
# The stop codon is a slot rather than a property of any part: stored parts
# hold residues only, so whichever parts a layout uses, exactly one '*' is
# appended and it is always the last character.

SLOT_KINDS = ("part", "vl", "vh", "light_constant", "stop")
STOP = "*"

# The layouts that ship with the package, as data rather than code. The part
# names in them must match headers in the parts FASTA; if the parts table does
# not hold them, a layout reports them as missing (see `missing_parts`) rather
# than building something incomplete.
#
# This file holds the defaults and travels with the code. Layouts the user
# builds are kept separately in data/presets.json and layered on top of these
# by bcr.presets.PresetStore.
BUILTIN_PRESET_FILE = Path(__file__).with_name("builtin_presets.json")


def load_builtin_presets(path=None) -> tuple[dict, str]:
    """Read the shipped layouts. Returns (presets, default key).

    A missing or malformed file is a broken installation rather than a
    recoverable condition, so it is raised instead of being papered over --
    silently falling back to no layouts would turn a packaging fault into a
    puzzling empty preset list.
    """
    path = Path(path or BUILTIN_PRESET_FILE)
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
    except FileNotFoundError as exc:
        raise RuntimeError(f"built-in preset file is missing: {path}") from exc
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"built-in preset file is not valid JSON: {path}: {exc}") from exc

    presets = data.get("presets")
    default = data.get("default")
    if not isinstance(presets, dict) or not presets:
        raise RuntimeError(f"{path} defines no presets")
    for key, value in presets.items():
        if not isinstance(value, dict) or not isinstance(value.get("slots"), list):
            raise RuntimeError(f"{path}: preset '{key}' has no slot list")
        value.setdefault("label", key)
        value.setdefault("description", "")
    if default not in presets:
        raise RuntimeError(f"{path}: default '{default}' is not one of its presets")
    return presets, default


PRESETS, DEFAULT_PRESET = load_builtin_presets()


def slot_label(slot: dict, subclass: str | None = None) -> str:
    """Display name for one slot, resolved against a subclass when known."""
    kind = slot["kind"]
    if kind == "part":
        return slot["part"]
    if kind == "vl":
        return "VL"
    if kind == "vh":
        return "VH"
    if kind == "light_constant":
        resolved = CL_FOR_SUBCLASS.get(subclass or "")
        return resolved or " / ".join(CL_FOR_SUBCLASS.values())
    if kind == "stop":
        return "Stop"
    return kind


def validate_slots(slots) -> list[str]:
    """Check a hand-built layout's structure. Problems, empty if sound.

    This checks only the shape of the layout, never whether the parts it names
    are currently loaded: a layout can legitimately be designed before its
    parts file exists. Availability is reported separately by `missing_parts`,
    which alerts in the UI and blocks the build, rather than preventing the
    layout from being saved.
    """
    if not isinstance(slots, list) or not slots:
        return ["a layout needs at least one component"]

    problems: list[str] = []
    for position, slot in enumerate(slots, 1):
        if not isinstance(slot, dict):
            problems.append(f"component {position}: not a slot")
            continue
        kind = slot.get("kind")
        if kind not in SLOT_KINDS:
            problems.append(f"component {position}: unknown kind {kind!r}")
            continue
        if kind == "part" and not slot.get("part"):
            problems.append(f"component {position}: no part named")

    kinds = [s.get("kind") for s in slots if isinstance(s, dict)]
    # The stop codon terminates translation, so nothing can follow it.
    if kinds.count("stop") > 1:
        problems.append("a layout can carry at most one stop codon")
    elif "stop" in kinds and kinds[-1] != "stop":
        problems.append("the stop codon must be the last component")
    for kind, label in (("vl", "VL"), ("vh", "VH")):
        if kinds.count(kind) > 1:
            problems.append(f"a layout can carry at most one {label}")
    if kinds.count("light_constant") > 1:
        problems.append("a layout can carry at most one light constant region")
    return problems


def describe_layout(
    preset: str = DEFAULT_PRESET,
    parts: dict[str, str] | None = None,
    subclass: str | None = None,
    slots=None,
) -> list[dict]:
    """Describe a layout's slots for display, without assembling anything.

    With `subclass` given, the light constant slot resolves to the actual part
    it will take; without it, the slot shows both candidates. This is what lets
    the layout be shown before an antibody has been chosen.

    Pass `slots` to describe a layout that is not a stored preset, such as one
    being edited in the preset builder.
    """
    parts = parts or {}
    if slots is None:
        slots = PRESETS[preset]["slots"]
    described = []
    for slot in slots:
        kind = slot["kind"]
        label = slot_label(slot, subclass)
        entry = {
            "kind": kind,
            "label": label,
            "length": None,
            "detail": "",
            "resolved": None,
        }
        if kind == "part":
            seq = parts.get(slot["part"])
            entry["resolved"] = slot["part"]
            if seq is not None:
                entry["length"] = len(seq.rstrip("*"))
                entry["detail"] = f"{entry['length']} aa"
            else:
                entry["detail"] = "not in the parts table"
        elif kind in ("vl", "vh"):
            entry["detail"] = "per antibody"
        elif kind == "light_constant":
            resolved = CL_FOR_SUBCLASS.get(subclass or "")
            entry["resolved"] = resolved
            if resolved:
                seq = parts.get(resolved)
                entry["length"] = len(seq.rstrip("*")) if seq else None
                entry["detail"] = f"{subclass}"
            else:
                entry["detail"] = "chosen per antibody"
        elif kind == "stop":
            entry["detail"] = "appended on build"
        described.append(entry)
    return described


@dataclass
class Segment:
    """One labelled stretch of the finished construct."""

    label: str
    seq: str
    start: int  # 1-based residue number in the construct
    end: int  # inclusive
    kind: str = ""  # the slot kind that produced it
    source: str = ""  # part name, or the antibody the V domain came from


@dataclass
class Construct:
    name: str
    seq: str = ""
    subclass: str | None = None
    segments: list[Segment] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors and bool(self.seq)

    @property
    def has_stop(self) -> bool:
        """Whether the construct ends in a stop codon, as its parts dictated."""
        return self.seq.endswith("*")

    @property
    def residue_count(self) -> int:
        """Translated length, not counting a trailing stop codon."""
        return len(self.seq.rstrip("*"))

    @property
    def size_note(self) -> str:
        return f"{self.residue_count} aa" + (" + stop" if self.has_stop else "")

    @property
    def layout(self) -> str:
        return " | ".join(f"{s.label}({len(s.seq)})" for s in self.segments)


def load_parts(path) -> tuple[dict[str, str], list[str]]:
    """Read misc.fasta and validate the parts against the knowledge base.

    Returns (parts, problems). A non-empty `problems` list means the part file
    cannot be trusted and nothing should be assembled from it.
    """
    from .fasta import read_fasta

    raw = {r.name: r.seq for r in read_fasta(path)}
    problems: list[str] = []
    parts: dict[str, str] = {}

    for name, seq in raw.items():
        # A trailing '*' in a source file is tolerated and dropped: parts are
        # stored as residues only, and the stop codon is contributed by the
        # layout instead (the `stop` slot), so it lands exactly once at the
        # C-terminus of the finished construct no matter which parts are used.
        core = seq[:-1] if seq.endswith("*") else seq
        if "*" in core:
            problems.append(
                f"'{name}' has a '*' before its C-terminus, which means the "
                f"sequence is truncated"
            )
        bad = sorted(set(core) - rules.AMINO_ACIDS)
        if bad:
            problems.append(f"'{name}' contains non-amino-acid characters: {''.join(bad)}")
        parts[name] = core

    for name in REQUIRED_PARTS:
        if name not in parts:
            problems.append(f"misc file is missing a '{name}' sequence")
        elif not parts[name]:
            problems.append(f"'{name}' is empty")

    for name, expected in EXPECTED_START.items():
        seq = parts.get(name)
        if seq and not seq.startswith(expected):
            problems.append(
                f"'{name}' should begin at its own first residue {expected} "
                f"but begins {seq[: len(expected)]}"
            )

    return parts, problems


def pair_chains(results: list[ChainResult]) -> dict[str, dict[str, ChainResult]]:
    """Group validated chain results by antibody name.

    `GL_VRC01_VH` and `GL_VRC01_VL` pair under `GL_VRC01`.
    """
    paired: dict[str, dict[str, ChainResult]] = {}
    for result in results:
        base, _ = parse_declared(result.name)
        slot = result.detected or "unknown"
        paired.setdefault(base, {})[slot] = result
    return paired


def build_from_layout(
    antibody_id: str,
    vh_seq: str,
    vl_seq: str,
    subclass: str | None,
    parts: dict[str, str],
    preset: str = DEFAULT_PRESET,
    slots=None,
) -> Construct:
    """Assemble one construct from plain sequences and a layout.

    This is the form the stored CSV tables feed: it takes sequences and a
    kappa/lambda call rather than ChainResult objects, so a construct can be
    rebuilt from the library without re-reading any FASTA.

    `slots` takes precedence over `preset`, which lets a saved user layout be
    built without registering it in PRESETS.
    """
    construct = Construct(name=antibody_id)

    if slots is None:
        if preset not in PRESETS:
            construct.errors.append(f"unknown layout preset '{preset}'")
            return construct
        slots = PRESETS[preset]["slots"]

    kinds = {s["kind"] for s in slots if isinstance(s, dict)}
    # Only demand a chain the layout actually places.
    if "vh" in kinds and not vh_seq:
        construct.errors.append("no heavy chain (VH) for this antibody")
    if "vl" in kinds and not vl_seq:
        construct.errors.append("no light chain (VL) for this antibody")

    cl_name = CL_FOR_SUBCLASS.get(subclass or "")
    needs_cl = any(s["kind"] == "light_constant" for s in slots)
    if needs_cl and cl_name is None:
        construct.errors.append(
            f"light chain is neither kappa nor lambda (subclass "
            f"{subclass or 'unknown'!r}), so no constant region can be chosen"
        )

    for slot in slots:
        if slot["kind"] == "part" and slot["part"] not in parts:
            hint = suggest_part(slot["part"], parts)
            construct.errors.append(
                f"the layout needs a part named '{slot['part']}' but the parts "
                f"table has no such entry"
                + (f" — did you mean '{hint}'?" if hint else "")
            )
    if needs_cl and cl_name and cl_name not in parts:
        construct.errors.append(
            f"'{cl_name}' is not in the parts table, so a {subclass} light "
            f"chain cannot be completed"
        )

    if construct.errors:
        return construct

    construct.subclass = subclass

    # The subclass call rests on the V region alone: a trimmed V domain has no
    # constant region left to cross-check against, which is exactly the
    # two-checkpoint rule section 3 asks for. The terminal K/L of the V domain
    # has no exceptions among functional human J genes (section 3b), so it is
    # a sound basis for picking the CL -- but the call is reported so it can be
    # checked against the original, untrimmed sequence.
    if cl_name:
        ends = vl_seq[-1] if vl_seq else "?"
        construct.warnings.append(
            f"{cl_name} chosen from the V-region kappa/lambda call "
            f"(VL ends {ends} -> {subclass}); no CL was present in the "
            f"trimmed V domain to cross-check"
        )

    pieces: list[str] = []
    cursor = 1
    for slot in slots:
        kind = slot["kind"]
        if kind == "part":
            seq, label, source = parts[slot["part"]], slot["part"], slot["part"]
        elif kind == "vl":
            seq, label, source = vl_seq, f"VL:{antibody_id}", antibody_id
        elif kind == "vh":
            seq, label, source = vh_seq, f"VH:{antibody_id}", antibody_id
        elif kind == "light_constant":
            seq, label, source = parts[cl_name], cl_name, cl_name
        elif kind == "stop":
            seq, label, source = STOP, "Stop", "layout"
        else:
            construct.errors.append(f"unknown slot kind '{kind}'")
            return construct

        end = cursor + len(seq) - 1
        construct.segments.append(
            Segment(label=label, seq=seq, start=cursor, end=end, kind=kind, source=source)
        )
        pieces.append(seq)
        cursor = end + 1

    construct.seq = "".join(pieces)
    _check_junctions(construct)
    return construct


def build(
    name: str,
    vh: ChainResult | None,
    vl: ChainResult | None,
    parts: dict[str, str],
    preset: str = DEFAULT_PRESET,
    slots=None,
) -> Construct:
    """Assemble from two freshly typed chains, as the `build` CLI does."""
    construct = Construct(name=name)
    if vh is None:
        construct.errors.append("no heavy chain (VH) for this antibody")
    if vl is None:
        construct.errors.append("no light chain (VL) for this antibody")
    if construct.errors:
        return construct

    construct = build_from_layout(
        name, vh.clean_seq, vl.clean_seq, vl.subclass, parts, preset, slots
    )

    # Carry forward anything the typing step was unsure about, so a weak
    # boundary call stays visible on the finished construct.
    for source in (vl, vh):
        if source.confidence != "high":
            construct.warnings.append(
                f"{source.name} was a {source.confidence}-confidence boundary "
                f"call; its junctions in this construct inherit that doubt"
            )
    return construct


def _check_junctions(construct: Construct) -> None:
    """Sanity-check the finished construct against the knowledge base."""
    seq = construct.seq
    labels = {s.label.split(":")[0]: s for s in construct.segments}

    # Section 6: no linker between the signal peptide and the mature
    # N-terminus. The VL must start immediately after the signal.
    signal = labels.get("Signal")
    vl_seg = next((s for s in construct.segments if s.label.startswith("VL:")), None)
    if signal and vl_seg and vl_seg.start != signal.end + 1:
        construct.errors.append("signal peptide is not flush against the VL N-terminus")
    if signal and vl_seg and seq[signal.end : signal.end + 2] == "GS":
        construct.warnings.append(
            "mature N-terminus begins GS, which abolishes N-terminal "
            "pyroglutamate formation (section 6) -- check the leader junction"
        )

    # The V/C junctions should read as the knowledge base says they do.
    for v_label, c_label in (("VL:", None), ("VH:", "CHandCT")):
        seg = next((s for s in construct.segments if s.label.startswith(v_label)), None)
        if seg is None:
            continue
        last = seg.seq[-1]
        if v_label == "VH:" and not seg.seq.endswith(("VTVSS", "VTVSA")):
            construct.warnings.append(
                f"VH does not end VTVSS/VTVSA (ends {seg.seq[-5:]})"
            )
        if v_label == "VL:" and last not in "KL":
            construct.warnings.append(f"VL does not end K or L (ends {last})")

    # A single-chain product translates straight through to its C-terminus, so
    # a stop is correct only as the final character. One in the middle means a
    # part was truncated and everything after it would be untranslated.
    if "*" in seq[:-1]:
        position = seq.index("*") + 1
        owner = next(
            (s.label for s in construct.segments if s.start <= position <= s.end),
            "unknown segment",
        )
        construct.errors.append(
            f"stop codon at residue {position}, inside {owner}; everything "
            f"downstream of it would not be translated"
        )

# The shipped layouts are checked at import, so a bad edit to
# builtin_presets.json is caught here rather than at build time.
for _key, _value in PRESETS.items():
    _problems = validate_slots(_value["slots"])
    if _problems:
        raise RuntimeError(
            f"built-in preset '{_key}' in {BUILTIN_PRESET_FILE.name} is not a "
            f"valid layout: {'; '.join(_problems)}"
        )
del _key, _value, _problems

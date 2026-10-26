"""Command line entry point.

    python3 -m bcr.cli clean sequences.fasta
    python3 -m bcr.cli build sequences_clean.fasta misc.fasta

`clean` checks each sequence really is the chain its header claims
(`<antibody>_VH` / `<antibody>_VL`), trims away any signal peptide and any
constant region, and writes the V domains to sequences_clean.fasta.

`build` assembles a membrane-bound BCR construct per antibody:

    Signal | VL | IGKC or IGLC2 | Linker | VH | CHandCT

It re-runs the `clean` checks on its input, so it accepts either a raw or an
already-trimmed FASTA.

Exit status is 1 if any record failed, so either command can gate a pipeline.
"""

import argparse
import sys

from . import assemble, chains
from .fasta import Record, read_fasta, write_fasta

DEFAULT_CLEAN_OUTPUT = "sequences_clean.fasta"
DEFAULT_BUILD_OUTPUT = "bcr_constructs.fasta"


def _report_chain(result: chains.ChainResult, verbose: bool) -> None:
    if result.ok:
        label = result.detected
        if result.subclass:
            label = f"{label}/{result.subclass}"
        trimmed = []
        if result.n_trimmed:
            trimmed.append(f"-{result.n_trimmed} N")
        if result.c_trimmed:
            trimmed.append(f"-{result.c_trimmed} C")
        trim_note = f"  [{', '.join(trimmed)}]" if trimmed else "  [already clean]"
        print(
            f"  OK    {result.name}: {label}, {len(result.clean_seq)} aa"
            f"{trim_note}  ({result.confidence} confidence)"
        )
    else:
        print(f"  FAIL  {result.name}:")
        for err in result.errors:
            print(f"          error: {err}")

    for warn in result.warnings:
        print(f"          warning: {warn}")
    if verbose and result.evidence:
        for item in result.evidence:
            print(f"          evidence: {item}")


def _load_and_type(path: str, verbose: bool):
    """Read a FASTA and run chain typing on every record."""
    try:
        records = read_fasta(path)
    except OSError as exc:
        print(f"error: cannot read {path}: {exc}", file=sys.stderr)
        return None

    if not records:
        print(f"error: no sequences found in {path}", file=sys.stderr)
        return None

    print(f"Read {len(records)} sequence(s) from {path}\n")
    results = [chains.process(rec.name, rec.seq) for rec in records]
    for result in results:
        _report_chain(result, verbose)
    return results


def cmd_clean(args) -> int:
    results = _load_and_type(args.input, args.verbose)
    if results is None:
        return 2

    passed = [r for r in results if r.ok]
    failed = [r for r in results if not r.ok]
    print(f"\n{len(passed)} passed, {len(failed)} failed")

    if failed and args.strict:
        print("--strict: nothing written because a record failed", file=sys.stderr)
        return 1

    if passed:
        out = [
            Record(
                name=r.name,
                seq=r.clean_seq,
                description=(
                    f"{r.detected}"
                    + (f"/{r.subclass}" if r.subclass else "")
                    + f" V domain, {len(r.clean_seq)} aa"
                ),
            )
            for r in passed
        ]
        write_fasta(args.output, out)
        print(f"Wrote {len(out)} clean sequence(s) to {args.output}")
    else:
        print("No clean sequences to write.")

    return 1 if failed else 0


def cmd_build(args) -> int:
    parts, problems = assemble.load_parts(args.misc)
    if problems:
        print(f"error: {args.misc} cannot be used:", file=sys.stderr)
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
        return 2
    print(
        f"Loaded {len(parts)} part(s) from {args.misc}: "
        + ", ".join(f"{k}({len(v)})" for k, v in parts.items())
        + "\n"
    )

    results = _load_and_type(args.input, args.verbose)
    if results is None:
        return 2

    usable = [r for r in results if r.ok]
    rejected = [r for r in results if not r.ok]
    if rejected:
        print(
            f"\n{len(rejected)} sequence(s) failed typing and are excluded "
            "from assembly"
        )

    layout = " | ".join(
        assemble.slot_label(s)
        for s in assemble.PRESETS[assemble.DEFAULT_PRESET]["slots"]
    )
    print(f"\nAssembling constructs: {layout}\n")
    constructs = []
    for name, chain_map in sorted(assemble.pair_chains(usable).items()):
        construct = assemble.build(
            name, chain_map.get("heavy"), chain_map.get("light"), parts
        )
        constructs.append(construct)
        # The default layout is the only one the CLI builds, but the name
        # still says so: the web UI can build several, into one file.
        record_name = assemble.construct_name(name, assemble.DEFAULT_PRESET)

        if construct.ok:
            print(
                f"  OK    {record_name}: {construct.size_note} "
                f"({construct.subclass})"
            )
            print(f"          layout: {construct.layout}")
            if args.verbose:
                for seg in construct.segments:
                    print(f"          {seg.start:>5}-{seg.end:<5} {seg.label}")
        else:
            print(f"  FAIL  {name}:")
            for err in construct.errors:
                print(f"          error: {err}")
        for warn in construct.warnings:
            print(f"          warning: {warn}")

    built = [c for c in constructs if c.ok]
    failed = [c for c in constructs if not c.ok]
    print(f"\n{len(built)} construct(s) built, {len(failed)} failed")

    if failed and args.strict:
        print("--strict: nothing written because a construct failed", file=sys.stderr)
        return 1

    if built:
        out = [
            Record(
                name=assemble.construct_name(c.name, assemble.DEFAULT_PRESET),
                seq=c.seq,
                description=f"{c.subclass} BCR, {c.size_note} | {c.layout}",
            )
            for c in built
        ]
        write_fasta(args.output, out)
        print(f"Wrote {len(out)} construct(s) to {args.output}")
    else:
        print("No constructs to write.")

    return 1 if (failed or rejected) else 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="bcr", description="Human antibody chain typing and BCR assembly."
    )
    sub = parser.add_subparsers(dest="command", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-v", "--verbose", action="store_true", help="show details")
    common.add_argument(
        "--strict", action="store_true", help="write nothing unless every record passes"
    )

    p_clean = sub.add_parser(
        "clean",
        parents=[common],
        help="validate _VH/_VL labels and trim to the V domain",
    )
    p_clean.add_argument("input", help="input FASTA")
    p_clean.add_argument("-o", "--output", default=DEFAULT_CLEAN_OUTPUT)
    p_clean.set_defaults(func=cmd_clean)

    p_build = sub.add_parser(
        "build", parents=[common], help="assemble full BCR constructs"
    )
    p_build.add_argument("input", help="FASTA of V regions (raw or trimmed)")
    p_build.add_argument("misc", help="FASTA of parts: Signal, Linker, CHandCT, IGKC, IGLC2")
    p_build.add_argument("-o", "--output", default=DEFAULT_BUILD_OUTPUT)
    p_build.set_defaults(func=cmd_build)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())

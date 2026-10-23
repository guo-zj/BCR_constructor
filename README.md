# BCR Constructor

Assembles B-cell-receptor amino acid sequences from antibody V regions.

Give it a FASTA of VH and VL sequences. It confirms each one really is the
chain its header claims, trims away any signal peptide or constant region to
leave the bare V domain, and then builds full constructs by splicing those
domains into a layout of parts — picking the right light constant region for
each antibody from its own κ/λ call.

> **Human antibodies only.** The typing and boundary rules are transcribed from
> [`antibody-chain-rules.md`](antibody-chain-rules.md), which covers human
> sequences; murine and other species differ at several of them.
>
> These rules are for fast automated calls and sanity-checking. For anything
> production-relevant, confirm with [ANARCI](https://github.com/oxpig/ANARCI)
> or [IMGT/DomainGapAlign](https://www.imgt.org/3Dstructure-DB/cgi/DomainGapAlign.cgi).
> A single wrong residue in a constant region survives into the plasmid.

## Requirements

- Python 3.10 or newer (developed on 3.14)
- [Flask](https://flask.palletsprojects.com/) 3.x — only for the web UI

Everything else is standard library. The command line needs no third-party
packages at all; there is no build step and no JavaScript bundling.

```bash
pip install flask        # only if you want the web UI
```

## Quick start

### Web UI

```bash
python3 -m bcr.web            # then open http://127.0.0.1:8765
```

The page has four cards for the main flow, and a **Utilities** banner along
the bottom that opens over the page when clicked:

| Where | What it does |
|---|---|
| Utilities → Sequence sanity check | Type and trim a VH/VL FASTA, review, then approve |
| Utilities → Construct parts | Import the parts FASTA, resolving any conflicts |
| Utilities → Assembly presets | Build layouts by dragging components into order |
| 1 · Antibody library | Stored V domains; tick the ones to assemble |
| 2 · Assembly layout | Pick a preset and see the order it will assemble in |
| 3 · Build plan | Which constant region each antibody gets; confirm to write |
| 4 · Constructs library | What has been built, with its segment map |

A typical first run: open Utilities, import your parts file, run the sequence
check and approve the antibodies, close the panel, tick the antibodies in card
1, review the plan in card 3, and build.

Nothing is stored until you approve it, and the two import flows always show
you what they are about to do first. The server binds `127.0.0.1` only.

### Command line

```bash
# Type and trim; writes sequences_clean.fasta
python3 -m bcr.cli clean sequences.fasta

# Assemble full constructs; writes bcr_constructs.fasta
python3 -m bcr.cli build sequences.fasta misc_2A.fasta
```

```
  OK    GL_VRC01_VH: heavy, 121 aa  [already clean]  (high confidence)
  OK    GL_VRC01_VL: light/kappa, 103 aa  [already clean]  (high confidence)
  OK    Test_VH: heavy, 116 aa  [already clean]  (high confidence)
  OK    Test_VL: light/kappa, 111 aa  [-1 C]  (high confidence)

Assembling constructs: SignalHC | VL | IGKC / IGLC2 | LinkerSC | VH | CHandCTSC | Stop

  OK    GL_VRC01_BCR: 806 aa + stop (kappa)
          layout: SignalHC(19) | VL:GL_VRC01(103) | IGKC(107) | LinkerSC(57) | VH:GL_VRC01(121) | CHandCTSC(399) | Stop(1)
  OK    Test_BCR: 809 aa + stop (kappa)
```

Both commands exit non-zero if any record failed, so they can gate a pipeline.

## Input files

### Sequences

A FASTA whose headers are `<antibody>_VH` or `<antibody>_VL`. The suffix is
checked against what the sequence actually reads as, and a contradiction is a
hard failure — nothing is written for that record.

```
>GL_VRC01_VH
QVQLVQSGAEVKKPGASVKVSCKAS...WGQGTLVTVSS
>GL_VRC01_VL
EIVLTQSPATLSLSPGERATLSC...FGQGTKLEIK
```

Sequences may carry a signal peptide, run on into the constant region, or
both — those are trimmed off. A record with no `_VH`/`_VL` suffix is still
typed, with a warning. An antibody is stored only when **both** a clean VH and
a clean VL are present.

### Parts

A FASTA of the pieces a construct is spliced from. **Names are yours to
choose** — they come from the file, not from the code — with two exceptions:
the light constant regions must be called `IGKC` and `IGLC2`, because a layout
cannot name one outright. Which of the two is used is decided per antibody
from its κ/λ call.

```
>SignalHC
MGWSCIILFLVATATGVHS
>LinkerSC
GGSSGSGSGSNWSHPQFEK...
>IGKC
RTVAAPSVFIFPPSDEQLK...
>IGLC2
GQPKAAPSVTLFPPSSEELQ...
```

Each part is stored as residues only: a trailing `*` is dropped on import and
one stop codon is appended at the end of the finished construct instead. A
constant region must begin at its own first residue (`RTVAAPSV`,
`GQPKAAPSV`), and a `*` anywhere but the end is reported as a truncation.

Re-importing a parts file compares each part against what is stored and asks
what to do about any that differ. The default is to keep what you have.

## Layouts

A layout is an ordered list of components, N- to C-terminus. The one that
ships is **Single chain BCR**:

```
SignalHC → VL → IGKC/IGLC2 → LinkerSC → VH → CHandCTSC → Stop
```

Components are of five kinds: a named `part` from the parts table, the
antibody's own `vl` or `vh`, a `light_constant` that resolves per antibody,
and the `stop` codon. Build your own in Utilities → Assembly presets by
dragging components into order; they are saved to `data/presets.json` and
become selectable alongside the built-in one.

If a layout names a part the table does not hold, you get an alert naming it —
with a "did you mean" suggestion if something close exists — and building is
blocked rather than producing a construct with a gap in it. A layout can be
designed before its parts exist; it just cannot be built yet.

Built-in layouts live in `bcr/builtin_presets.json`. Editing that file changes
what ships, with no code change.

## Stored data

| File | Holds |
|---|---|
| `data/antibodies.csv` | trimmed VH/VL, with chain type, κ/λ, confidence and how much was trimmed |
| `data/parts.csv` | construct parts, with the role inferred from each sequence |
| `data/constructs.csv` | assembled constructs and their segment maps |
| `data/presets.json` | layouts you have built |

Plain CSV, so they open in a spreadsheet or feed any other tool. Writes are
atomic — rendered to a temporary file and moved into place — so an interrupted
write cannot truncate a table. Re-processing an antibody replaces its row
rather than adding a duplicate.

## How the calls are made

**Heavy or light** is decided by scoring both hypotheses and taking the
stronger, rather than checking one signal at a time: the framework 4 motif
(`WGxG` with Trp for heavy, `FGxG` with Phe for light), the N-terminal
`SG`/`TQ` motif, the conserved FR1 Cys, identity to germline FR1, the constant
region opening if present, and the CDR3 lead-in.

**Where to trim.** The C-terminal boundary is exact: the VH domain ends at
Trp + 10, the VL domain at Phe + 9. The N-terminal boundary has no equivalent
counting rule, so every start position that yields a plausible V length is
scored against four independent FR1 features and the best is taken. That score
is reported as high, medium or low confidence, so a weak call is visible rather
than silent.

**κ or λ** uses one V-region and one C-region checkpoint where both are
available, because a Vλ-on-Cκ hybrid is a real engineered format that a
V-region-only call would mislabel. For a stored V domain there is no constant
region left to check against, so the call rests on the terminal K or L alone
and every construct carries a warning saying so.

## Tests

```bash
python3 -m unittest discover -s tests -t .
```

154 tests, no third-party test runner. They cover the counting rules against
all functional human J segments, the κ/λ hybrid trap, boundary trimming,
construct assembly, the CSV and JSON stores, and every web endpoint. Test
inputs are pinned in `tests/fixtures/`, so editing your own working files will
not break the suite.

## Limitations

- The λ path is exercised only by a synthetic fixture built from germline
  IGLV1 frameworks; no real λ sequence has been through it.
- N-terminal boundary scoring is calibrated on the germline FR1s listed in the
  knowledge base — three per chain, with only one IGLV among them. λ
  N-terminal calls are the least-tested part of typing.
- Building replaces the output FASTA rather than appending to it. The UI says
  so and reports how many records were displaced.
- Nothing here substitutes for ANARCI or IMGT. Read well past any boundary
  before trusting it.

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
| Utilities → Sequence sanity check | Drop or point at a VH/VL FASTA; type and trim it, review, then approve |
| Utilities → Construct parts | Drop or point at the parts FASTA, resolving any conflicts |
| Utilities → Assembly presets | Build layouts by dragging components into order |
| 1 · Antibody library | Stored V domains; tick the ones to assemble, or edit them by hand |
| 2 · Assembly layout | Pick a preset and see the order it will assemble in |
| 2 · Express mode *(turn the card over)* | Each antibody's whole light and heavy chains, and the four express layouts |
| 3 · Build plan | Which constant region each antibody gets; confirm to write |
| 4 · Constructs library | What has been built, with its segment map; rename and annotate |

Both input tabs accept a **`.fa` / `.fasta` file dragged onto them** — the
file is copied to `data/uploads/` and processed from there — or a path to a
file already on the machine. The sanity-check tab shows the header naming
rule as a small diagram, if you would rather look than read.

A typical first run: open Utilities, import your parts file, run the sequence
check and approve the antibodies, close the panel, tick the antibodies in card
1, review the plan in card 3, and build.

Nothing is stored until you approve it, and the two import flows always show
you what they are about to do first. **Trimming is offered, not applied:** any
record the rules would shorten gets a trim-or-keep choice of its own, and
approving leaves it alone until you answer. The server binds `127.0.0.1` only.

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
  OK    C05_VH: heavy, 138 aa  [already clean]  (high confidence)
          warning: V domain is 138 aa, outside the usual 98-128 for heavy

Assembling constructs: SignalHC | Light chain (VL+CL) | LinkerSC | Heavy chain IgG1 | Stop

  OK    GL_VRC01_BCR: 806 aa + stop (kappa)
          layout: SignalHC(19) | LC:GL_VRC01(VL+IGKC)(210) | LinkerSC(57) | HC:GL_VRC01(VH+CHTMCTIgG1)(520) | Stop(1)
  OK    C05_BCR: 827 aa + stop (kappa)
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
both. Those can be trimmed off — but you are asked first, per record, and
"keep as given" stores the sequence exactly as it came, with the row saying
so. A record with no `_VH`/`_VL` suffix is still
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

## Editing the library

Every table can be edited in place: **Edit** on a row, or **+ Add by hand** to
type a record in from scratch.

A hand-typed antibody is held to exactly the same rules as an imported one.
Paste in a sequence — or a whole FASTA record, header and line breaks
included — and on save it is typed, trimmed to its V domain, and κ/λ called,
with the result reported back. A sequence that does not type as a VH/VL pair
is refused with the reasons, rather than stored as a row nothing else can
reason about. Parts go through the same checks the parts file does: a trailing
`*` is dropped, an internal one is a truncation, and a light constant must
still begin at its own first residue.

Trimming a pasted sequence is likewise offered rather than done: if the rules
would shorten it, nothing is saved until you pick trim or keep, with the
residues in question shown.

Renames follow through. Renaming an antibody repoints the constructs built
from it; renaming a part rewrites the slots of every layout you have saved,
and tells you if a **built-in** layout still asks for the old name. `IGKC` and
`IGLC2` cannot be renamed or deleted, because a `light_constant` slot finds
them by name — but a corrected sequence can be saved over them, which is
usually what was wanted.

A construct is assembled rather than typed in, so only its id and your notes
are editable; the sequence and segment map it was built with are shown
read-only. Rebuild to change them.

## Express mode

Card 2 turns over. The back shows, for every antibody in the library, the three
whole chains it implies:

| Column | What it is |
|---|---|
| `light_chain` | VL + `IGKC` or `IGLC2`, chosen from that antibody's own κ/λ call |
| `heavy_chain_igg1` | VH + `CHTMCTIgG1` |
| `heavy_chain_igm` | VH + `CHTMCTIgM` |

They are kept in `antibodies.csv` as well, so the table is useful on its own.
They are derived, not authored: the same function that fills them resolves a
whole-chain slot at build time, so what you see is what gets built. If a row's
columns are missing or disagree with its own V domains, the panel says so and
offers to write them — it will not rewrite the CSV behind your back.

Four layouts come with express mode:

```
Single chain BCR_IGG1   SignalHC → light chain → LinkerSC → heavy chain IgG1 → Stop
Single chain BCR_IGM    SignalHC → light chain → LinkerSC → heavy chain IgM  → Stop
2Alinked_IGG1           SignalHC → heavy chain IgG1 → Furin → F2A → SignalLC
                                 → light chain → T2A → EGFP → Stop
2Alinked_IGM            as above, with the IgM heavy chain
```

A whole-chain component is shorthand, not a different construct: it assembles
exactly what naming VL, the light constant, VH and the CH part by hand would,
and a test holds the two byte-identical. Both are draggable in the preset
builder, so your own layouts can use either form.

## Intrinsic parts

Four parts ship with the package and cannot be edited, renamed or deleted:

| Name | Length | Source |
|---|---|---|
| `IGKC` | 107 aa | §5 of the knowledge base |
| `IGLC2` | 106 aa | §5 of the knowledge base |
| `CHTMCTIgG1` | 399 aa | IgG1 CH1–CH3 + transmembrane and cytoplasmic tail |
| `CHTMCTIgM` | 474 aa | IgM CH1–CH4 + transmembrane and cytoplasmic tail |

Every layout resolves to them — by name for the heavy constants, per antibody
for the light ones — so they are held outside the parts table entirely rather
than stored as rows that merely refuse to change. Their four names are
**reserved**: if a parts file you import uses one, you are offered a different
name for that record instead, and the intrinsic sequence stays as it is. A
record whose sequence is identical to the intrinsic one needs no rename, since
it would override nothing.

## Layouts

A layout is an ordered list of components, N- to C-terminus. Four ship, all of
them express layouts; **Single chain BCR_IGG1** is the default:

```
SignalHC → light chain → LinkerSC → heavy chain IgG1 → Stop
```

Components are of seven kinds: a named `part` from the parts table, the
antibody's own `vl` or `vh`, a `light_constant` that resolves per antibody, a
whole `light_chain` or `heavy_chain` (the shorthand express mode uses), and the
`stop` codon. Build your own in Utilities → Assembly presets by
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
| `data/uploads/` | FASTA files you dropped onto the page, kept exactly as given |

`antibodies.csv` also holds the three whole-chain columns described under
[Express mode](#express-mode).

Each table has two note columns: `notes`, carrying whatever the rules warned
about, and `user_notes`, which is yours. Editing a record touches only the
second, so your own remarks and the rules' never overwrite each other.

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

236 tests, no third-party test runner. They cover the counting rules against
all functional human J segments, the κ/λ hybrid trap, boundary trimming,
construct assembly, the CSV and JSON stores, and every web endpoint. Test
inputs are pinned in `tests/fixtures/`, so editing your own working files will
not break the suite.

## Limitations

- The λ path is exercised only by a synthetic fixture built from germline
  IGLV1 frameworks; no real λ sequence has been through it.
- A VH is scored for start positions up to 150 residues, wider than the
  knowledge base's 98–128, because CDR-H3 length is effectively unbounded —
  C05, at 138 aa, was having 15 residues of FR1 cut off to fit a tighter
  window. A length outside 98–128 is still reported as unusual.
- N-terminal boundary scoring is calibrated on the germline FR1s listed in the
  knowledge base — three per chain, with only one IGLV among them. λ
  N-terminal calls are the least-tested part of typing.
- Building replaces the output FASTA rather than appending to it. The UI says
  so and reports how many records were displaced.
- Nothing here substitutes for ANARCI or IMGT. Read well past any boundary
  before trusting it.

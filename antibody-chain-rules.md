# Antibody Sequence Rules — Human

Quick reference for chain typing and V/C boundary assignment.
**Scope: human antibodies only.** Murine and other species differ at several of these
rules (noted where the difference is a common trap).

**Always confirm with [ANARCI](https://github.com/oxpig/ANARCI) or
[IMGT/DomainGapAlign](https://www.imgt.org/3Dstructure-DB/cgi/DomainGapAlign.cgi)
before anything production-relevant.** The rules below are for fast manual calls and
for sanity-checking automated output — not a substitute for it.

---

## 1. Heavy chain or light chain?

Three independent signals. Any one usually settles it; use two to be safe.

### 1a. Framework 4 motif — the strongest single marker

| Chain | FR4 motif | Conserved residue |
|---|---|---|
| **Heavy** | `WGxG` | **Trp** (Kabat 103 / IMGT 118) |
| **Light** | `FGxG` | **Phe** (Kabat 98 / IMGT 118) |

Trp → heavy. Phe → light. This holds across all functional human J genes.

### 1b. N-terminal motif

Look at residues ~4–8:

| Chain | Pattern | Examples |
|---|---|---|
| **Heavy** | contains **`SG`** | `QVQLVQ**SG**AEVKKPGASVKVSCKAS` (IGHV1)<br>`EVQLVE**SG**GGLVQPGGSLRLSCAAS` (IGHV3)<br>`QVQLQE**SG**PGLVKPSETLSLTCTVS` (IGHV4) |
| **Light** | contains **`TQ`** | `DIQM**TQ**SPSSLSASVGDRVTITC` (IGKV1)<br>`EIVL**TQ**SPGTLSLSPGERATLSC` (IGKV3)<br>`QSVL**TQ**PPSVSGAPGQRVTISC` (IGLV1) |

Mnemonic: **S-G = heavy, T-Q = light.**

> Germline FR1 as shown. Somatic hypermutation and degenerate cloning primers both
> alter the first few residues routinely, so treat N-terminal matching as supporting
> evidence, not proof.

### 1c. C-terminus of the V domain

| Chain | Ends |
|---|---|
| **Heavy** | `...VTVSS` |
| **Light, kappa** | `...EIK` / `...DIK` (terminal **K**) |
| **Light, lambda** | `...TVL` (terminal **L**) |

### 1d. Supporting checks

- **CDR3 lead-in:** heavy is usually `YYCAR` / `YYCAK` / `YYCTR`. Light is `YYC` /
  `YFC` / `FYC` followed by a non-AR residue.
- **Length:** VH ≈ 98–128 linear residues (varies most with CDR-H3). VL ≈ 107–115.
  Full IgG1 heavy ≈ 450; full light ≈ 214–216.
- **Constant region start:** IgG1 CH1 opens `ASTKGPSVFPLAP`. Nothing in a light
  chain resembles this.

---

## 2. Splicing the V region — heavy chain

### The rule

Find the conserved **Trp** of `WGxG` in FR4. **The VH domain ends at Trp + 10.**

```
... Y Y C A R <CDR-H3> W G Q G T L V T V S S | A S T K G P S V F P L A P ...
                                ^Trp                    ^+10  ^CH1 starts
```

Count: `W`(1) `G`(2) `x`(3) `G`(4) `T`(5) `x`(6) `V`(7) `T`(8) `V`(9) `S`(10) `S`(11)
— the domain's last residue is the second `S`, which is Trp+10.

### Human IGHJ framework 4 — all six functional genes

| Gene | FR4 |
|---|---|
| IGHJ1 | `WGQG`**`TLVTVSS`** |
| IGHJ2 | `WGRG`**`TLVTVSS`** |
| IGHJ3 | `WGQG`**`TMVTVSS`** |
| IGHJ4 | `WGQG`**`TLVTVSS`** |
| IGHJ5 | `WGQG`**`TLVTVSS`** |
| IGHJ6 | `WGQG`**`TTVTVSS`** |

Variation sits at the `TLV` / `TMV` / `TTV` position and in `WGQG` vs `WGRG`.
**The terminal `VTVSS` is invariant across human IGHJ.**

### The junction

```
...WGQGTLVTVSS | ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSW...
      VH end        IgG1 CH1 begins
```

### Notes

- **Isotype is irrelevant to this boundary.** The VH terminus is set by the IGHJ gene.
  IgG1, IgG4, IgM, IgA — same VH, same `...VTVSS`, then whichever constant region follows.
- **Search `VTV[SA]`, not `TVSS`,** if your pipeline may see non-human sequences.
  Murine IGHJ3 ends `WGQGTLVTVS`**`A`** and will be silently dropped by a literal
  `TVSS` match.
- **FR4 can be somatically mutated.** Rare — FR4 is well conserved and mutation there
  is selected against — but it exists in matured repertoires.
- **Engineered formats modify it.** scFvs splice `...VTVSS` directly into `(G4S)3`;
  some library formats truncate FR4.
- **Camelid VHH** also generally ends `...VTVSS`, so the marker carries into nanobody work.
- **IgG1 allotypes** (G1m17,1 / G1m3,1 / G1m3 / G1m17) differ at a few positions
  including EU 214, 356, 358, 431. They do not affect the VH/CH1 boundary, but will
  show up as mismatches when diffing against a reference. Pull the exact allotype
  from IMGT/GENE-DB rather than from memory.

---

## 3. Kappa or lambda?

Four checkpoints. **Use at least two — one from the V region and one from the C
region** (see the warning at the end of this section).

### 3a. N-terminal motif

| | Typical FR1 start |
|---|---|
| **Kappa** | `DIQMTQ...` `DIVMTQ...` `EIVLTQ...` `EIVMTQ...` |
| **Lambda** | `QSVLTQ...` `QSALTQ...` `SYELTQ...` `NFMLTQ...` `QTVVTQ...` `QPVLTQ...` `QLVLTQ...` |

Rough heuristic: kappa tends to open `D`/`E`, lambda tends to open `Q`/`S`/`N`.
Not absolute — use with the checks below.

### 3b. Last residue of the V domain — single-letter discriminator

| | V ends in | J segments |
|---|---|---|
| **Kappa** | **`K`** | `...KVEIK` `...KLEIK` `...KVDIK` `...RLEIK` |
| **Lambda** | **`L`** | `...KVTVL` `...KLTVL` `...QLTVL` |

**No exceptions among functional human J genes.**

### 3c. Start of the constant region

| | First residues of CL |
|---|---|
| **IGKC** | `RTVAAPSV`**`FIFPPS`**`DEQLK...` |
| **IGLC2** | `GQPK`**`AAPSV`**`TLFPPS`**`SEELQ...` |

⚠️ These share `AAPSV` and `FPPS` — which is exactly how a κ/λ mix-up hides in plain
sight. The discriminating stretch is what sits between:

```
kappa    R T V A A P S V   F I F P P S D E Q L K ...
lambda   G Q P K A A P S V T L F P P S S E E L Q ...
         ^^^^ lambda-only 4-residue extension
```

- kappa: `AAPSV` → `FIFPPS`
- lambda: `AAPSV` → `TLFPPS`
- lambda carries an N-terminal `GQPK` that kappa lacks

### 3d. C-terminus

| | Ends | Interchain Cys |
|---|---|---|
| **Kappa** | `...VTKSFNRGE`**`C`** | **C107** — the terminal residue |
| **Lambda** | `...VEKTVAPTE`**`C`**`S` | **C105** — penultimate, trailing Ser |

Both bond to Cys220 (EU) in the IgG1 upper hinge. The asymmetry matters for disulfide
constraints in modelling and for interpreting non-reduced LC-MS.

### ⚠️ Read 20+ residues into the constant region

A boundary call alone does **not** confirm both halves match. Hybrid constructs exist —
Vλ grafted onto Cκ is a real engineered format (common-light-chain and κλ-body
bispecific platforms use a shared Cκ to simplify purification). A chain can begin with
nine genuine lambda CL residues (`GQPKAAPSV`) and then switch to kappa from `FIFPPS`
onward.

If you only read the first few residues of CL, you will call it lambda and be wrong.

### Downstream consequences of the call

| | Kappa | Lambda |
|---|---|---|
| CL length | 107 aa | 106 aa |
| Affinity capture | KappaSelect | LambdaFabSelect |
| Protein L | binds — but depends on the **V region** (VκI/III/IV framework), **not** CL | no binding |
| pI / charge profile | — | differs; shifts IEX/cIEF and HIC retention |
| Developability | — | somewhat higher aggregation/polyreactivity liability in large datasets |

Note the Protein L row: a **Vλ-on-Cκ** construct gains KappaSelect capture but still
will **not** bind Protein L, because that interaction reads the V region.

---

## 4. Splicing the V region — light chain

### The rule

Find the conserved **Phe** of `FGxG` in FR4. **The VL domain ends at Phe + 9.**
Identical rule for kappa and lambda.

```
kappa   ... Y Y C Q Q Y N S Y P L T  F G Q G T K V E I K | R T V A A P S V ...
                                      ^Phe          ^+9    ^CL starts

lambda  ... Y Y C A A W D D S L K G A V  F G G G T Q L T V L | G Q P K A A P S V ...
                                          ^Phe          ^+9    ^CL starts
```

Count: `F`(1) `G`(2) `x`(3) `G`(4) `T`(5) `x`(6) `V`(7) `x`(8) `x`(9) `x`(10) — the
domain's last residue is the tenth, i.e. Phe+9.

### Human light-chain J segments

| Kappa | FR4 | | Lambda | FR4 |
|---|---|---|---|---|
| IGKJ1 | `WT`**`FGQGTKVEIK`** | | IGLJ1 | `YV`**`FGTGTKVTVL`** |
| IGKJ2 | `YT`**`FGQGTKLEIK`** | | IGLJ2 | `VV`**`FGGGTKLTVL`** |
| IGKJ3 | `FT`**`FGPGTKVDIK`** | | IGLJ3 | `VV`**`FGGGTKLTVL`** |
| IGKJ4 | `LT`**`FGGGTKVEIK`** | | IGLJ6 | `NV`**`FGSGTKVTVL`** |
| IGKJ5 | `IT`**`FGQGTRLEIK`** | | IGLJ7 | `AV`**`FGGGTQLTVL`** |

Note IGKJ3 uses `FGPGT` rather than `FGQGT`/`FGGGT` — still matches `FGxG`.

### The junctions

```
kappa:   ...FGQGTKVEIK | RTVAAPSVFIFPPSDEQLKSGTASVVC...
lambda:  ...FGGGTQLTVL | GQPKAAPSVTLFPPSSEELQANKATLVC...
```

### Cross-checks

- **Conserved V Cys** at Kabat 88 / IMGT 104, in the `YYC` that opens CDR-L3,
  ~15–20 residues upstream of `FGxG`. Disulfides to Cys23. If you can locate `YYC`
  and `FGxG`, CDR-L3 is everything between them.
- **CL intradomain Cys:** IGKC position 27 (→80); IGLC position 28 (→81).
- **Total length:** a light chain that isn't ≈210–220 residues has something unusual
  going on.

---

## 5. Reference sequences

### IGKC — human kappa constant (107 aa)

```
RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC
```

### IGLC2 — human lambda constant (106 aa)

```
GQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS
```

IGLC2 is the standard choice for recombinant work — most highly expressed lambda
constant gene, and what most therapeutic lambda antibodies use. IGLC1, IGLC3 and
IGLC7 are also functional and differ at a handful of positions (IGLC7 most). Pull
those from IMGT/GENE-DB rather than from memory — a single wrong residue in a
constant region survives into the plasmid.

Kappa has Km allotypes (Km1 / Km1,2 / Km3) differing at IGKC positions 45 and 101.
The sequence above is the common Km3.

### IgG1 CH1 — first residues (for boundary identification)

```
ASTKGPSVFPLAP...
```

---

## 6. General cautions

**Linear position ≠ Kabat position.** CDR lengths shift everything downstream. A light
chain whose V domain ends at linear residue 110 is normal — don't assume the boundary
lands at Kabat 107 in raw coordinates. Under IMGT the V-DOMAIN ends at 128.

**The counting rules locate the junction; they do not validate the rest of the chain.**
Always read well past the boundary.

**Splice whole domains at defined boundaries.** Constant-region swapping is reliable —
it's the basis of chimeric and humanized antibodies. Swapping *part* of a domain is
what causes problems. Keep the J segment with the V domain; start the constant region
at its own first residue (`RTVAAPSV`, `GQPKAAPSV`, `ASTKGPSVFPLAP`).

**Signal peptides butt directly against the mature N-terminus — no GS linker.**
Signal peptidase reads the local sequence (small neutral residues at −1 and −3).
A linker destroys that context and, if cleavage still works, leaves `GS...` on the
mature protein — which abolishes N-terminal pyroglutamate formation and changes the
charge-variant profile. Validate the actual junction with SignalP 6.0; if cleavage
looks poor, swap the leader rather than padding it.

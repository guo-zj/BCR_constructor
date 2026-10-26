"""Chain typing and V-domain boundary assignment.

Two jobs, in order:

1. Decide whether a sequence really is a VH or a VL, and check that against
   the chain declared in its FASTA header (`<name>_VH` / `<name>_VL`).
2. Trim it to the V domain alone, dropping any signal peptide on the N side
   and any constant region on the C side.

The C-terminal boundary is exact: the knowledge base gives Trp + 10 for heavy
and Phe + 9 for light. The N-terminal boundary has no equivalent counting
rule, so it is found by scoring every start position that yields a plausible
V length against four independent FR1 features (first residue, the SG/TQ
N-terminal motif at its canonical offset, the conserved FR1 Cys, and identity
to germline FR1). The winning score is reported so a weak call is visible
rather than silent.
"""

from dataclasses import dataclass, field

from . import rules


@dataclass
class StartCall:
    """Where the mature V domain begins, and why we think so."""

    start: int
    score: float
    reasons: list[str] = field(default_factory=list)


@dataclass
class Hypothesis:
    """One chain-class hypothesis evaluated against a sequence."""

    chain: str  # "heavy" / "light"
    fr4_index: int  # index of the conserved Trp / Phe
    strict: bool  # matched the strict FR4 pattern
    v_start: int
    v_end: int  # inclusive
    score: float
    reasons: list[str] = field(default_factory=list)

    @property
    def length(self) -> int:
        return self.v_end - self.v_start + 1


@dataclass
class ChainResult:
    """Everything the CLI needs to report and write one record."""

    name: str
    declared: str | None  # "heavy" / "light" / None if header unlabelled
    detected: str | None
    subclass: str | None  # "kappa" / "lambda" / None
    clean_seq: str = ""
    # The sequence as given, minus a trailing stop codon: what gets stored if
    # the user declines the proposed trim and keeps the record as it came.
    input_seq: str = ""
    v_start: int = 0
    v_end: int = 0
    n_trimmed: int = 0
    c_trimmed: int = 0
    score: float = 0.0
    confidence: str = "none"
    evidence: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    # True when the rules rejected this record and the user chose to store it
    # anyway. The sequence is then whatever came in, untrimmed and untyped;
    # see accept_as_given.
    accepted_as_given: bool = False

    @property
    def ok(self) -> bool:
        return not self.errors and bool(self.clean_seq)


# --- header parsing -------------------------------------------------------


def parse_declared(name: str) -> tuple[str, str | None]:
    """Split `antibody_VH` into ("antibody", "heavy").

    Returns (base_name, None) when the header carries no recognised suffix.
    """
    if "_" in name:
        base, _, suffix = name.rpartition("_")
        chain = rules.DECLARED.get(suffix.upper())
        if chain:
            return base, chain
    return name, None


# --- N-terminal boundary --------------------------------------------------


def _germline_identity(seq: str, start: int, germlines) -> float:
    """Best fractional identity of seq[start:] against any germline FR1."""
    best = 0.0
    for germline in germlines:
        span = min(len(germline), 25, len(seq) - start)
        if span < 10:
            continue
        matches = sum(1 for i in range(span) if seq[start + i] == germline[i])
        best = max(best, matches / span)
    return best


def find_v_start(seq: str, v_end: int, spec: rules.ChainSpec) -> StartCall | None:
    """Score every plausible mature N-terminus and return the best.

    `v_end` is the inclusive index of the last V-domain residue.
    """
    min_len, max_len = spec.len_range
    typical_min, typical_max = spec.len_typical
    lowest = max(0, v_end - max_len + 1)
    highest = v_end - min_len + 1
    best: StartCall | None = None

    for start in range(lowest, highest + 1):
        length = v_end - start + 1
        score = 0.0
        reasons: list[str] = []

        if typical_min <= length <= typical_max:
            score += 1.5

        if seq[start] in spec.first_residues:
            score += 2.0
            reasons.append(f"FR1 opens {seq[start]}")

        # N-terminal motif (section 1b) at its canonical germline offset.
        off = spec.nterm_offset
        motif = spec.nterm_motif
        if seq[start + off : start + off + len(motif)] == motif:
            score += 3.5
            reasons.append(f"{motif} at residue {off + 1}")
        else:
            w0, w1 = spec.nterm_window
            if motif in seq[start + w0 : start + w1 + len(motif)]:
                score += 1.5
                reasons.append(f"{motif} near residue {off + 1}")

        # Conserved FR1 Cys (Kabat 23; section 4 cross-checks).
        cys_hit = next(
            (
                o
                for o in spec.cys_offsets
                if start + o <= v_end and seq[start + o] == "C"
            ),
            None,
        )
        if cys_hit is not None:
            score += 3.0
            reasons.append(f"FR1 Cys at residue {cys_hit + 1}")
        else:
            c0, c1 = spec.cys_window
            if "C" in seq[start + c0 : start + c1 + 1]:
                score += 1.5
                reasons.append("FR1 Cys near expected position")

        identity = _germline_identity(seq, start, spec.fr1_germline)
        score += 5.0 * identity
        if identity >= 0.5:
            reasons.append(f"germline FR1 identity {identity:.0%}")

        if best is None or score > best.score:
            best = StartCall(start=start, score=score, reasons=reasons)

    return best


# --- hypothesis evaluation ------------------------------------------------


def _constant_evidence(seq: str, v_end: int) -> list[str]:
    """Which constant-region opening, if any, follows the V domain."""
    tail = seq[v_end + 1 :]
    found = []
    if tail.startswith(rules.CH1_IGG1_START[:8]) or rules.CH1_IGG1_START in tail:
        found.append("heavy")
    if tail.startswith(rules.CL_KAPPA_START[:8]) or rules.CL_KAPPA_START in tail:
        found.append("kappa")
    if tail.startswith(rules.CL_LAMBDA_START[:9]) or rules.CL_LAMBDA_START in tail:
        found.append("lambda")
    return found


def evaluate(seq: str, chain: str) -> Hypothesis | None:
    """Best hypothesis that `seq` contains a V domain of class `chain`.

    Every FR4 anchor match is tried, not just the first: a full-length heavy
    chain can carry an incidental FGxG-like match in its Fc, and a sequence
    read past the V domain can carry more than one. Scoring the anchors
    against their implied FR1 lets the real one win.
    """
    spec = rules.SPECS[chain]
    candidates = [(m, True) for m in spec.fr4_strict.finditer(seq)]
    if not candidates:
        candidates = [(m, False) for m in spec.fr4_loose.finditer(seq)]
    if not candidates:
        return None

    best: Hypothesis | None = None
    for match, strict in candidates:
        fr4_index = match.start()
        v_end = fr4_index + spec.v_end_offset
        if v_end >= len(seq):
            continue  # FR4 runs off the end; domain is truncated
        start_call = find_v_start(seq, v_end, spec)
        if start_call is None:
            continue

        score = start_call.score
        reasons = list(start_call.reasons)

        if strict:
            score += 4.0
            reasons.insert(0, f"FR4 {seq[fr4_index:fr4_index + 4]} (canonical)")
        else:
            score += 1.0
            reasons.insert(0, f"FR4 {seq[fr4_index:fr4_index + 4]} (mutated FR4?)")

        constants = _constant_evidence(seq, v_end)
        if chain == "heavy" and "heavy" in constants:
            score += 4.0
            reasons.append("IgG1 CH1 begins at the boundary")
        elif chain == "light" and ("kappa" in constants or "lambda" in constants):
            score += 4.0
            reasons.append("light CL begins at the boundary")

        if chain == "heavy" and rules.HEAVY_CDR3_LEADIN.search(
            seq[max(0, v_end - 60) : v_end]
        ):
            score += 2.0
            reasons.append("heavy CDR3 lead-in (YYCAR/AK/TR)")

        if best is None or score > best.score:
            best = Hypothesis(
                chain=chain,
                fr4_index=fr4_index,
                strict=strict,
                v_start=start_call.start,
                v_end=v_end,
                score=score,
                reasons=reasons,
            )
    return best


# --- kappa / lambda -------------------------------------------------------


def call_light_subclass(seq: str, v_end: int, v_start: int) -> tuple[str | None, list[str], list[str]]:
    """Kappa or lambda, using one V-region and one C-region checkpoint.

    Section 3 of the knowledge base insists on at least two checkpoints, one
    from each region, because Vlambda-on-Ckappa hybrids are a real engineered
    format and a V-region-only call would mislabel them.
    """
    evidence: list[str] = []
    warnings: list[str] = []

    # 3b: last residue of the V domain. No exceptions among functional human J.
    last = seq[v_end]
    v_call = {"K": "kappa", "L": "lambda"}.get(last)
    if v_call:
        evidence.append(f"V ends {last} -> {v_call}")

    # 3c/3d: read into the constant region, past the shared AAPSV/FPPS.
    tail = seq[v_end + 1 :]
    c_call = None
    if rules.KAPPA_CL_DISCRIMINATOR in tail[:30]:
        c_call = "kappa"
        evidence.append("CL carries FIFPPS -> kappa")
    elif rules.LAMBDA_CL_DISCRIMINATOR in tail[:30]:
        c_call = "lambda"
        evidence.append("CL carries TLFPPS -> lambda")
    elif rules.CL_SHARED_MOTIF in tail[:20]:
        warnings.append(
            "CL present but too short to pass the shared AAPSV/FPPS stretch; "
            "kappa/lambda not confirmed from the C region"
        )

    # 3a: FR1 opening, weakest of the three.
    fr1 = seq[v_start : v_start + 6]
    if fr1 in rules.KAPPA_FR1_STARTS:
        evidence.append(f"FR1 {fr1} -> kappa")
    elif fr1 in rules.LAMBDA_FR1_STARTS:
        evidence.append(f"FR1 {fr1} -> lambda")

    if v_call and c_call and v_call != c_call:
        warnings.append(
            f"hybrid light chain: V region reads {v_call} but CL reads {c_call} "
            f"(V{v_call[0]}-on-C{c_call[0]} constructs are a real format) — "
            "V domain trimmed as given, but confirm the intended format"
        )
        return v_call, evidence, warnings

    if not c_call and v_call:
        warnings.append(
            "kappa/lambda called from the V region only; no constant region in "
            "the input to cross-check against"
        )
    return v_call or c_call, evidence, warnings


# --- top level ------------------------------------------------------------


def _confidence(hyp: Hypothesis) -> str:
    if hyp.strict and hyp.score >= 13.0:
        return "high"
    if hyp.score >= 9.0:
        return "medium"
    return "low"


def process(name: str, seq: str) -> ChainResult:
    """Type one sequence, check it against its header, and trim it."""
    base, declared = parse_declared(name)
    result = ChainResult(name=name, declared=declared, detected=None, subclass=None)

    # A trailing stop codon is not a residue. Trimming to the V domain would
    # discard it in any case; dropping it here keeps it out of the residue
    # validation below. An internal '*' is left in place so it is reported as
    # the truncation it represents.
    seq = seq.rstrip("*")
    result.input_seq = seq

    bad = sorted(set(seq) - rules.AMINO_ACIDS)
    if bad:
        result.errors.append(f"non-amino-acid characters in sequence: {''.join(bad)}")
        return result
    if not seq:
        result.errors.append("empty sequence")
        return result
    if declared is None:
        result.warnings.append(
            f"header has no _VH/_VL suffix; typed from sequence alone"
        )

    hypotheses = [h for h in (evaluate(seq, "heavy"), evaluate(seq, "light")) if h]
    if not hypotheses:
        result.errors.append(
            "no framework 4 anchor found (neither WGxG...VTV[SA] nor "
            "FGxGT...[KL]); this does not look like an antibody V domain"
        )
        return result

    hypotheses.sort(key=lambda h: h.score, reverse=True)
    best = hypotheses[0]
    result.detected = best.chain
    result.score = best.score
    result.evidence = best.reasons
    result.confidence = _confidence(best)

    if len(hypotheses) > 1:
        runner_up = hypotheses[1]
        if best.score - runner_up.score < 3.0:
            result.warnings.append(
                f"heavy/light call is close ({best.chain} {best.score:.1f} vs "
                f"{runner_up.chain} {runner_up.score:.1f}); verify with ANARCI"
            )

    if declared and result.detected != declared:
        result.errors.append(
            f"header declares {declared} but the sequence reads as "
            f"{result.detected} ({'; '.join(best.reasons)})"
        )
        return result

    if best.chain == "light":
        subclass, sub_evidence, sub_warnings = call_light_subclass(
            seq, best.v_end, best.v_start
        )
        result.subclass = subclass
        result.evidence.extend(sub_evidence)
        result.warnings.extend(sub_warnings)

    result.clean_seq = seq[best.v_start : best.v_end + 1]
    result.v_start = best.v_start
    result.v_end = best.v_end
    result.n_trimmed = best.v_start
    result.c_trimmed = len(seq) - best.v_end - 1

    spec = rules.SPECS[best.chain]
    typical_min, typical_max = spec.len_typical
    if not typical_min <= best.length <= typical_max:
        result.warnings.append(
            f"V domain is {best.length} aa, outside the usual "
            f"{typical_min}-{typical_max} for {best.chain}"
        )
    if result.confidence == "low":
        result.warnings.append(
            "weak N-terminal evidence; the mature start may be misplaced — "
            "confirm with ANARCI before using this sequence"
        )
    if result.n_trimmed and result.confidence != "high":
        result.warnings.append(
            f"{result.n_trimmed} residue(s) removed from the N-terminus on a "
            f"{result.confidence}-confidence boundary call"
        )

    return result


# --- storing a record the rules rejected ----------------------------------
#
# Typing is a judgement, and a judgement can be wrong: an engineered format,
# a chimera, a fragment the FR4 regexes do not recognise. So a failure is an
# answer the user may overrule -- but only in one direction. Overruling stores
# the sequence exactly as it came in, with no boundary call and no trim, and
# the row says so, because assembly splices whatever is stored verbatim.


def can_accept_as_given(result: ChainResult) -> tuple[str | None, str]:
    """Whether a rejected record could still be stored verbatim.

    Returns (chain, "") with the chain it would be stored as, or (None, why)
    when it could not be. Two failures are not a judgement and cannot be
    overruled: a sequence holding characters that are not amino acids is not a
    sequence, and an empty one has nothing to store. For the rest, something
    still has to say which chain it is -- the header's _VH/_VL suffix, or
    failing that the class the sequence read as.
    """
    if result.ok:
        return None, "this record passed the sanity check"
    if not result.input_seq:
        return None, "there is no sequence to store"
    bad = sorted(set(result.input_seq) - rules.AMINO_ACIDS)
    if bad:
        return None, (
            f"the sequence holds characters that are not amino acids "
            f"({''.join(bad)}); correct those first"
        )
    chain = result.declared or result.detected
    if not chain:
        return None, (
            "nothing says whether this is a heavy or a light chain: the header "
            "carries no _VH/_VL suffix and the sequence read as neither"
        )
    return chain, ""


#: Kappa or lambda has to come from the user for a light chain stored as
#: given: both checkpoints in section 3 need a V end, and an overruled record
#: has no boundary at all. Without it the row could not name a CL, so the
#: whole light chain and every express layout would be unbuildable.
LIGHT_SUBCLASSES = ("kappa", "lambda")


def accept_as_given(result: ChainResult, subclass: str | None = None) -> ChainResult:
    """Accept a rejected record exactly as it came in.

    Raises ValueError when the failure is not one that can be overruled. The
    chain is the one the header declared where it declared one -- the user is
    insisting on their own label, so it is their label that is honoured, not
    the class the sequence read as.

    A light chain needs `subclass`: "kappa" or "lambda". Nothing can be read
    off the sequence without a V end, so the user states it, and the row
    records that they stated it rather than that it was derived.
    """
    chain, reason = can_accept_as_given(result)
    if not chain:
        raise ValueError(f"'{result.name}' cannot be stored as given: {reason}.")

    subclass = (subclass or "").strip().lower() or None
    if chain == "light":
        if subclass not in LIGHT_SUBCLASSES:
            raise ValueError(
                f"'{result.name}' is a light chain, so storing it as given "
                f"needs you to say whether it is kappa or lambda: nothing can "
                f"be read off a sequence whose V domain was never found"
                + (f" (not {subclass!r})" if subclass else "")
            )
    elif subclass:
        raise ValueError(
            f"'{result.name}' is a heavy chain; kappa/lambda does not apply to it"
        )

    seq = result.input_seq
    accepted = ChainResult(
        name=result.name,
        declared=result.declared,
        detected=chain,
        # Stated, not derived: no boundary was established, so neither of
        # section 3's checkpoints could be read.
        subclass=subclass,
        clean_seq=seq,
        input_seq=seq,
        v_start=0,
        v_end=len(seq) - 1,
        n_trimmed=0,
        c_trimmed=0,
        score=result.score,
        confidence="none",
        evidence=list(result.evidence),
        accepted_as_given=True,
    )
    accepted.warnings.append(
        f"stored as given at your request; it did not pass the sanity check: "
        f"{'; '.join(result.errors)}"
    )
    accepted.warnings.append(
        f"not confirmed to be a {chain} V domain and not trimmed, so all "
        f"{len(seq)} residues are stored and assembly will splice them verbatim"
    )
    if subclass:
        accepted.warnings.append(
            f"{subclass} stated by you, not derived: neither the V-region nor "
            f"the C-region checkpoint of section 3 could be read, so the CL "
            f"this chain is given rests on your word alone"
        )
    accepted.warnings.extend(result.warnings)
    return accepted

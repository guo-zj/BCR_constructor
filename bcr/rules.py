"""Motifs and reference sequences transcribed from antibody-chain-rules.md.

Human antibodies only. Every constant below names the section of that document
it comes from, so the two can be diffed when the rules are revised.

These rules locate junctions for fast automated calls. Per the knowledge base,
anything production-relevant should still be confirmed with ANARCI or
IMGT/DomainGapAlign.
"""

import re

# --- Sections 2 and 4: framework 4 anchors --------------------------------
#
# Heavy: conserved Trp of WGxG; the VH domain ends at Trp + 10.
# Light: conserved Phe of FGxG; the VL domain ends at Phe + 9 (kappa and
# lambda alike).
#
# The heavy pattern ends VTV[SA] rather than TVSS, per the section 2 note that
# a literal TVSS silently drops murine IGHJ3 (...VTVSA). All six functional
# human IGHJ genes carry T at FR4 position 5 (WGxG-T), and every human light
# J segment carries T at the same spot, so the strict patterns pin that down.
# The loose patterns relax it for the somatically mutated FR4 the knowledge
# base warns is rare but real.

HEAVY_FR4_STRICT = re.compile(r"WG[A-Z]GT[A-Z]V[ST]V[SA]")
HEAVY_FR4_LOOSE = re.compile(r"WG[A-Z]G[A-Z]{2}V[ST]V[SA]")

# Light FR4 is FGxGT + 4 + terminal K (kappa) or L (lambda); see section 3b,
# which notes there are no exceptions among functional human J genes.
LIGHT_FR4_STRICT = re.compile(r"FG[A-Z]GT[A-Z]{4}[KL]")
LIGHT_FR4_LOOSE = re.compile(r"FG[A-Z]G[A-Z]{5}[KL]")

HEAVY_V_END_OFFSET = 10  # Trp + 10
LIGHT_V_END_OFFSET = 9  # Phe + 9

# --- Sections 1d, 2, 3c, 5: constant region openings ----------------------
# Presence of any of these means the input runs past the V domain.

CH1_IGG1_START = "ASTKGPSVFPLAP"  # section 5; nothing in a light chain resembles it
CL_KAPPA_START = "RTVAAPSVFIFPPS"
CL_LAMBDA_START = "GQPKAAPSVTLFPPS"

# IGHM CH1. NOT from the knowledge base -- that file covers typing and the
# V/C boundary and says nothing about IgM. Taken from the IGHM constant region
# shipped as an intrinsic part, and used only to recognise an IgM heavy
# constant region by its opening. Verify against IMGT before relying on it for
# anything more than a role label.
CH1_IGM_START = "GSASAPTLFPLVSC"

# Section 3c: the two CL openings share AAPSV and FPPS. What discriminates is
# the stretch between them, plus lambda's N-terminal GQPK extension.
CL_SHARED_MOTIF = "AAPSV"
KAPPA_CL_DISCRIMINATOR = "FIFPPS"
LAMBDA_CL_DISCRIMINATOR = "TLFPPS"

# --- Section 1b: germline FR1, used to score the mature N-terminus --------
# Shown as germline. Somatic hypermutation and degenerate cloning primers
# routinely alter the first few residues, so identity against these is scored
# as supporting evidence, never required.

HEAVY_FR1_GERMLINE = (
    "QVQLVQSGAEVKKPGASVKVSCKAS",  # IGHV1
    "EVQLVESGGGLVQPGGSLRLSCAAS",  # IGHV3
    "QVQLQESGPGLVKPSETLSLTCTVS",  # IGHV4
)
LIGHT_FR1_GERMLINE = (
    "DIQMTQSPSSLSASVGDRVTITC",  # IGKV1
    "EIVLTQSPGTLSLSPGERATLSC",  # IGKV3
    "QSVLTQPPSVSGAPGQRVTISC",  # IGLV1
)

# --- Section 3a: FR1 openings, for the kappa/lambda call ------------------
KAPPA_FR1_STARTS = ("DIQMTQ", "DIVMTQ", "EIVLTQ", "EIVMTQ")
LAMBDA_FR1_STARTS = ("QSVLTQ", "QSALTQ", "SYELTQ", "NFMLTQ", "QTVVTQ", "QPVLTQ", "QLVLTQ")

# --- Section 1d: CDR3 lead-in --------------------------------------------
# Heavy is usually YYCAR / YYCAK / YYCTR. Light is YYC/YFC/FYC followed by a
# non-AR residue, which is weaker, so only the heavy form is scored.
HEAVY_CDR3_LEADIN = re.compile(r"Y[YFH]C(?:AR|AK|TR)")

# --- Section 1d: expected lengths ----------------------------------------
# VH 98-128 linear residues (varies most with CDR-H3); VL 107-115. _TYPICAL is
# the knowledge base's own range and is what gets scored and warned about.
# _RANGE is only the window of start positions find_v_start is allowed to
# consider, and is deliberately wider: a start outside it is never scored at
# all, so a cap set too tightly does not merely lower confidence -- it forces
# a trim, silently amputating FR1 to make the domain fit.
#
# The heavy cap is 150 rather than the section's 128 because CDR-H3 length is
# unbounded in practice. C05 (anti-influenza, 26-residue CDR-H3) is a 138 aa
# VH that is already clean: at a cap of 135 its true N-terminus could not be
# considered, and the best start inside the window cut 15 residues and lost
# both the SG motif and the FR1 Cys. At 150 it types high-confidence untrimmed
# and still carries the "outside the expected 98-128" warning.
HEAVY_LEN_RANGE = (95, 150)
HEAVY_LEN_TYPICAL = (98, 128)
LIGHT_LEN_RANGE = (100, 125)
LIGHT_LEN_TYPICAL = (107, 115)


class ChainSpec:
    """Per-chain parameters for anchoring the mature N-terminus.

    nterm_motif sits at a strikingly fixed offset in germline FR1: for heavy,
    the SG of section 1b is at residues 7-8 (offset 6); for light, the TQ is at
    residues 5-6 (offset 4) in all three germline examples, kappa and lambda
    alike. The conserved FR1 Cys (Kabat 23, section 4 cross-checks) lands at
    residue 22 in every heavy germline above, and at 22-23 in light, since
    light FR1 length varies between IGKV and IGLV.
    """

    def __init__(
        self,
        name,
        fr4_strict,
        fr4_loose,
        v_end_offset,
        first_residues,
        nterm_motif,
        nterm_offset,
        nterm_window,
        cys_offsets,
        cys_window,
        fr1_germline,
        len_range,
        len_typical,
    ):
        self.name = name
        self.fr4_strict = fr4_strict
        self.fr4_loose = fr4_loose
        self.v_end_offset = v_end_offset
        self.first_residues = first_residues
        self.nterm_motif = nterm_motif
        self.nterm_offset = nterm_offset
        self.nterm_window = nterm_window
        self.cys_offsets = cys_offsets
        self.cys_window = cys_window
        self.fr1_germline = fr1_germline
        self.len_range = len_range
        self.len_typical = len_typical


HEAVY = ChainSpec(
    name="heavy",
    fr4_strict=HEAVY_FR4_STRICT,
    fr4_loose=HEAVY_FR4_LOOSE,
    v_end_offset=HEAVY_V_END_OFFSET,
    first_residues="QE",
    nterm_motif="SG",
    nterm_offset=6,
    nterm_window=(3, 9),
    cys_offsets=(21,),
    cys_window=(20, 23),
    fr1_germline=HEAVY_FR1_GERMLINE,
    len_range=HEAVY_LEN_RANGE,
    len_typical=HEAVY_LEN_TYPICAL,
)

LIGHT = ChainSpec(
    name="light",
    fr4_strict=LIGHT_FR4_STRICT,
    fr4_loose=LIGHT_FR4_LOOSE,
    v_end_offset=LIGHT_V_END_OFFSET,
    first_residues="DEQSN",
    nterm_motif="TQ",
    nterm_offset=4,
    nterm_window=(2, 7),
    cys_offsets=(21, 22),
    cys_window=(20, 24),
    fr1_germline=LIGHT_FR1_GERMLINE,
    len_range=LIGHT_LEN_RANGE,
    len_typical=LIGHT_LEN_TYPICAL,
)

SPECS = {"heavy": HEAVY, "light": LIGHT}

# Declared header suffix -> chain class
DECLARED = {"VH": "heavy", "VL": "light"}

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

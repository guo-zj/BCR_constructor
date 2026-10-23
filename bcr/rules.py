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
# VH 98-128 linear residues (varies most with CDR-H3); VL 107-115. Widened
# here so a genuinely unusual chain is reported rather than silently rejected;
# anything outside the knowledge base's own range is flagged.
HEAVY_LEN_RANGE = (95, 135)
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

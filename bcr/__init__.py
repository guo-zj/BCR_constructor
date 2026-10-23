"""BCR constructor — antibody chain typing, trimming and assembly.

Human antibodies only. Rules are transcribed from antibody-chain-rules.md.
"""

from .chains import ChainResult, process
from .fasta import Record, read_fasta, write_fasta

__all__ = ["ChainResult", "process", "Record", "read_fasta", "write_fasta"]

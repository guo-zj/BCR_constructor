"""Minimal FASTA reading and writing. Stdlib only."""

from dataclasses import dataclass

# A header's trailing token is taken as sequence, not description, when it is
# long and made only of residue letters. Some files put the sequence on the
# header line separated by a tab; a genuine description carries spaces,
# punctuation or digits and will not reach this threshold.
_SEQ_IN_HEADER_MIN = 20
_RESIDUE_LETTERS = set("ACDEFGHIKLMNPQRSTVWYXBZJU*-")


def _looks_like_sequence(token: str) -> bool:
    return (
        len(token) >= _SEQ_IN_HEADER_MIN
        and not set(token.upper()) - _RESIDUE_LETTERS
    )


@dataclass
class Record:
    name: str
    seq: str
    description: str = ""

    @property
    def header(self) -> str:
        return f"{self.name} {self.description}".strip()


def read_fasta(path) -> list[Record]:
    """Parse a FASTA file.

    Tolerates sequences wrapped over several lines, blank lines, and lowercase
    residues. Whitespace and '-' alignment gaps are dropped.

    A '*' stop character is kept, so a part file that ends its coding sequence
    with one carries it through into anything built from that part. Callers
    that want residues only should strip it themselves.
    """
    records: list[Record] = []
    name = None
    description = ""
    chunks: list[str] = []

    def flush():
        if name is not None:
            seq = "".join(chunks).upper().replace("-", "")
            records.append(Record(name=name, seq=seq, description=description))

    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith((";", "#")):
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                parts = header.split(None, 1)
                name = parts[0] if parts else ""
                description = ""
                chunks = []
                if len(parts) > 1:
                    trailing = "".join(parts[1].split())
                    if _looks_like_sequence(trailing):
                        chunks.append(trailing)
                    else:
                        description = parts[1]
            else:
                chunks.append("".join(line.split()))
    flush()
    return records


def write_fasta(path, records: list[Record], width: int = 0) -> None:
    """Write records, one sequence per line.

    `width` wraps sequence lines when set to a positive number; the default of
    0 writes each sequence on a single line, which keeps records greppable and
    easy to diff.
    """
    with open(path, "w", encoding="utf-8") as fh:
        for rec in records:
            fh.write(f">{rec.header}\n")
            seq = rec.seq
            if width and width > 0:
                for i in range(0, len(seq), width):
                    fh.write(seq[i : i + width] + "\n")
            else:
                fh.write(seq + "\n")

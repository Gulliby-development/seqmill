from dataclasses import dataclass


class SeqError(Exception):
    """An error that can be reported to a seqmill caller."""

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(message)


class FormatError(Exception):
    """A sequence-format error associated with a one-based input line."""

    def __init__(self, line: int, message: str) -> None:
        self.line = line
        self.message = message
        super().__init__(message)


@dataclass
class Record:
    name: str
    description: str
    seq: str
    qual: str | None = None


CODONS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

_BASE_INDEX = {"T": 0, "C": 1, "A": 2, "G": 3}
_COMPLEMENTS = str.maketrans("ACGTN", "TGCAN")


def show(text: str) -> str:
    """Return printable user text in quotes, replacing non-ASCII-printable chars."""
    printable = "".join(char if 33 <= ord(char) <= 126 else "?" for char in text)
    return f'"{printable}"'


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of an uppercase ACGTN sequence."""
    return seq.translate(_COMPLEMENTS)[::-1]


def translate(seq: str, to_stop: bool = False) -> str:
    """Translate complete codons in an uppercase ACGTN sequence."""
    amino_acids: list[str] = []
    for start in range(0, len(seq) - 2, 3):
        codon = seq[start : start + 3]
        if "N" in codon:
            amino_acid = "X"
        else:
            index = (
                _BASE_INDEX[codon[0]] * 16
                + _BASE_INDEX[codon[1]] * 4
                + _BASE_INDEX[codon[2]]
            )
            amino_acid = CODONS[index]
        if to_stop and amino_acid == "*":
            break
        amino_acids.append(amino_acid)
    return "".join(amino_acids)

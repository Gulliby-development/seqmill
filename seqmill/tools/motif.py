from seqmill.core.options import OptionSpec, parse_options
from seqmill.core.seq import SeqError, show
from seqmill.core.text import render_table
from seqmill.formats.files import Io, check_stdin_once, load_records

_IUPAC_BASES: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGTN"),
}
_IUPAC_LETTERS = frozenset(_IUPAC_BASES)
_IUPAC_COMPLEMENT = str.maketrans(
    {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "V": "B",
        "D": "H",
        "H": "D",
        "N": "N",
    }
)
_HEADER = ["name", "strand", "start", "end", "mismatches", "match"]


def _valid_pattern(pattern: str) -> bool:
    return (
        1 <= len(pattern) <= 64
        and pattern.isascii()
        and all(letter.upper() in _IUPAC_LETTERS for letter in pattern)
    )


def _reverse_complement_pattern(pattern: str) -> str:
    return pattern.translate(_IUPAC_COMPLEMENT)[::-1]


def _find_matches(
    sequence: str, pattern: str, max_mismatches: int, both_strands: bool
) -> list[tuple[int, int, list[str]]]:
    patterns = [("+", 0, pattern)]
    if both_strands:
        patterns.append(("-", 1, _reverse_complement_pattern(pattern)))

    matches: list[tuple[int, int, list[str]]] = []
    for strand, strand_order, oriented_pattern in patterns:
        for start in range(len(sequence) - len(oriented_pattern) + 1):
            window = sequence[start : start + len(oriented_pattern)]
            mismatch_count = 0
            for base, symbol in zip(window, oriented_pattern):
                if base not in _IUPAC_BASES[symbol]:
                    mismatch_count += 1
                    if mismatch_count > max_mismatches:
                        break
            if mismatch_count <= max_mismatches:
                matches.append(
                    (
                        start,
                        strand_order,
                        [
                            strand,
                            str(start + 1),
                            str(start + len(oriented_pattern)),
                            str(mismatch_count),
                            window,
                        ],
                    )
                )

    matches.sort(key=lambda match: (match[0], match[1]))
    return matches


def run(argv: list[str], io: Io) -> int:
    parsed = parse_options(
        argv,
        [
            OptionSpec("--mismatches", "int", default=0, lo=0, hi=8),
            OptionSpec("--both-strands", "flag", default=False),
        ],
    )

    args = parsed.args
    if not args:
        raise SeqError("missing PATTERN")
    pattern = args[0]
    if not _valid_pattern(pattern):
        raise SeqError(
            f'invalid pattern {show(pattern)} '
            "(expected 1 to 64 letters from ACGTNRYSWKMBDHV)"
        )
    pattern = pattern.upper()

    paths = args[1:]
    if not paths:
        raise SeqError("missing FILE")
    check_stdin_once(paths)

    max_mismatches = parsed.values["--mismatches"]
    both_strands = parsed.values["--both-strands"]
    rows: list[list[str]] = []
    for path in paths:
        loaded = load_records(io, path)
        for record in loaded.records:
            matches = _find_matches(
                record.seq, pattern, max_mismatches, both_strands
            )
            for _, _, fields in matches:
                rows.append([record.name, *fields])

    if not rows:
        return 1

    lines = render_table(_HEADER, rows, {2, 3, 4})
    io.stdout.write("\n".join(lines) + "\n")
    return 0

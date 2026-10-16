from array import array

from seqmill.core.options import OptionSpec, parse_options
from seqmill.core.seq import SeqError, show
from seqmill.core.text import fixed1
from seqmill.formats.files import Io, check_stdin_once, load_records

_SPECS = [
    OptionSpec("--mode", "choice", "global", choices=("global", "local")),
    OptionSpec("--match", "int", 1, 1, 100),
    OptionSpec("--mismatch", "int", -1, -100, 0),
    OptionSpec("--gap", "int", -2, -100, -1),
    OptionSpec("--width", "int", 60, 1, 1000),
]


def _load_first(io: Io, path: str) -> tuple[str, str]:
    loaded = load_records(io, path)
    if not loaded.records:
        raise SeqError(f"{loaded.display}: no sequences")

    record = loaded.records[0]
    if not record.seq:
        raise SeqError(f"{loaded.display}: sequence {show(record.name)} is empty")
    if len(record.seq) > 2000:
        raise SeqError(f"{loaded.display}: sequence {show(record.name)} is longer than 2000 bases")
    return record.name, record.seq


def _score_matrix(
    a_seq: str,
    b_seq: str,
    mode: str,
    match: int,
    mismatch: int,
    gap: int,
) -> tuple[list[array], int, int, int]:
    rows = len(a_seq) + 1
    columns = len(b_seq) + 1
    scores = [array("q", [0]) * columns for _ in range(rows)]

    if mode == "global":
        scores[0] = array("q", (column * gap for column in range(columns)))
        for row in range(1, rows):
            scores[row][0] = row * gap

    best_score = 0
    best_row = 0
    best_column = 0
    for row in range(1, rows):
        for column in range(1, columns):
            substitution = match if a_seq[row - 1] == b_seq[column - 1] else mismatch
            diagonal = scores[row - 1][column - 1] + substitution
            up = scores[row - 1][column] + gap
            left = scores[row][column - 1] + gap
            score = max(diagonal, up, left)
            if mode == "local":
                score = max(0, score)
            scores[row][column] = score

            if mode == "local" and score > best_score:
                best_score = score
                best_row = row
                best_column = column

    if mode == "global":
        best_row = len(a_seq)
        best_column = len(b_seq)
        best_score = scores[best_row][best_column]
    return scores, best_score, best_row, best_column


def _traceback(
    a_seq: str,
    b_seq: str,
    scores: list[array],
    mode: str,
    match: int,
    mismatch: int,
    gap: int,
    end_row: int,
    end_column: int,
) -> tuple[str, str, int, int]:
    row = end_row
    column = end_column
    aligned_a: list[str] = []
    aligned_b: list[str] = []

    while row > 0 or column > 0:
        current = scores[row][column]
        if mode == "local" and current == 0:
            break

        if row > 0 and column > 0:
            substitution = match if a_seq[row - 1] == b_seq[column - 1] else mismatch
            if current == scores[row - 1][column - 1] + substitution:
                aligned_a.append(a_seq[row - 1])
                aligned_b.append(b_seq[column - 1])
                row -= 1
                column -= 1
                continue

        if row > 0 and current == scores[row - 1][column] + gap:
            aligned_a.append(a_seq[row - 1])
            aligned_b.append("-")
            row -= 1
            continue

        aligned_a.append("-")
        aligned_b.append(b_seq[column - 1])
        column -= 1

    return "".join(reversed(aligned_a)), "".join(reversed(aligned_b)), row, column


def _render(
    a_name: str,
    b_name: str,
    score: int,
    aligned_a: str,
    aligned_b: str,
    a_start: int,
    b_start: int,
    width: int,
) -> str:
    identity = sum(
        1 for a_base, b_base in zip(aligned_a, aligned_b) if a_base != "-" and a_base == b_base
    )
    length = len(aligned_a)
    markers = "".join(
        "|" if a_base == b_base and a_base != "-" else "-" if "-" in (a_base, b_base) else "x"
        for a_base, b_base in zip(aligned_a, aligned_b)
    )
    lines = [
        f"score: {score}",
        f"length: {length}",
        f"identity: {identity}/{length} ({fixed1(100 * identity, length)}%)",
        f"A: {a_name} {a_start}-{a_start + sum(base != '-' for base in aligned_a) - 1}",
        f"B: {b_name} {b_start}-{b_start + sum(base != '-' for base in aligned_b) - 1}",
    ]
    for start in range(0, length, width):
        end = start + width
        lines.extend(("", aligned_a[start:end], markers[start:end], aligned_b[start:end]))
    return "\n".join(lines) + "\n"


def run(argv: list[str], io: Io) -> int:
    parsed = parse_options(argv, _SPECS)
    if not parsed.args:
        raise SeqError("missing A")
    if len(parsed.args) == 1:
        raise SeqError("missing B")
    if len(parsed.args) > 2:
        raise SeqError("too many arguments")
    check_stdin_once(parsed.args)

    mode = parsed.values["--mode"]
    match = parsed.values["--match"]
    mismatch = parsed.values["--mismatch"]
    gap = parsed.values["--gap"]
    width = parsed.values["--width"]

    name_a, sequence_a = _load_first(io, parsed.args[0])
    name_b, sequence_b = _load_first(io, parsed.args[1])

    scores, score, end_row, end_column = _score_matrix(
        sequence_a, sequence_b, mode, match, mismatch, gap
    )
    if mode == "local" and score == 0:
        io.stdout.write("score: 0\nno alignment\n")
        return 1

    aligned_a, aligned_b, start_row, start_column = _traceback(
        sequence_a,
        sequence_b,
        scores,
        mode,
        match,
        mismatch,
        gap,
        end_row,
        end_column,
    )
    output = _render(
        name_a,
        name_b,
        score,
        aligned_a,
        aligned_b,
        start_row + 1,
        start_column + 1,
        width,
    )
    io.stdout.write(output)
    return 0

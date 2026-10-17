from seqmill.core.options import OptionSpec, parse_options
from seqmill.core.seq import Record, SeqError, reverse_complement
from seqmill.core.text import render_table
from seqmill.formats.files import Io, check_stdin_once, load_records

_STOPS = {"TAA", "TAG", "TGA"}
_HEADER = ["name", "frame", "start", "end", "aa"]


def _scan_frame(
    sequence: str, offset: int, frame: int, minimum: int
) -> list[tuple[int, int, int, int]]:
    found: list[tuple[int, int, int, int]] = []
    open_start: int | None = None

    for codon_start in range(offset, len(sequence) - 2, 3):
        codon = sequence[codon_start : codon_start + 3]
        if open_start is None:
            if codon == "ATG":
                open_start = codon_start
            continue

        if codon not in _STOPS:
            continue

        amino_acids = (codon_start - open_start) // 3
        if amino_acids >= minimum:
            if frame > 0:
                start = open_start + 1
                end = codon_start + 3
            else:
                start = len(sequence) - (codon_start + 3) + 1
                end = len(sequence) - (open_start + 1) + 1
            found.append((start, end, frame, amino_acids))
        open_start = None

    return found


def _record_orfs(record: Record, minimum: int, strand: str) -> list[list[str]]:
    found: list[tuple[int, int, int, int]] = []
    scan_forward = strand != "reverse"
    scan_reverse = strand != "forward"
    reverse = reverse_complement(record.seq) if scan_reverse else ""

    for offset in range(3):
        if scan_forward:
            found.extend(_scan_frame(record.seq, offset, offset + 1, minimum))
        if scan_reverse:
            found.extend(_scan_frame(reverse, offset, -(offset + 1), minimum))

    found.sort(key=lambda item: (item[0], item[1], item[2]))
    return [
        [record.name, f"{frame:+d}", str(start), str(end), str(amino_acids)]
        for start, end, frame, amino_acids in found
    ]


def run(argv: list[str], io: Io) -> int:
    parsed = parse_options(
        argv,
        [
            OptionSpec("--min", "int", default=30, lo=1, hi=100000),
            OptionSpec(
                "--strand", "choice", default="both", choices=("both", "forward", "reverse")
            ),
        ],
    )
    paths = parsed.args
    if not paths:
        raise SeqError("missing FILE")
    check_stdin_once(paths)

    rows: list[list[str]] = []
    minimum = parsed.values["--min"]
    strand = parsed.values["--strand"]
    for path in paths:
        loaded = load_records(io, path)
        for record in loaded.records:
            rows.extend(_record_orfs(record, minimum, strand))

    if not rows:
        return 1

    lines = render_table(_HEADER, rows, {1, 2, 3, 4})
    io.stdout.write("\n".join(lines) + "\n")
    return 0

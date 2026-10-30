from seqmill.core.options import OptionSpec, parse_options
from seqmill.core.seq import SeqError
from seqmill.core.text import fixed1, render_table
from seqmill.formats.files import Io, check_stdin_once, load_records

_HEADER = ["file", "records", "bases", "min", "max", "mean", "n50", "gc%", "n"]


def run(argv: list[str], io: Io) -> int:
    parsed = parse_options(argv, [OptionSpec("--no-header", "flag", False)])
    if not parsed.args:
        raise SeqError("missing FILE")
    check_stdin_once(parsed.args)

    loaded_files = [load_records(io, path) for path in parsed.args]
    rows: list[list[str]] = []
    for loaded in loaded_files:
        lengths = [len(record.seq) for record in loaded.records]
        record_count = len(loaded.records)
        bases = sum(lengths)
        n_bases = sum(record.seq.count("N") for record in loaded.records)
        gc_bases = sum(
            record.seq.count("G") + record.seq.count("C")
            for record in loaded.records
        )
        canonical_bases = bases - n_bases

        if record_count:
            minimum = str(min(lengths))
            maximum = str(max(lengths))
            mean = fixed1(bases, record_count)
        else:
            minimum = maximum = mean = "-"

        if bases:
            cumulative = 0
            n50 = "-"
            for length in sorted(lengths, reverse=True):
                cumulative += length
                if cumulative * 2 >= bases:
                    n50 = str(length)
                    break
        else:
            n50 = "-"

        gc_percent = fixed1(100 * gc_bases, canonical_bases) if canonical_bases else "-"
        rows.append(
            [
                loaded.display,
                str(record_count),
                str(bases),
                minimum,
                maximum,
                mean,
                n50,
                gc_percent,
                str(n_bases),
            ]
        )

    lines = render_table(_HEADER, rows, set(range(1, len(_HEADER))))
    if parsed.values["--no-header"]:
        lines = lines[1:]
    output = "\n".join(lines) + "\n"
    io.stdout.write(output)
    return 0

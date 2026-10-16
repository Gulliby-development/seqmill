from seqmill.core.options import OptionSpec, parse_options
from seqmill.core.seq import Record, SeqError, reverse_complement, translate
from seqmill.formats.fasta import format_fasta
from seqmill.formats.fastq import format_fastq
from seqmill.formats.files import Io, check_stdin_once, load_records


def run_revcomp(argv: list[str], io: Io) -> int:
    specs = [OptionSpec("--width", "int", 60, 0, 1000)]
    parsed = parse_options(argv, specs)
    if not parsed.args:
        raise SeqError("missing FILE")
    check_stdin_once(parsed.args)

    loaded_files = [load_records(io, path) for path in parsed.args]
    width = int(parsed.values["--width"])
    output: list[str] = []

    for loaded in loaded_files:
        records = [
            Record(
                record.name,
                record.description,
                reverse_complement(record.seq),
                record.qual[::-1] if record.qual is not None else None,
            )
            for record in loaded.records
        ]
        if loaded.kind == "fastq":
            output.append(format_fastq(records))
        else:
            output.append(format_fasta(records, width))

    io.stdout.write("".join(output))
    return 0


def run_translate(argv: list[str], io: Io) -> int:
    specs = [
        OptionSpec("--frame", "choice", "1", choices=("1", "2", "3", "-1", "-2", "-3")),
        OptionSpec("--to-stop", "flag", False),
        OptionSpec("--width", "int", 60, 0, 1000),
    ]
    parsed = parse_options(argv, specs)
    if not parsed.args:
        raise SeqError("missing FILE")
    check_stdin_once(parsed.args)

    loaded_files = [load_records(io, path) for path in parsed.args]
    frame = str(parsed.values["--frame"])
    offset = abs(int(frame)) - 1
    to_stop = bool(parsed.values["--to-stop"])
    width = int(parsed.values["--width"])
    translated_records: list[Record] = []

    for loaded in loaded_files:
        for record in loaded.records:
            sequence = record.seq
            if frame.startswith("-"):
                sequence = reverse_complement(sequence)
            protein = translate(sequence[offset:], to_stop=to_stop)
            translated_records.append(Record(record.name, record.description, protein))

    io.stdout.write(format_fasta(translated_records, width))
    return 0

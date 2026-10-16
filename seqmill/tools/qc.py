from typing import cast

from seqmill.core.options import OptionSpec, parse_options
from seqmill.core.seq import Record, SeqError
from seqmill.formats.fastq import format_fastq
from seqmill.formats.files import Io, check_stdin_once, load_records

_OPTION_SPECS = [
    OptionSpec("--trim-end", "int", default=0, lo=0, hi=93),
    OptionSpec("--min-length", "int", default=1, lo=0, hi=1_000_000),
    OptionSpec("--min-quality", "int", default=0, lo=0, hi=93),
]


def run(argv: list[str], io: Io) -> int:
    parsed = parse_options(argv, _OPTION_SPECS)
    paths = parsed.args
    if not paths:
        raise SeqError("missing FILE")
    check_stdin_once(paths)

    trim_end = cast(int, parsed.values["--trim-end"])
    min_length = cast(int, parsed.values["--min-length"])
    min_quality = cast(int, parsed.values["--min-quality"])

    kept_records: list[Record] = []
    total_reads = 0
    too_short = 0
    low_quality = 0
    trimmed_bases = 0

    for path in paths:
        loaded = load_records(io, path)
        if loaded.records and loaded.kind != "fastq":
            raise SeqError(f"{loaded.display}: not a FASTQ file")

        for record in loaded.records:
            total_reads += 1
            quality = record.qual if record.qual is not None else ""
            retained_length = len(record.seq)
            while (
                retained_length > 0
                and ord(quality[retained_length - 1]) - 33 < trim_end
            ):
                retained_length -= 1
            trimmed_bases += len(record.seq) - retained_length

            sequence = record.seq[:retained_length]
            retained_quality = quality[:retained_length]
            if retained_length < min_length:
                too_short += 1
            elif sum(ord(character) - 33 for character in retained_quality) < (
                min_quality * retained_length
            ):
                low_quality += 1
            else:
                kept_records.append(
                    Record(record.name, record.description, sequence, retained_quality)
                )

    io.stdout.write(format_fastq(kept_records))
    io.stderr.write(
        f"kept {len(kept_records)} of {total_reads} reads: {too_short} too short, "
        f"{low_quality} low quality, {trimmed_bases} bases trimmed\n"
    )
    return 0

from seqmill.core.seq import FormatError, Record, show
from seqmill.formats.fastq import _header_fields

_VALID_BASES = frozenset("ACGTNacgtn")


def parse_fasta(lines: list[str]) -> list[Record]:
    records: list[Record] = []
    name = ""
    description = ""
    sequence_parts: list[str] = []

    for line_number, line in enumerate(lines, start=1):
        if not line.strip(" \t"):
            continue

        if line.startswith(">"):
            next_name, next_description = _header_fields(line)
            if not next_name:
                raise FormatError(line_number, "empty name")
            if name:
                records.append(Record(name, description, "".join(sequence_parts)))
            name = next_name
            description = next_description
            sequence_parts = []
            continue

        if not name:
            raise FormatError(line_number, "sequence data before the first header")

        for character in line:
            if character not in _VALID_BASES:
                raise FormatError(
                    line_number,
                    f"invalid character {show(character)} in sequence",
                )
        sequence_parts.append(line.upper())

    if name:
        records.append(Record(name, description, "".join(sequence_parts)))
    return records


def format_fasta(records: list[Record], width: int) -> str:
    output: list[str] = []
    for record in records:
        header = ">" + record.name
        if record.description:
            header += " " + record.description
        output.append(header)
        if width == 0:
            if record.seq:
                output.append(record.seq)
        else:
            for start in range(0, len(record.seq), width):
                output.append(record.seq[start : start + width])
    if not output:
        return ""
    return "\n".join(output) + "\n"

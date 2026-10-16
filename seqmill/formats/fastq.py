from seqmill.core.seq import FormatError, Record, show

_VALID_BASES = frozenset("ACGTNacgtn")


def _header_fields(header: str) -> tuple[str, str]:
    text = header[1:]
    separator = min(
        (index for index in (text.find(" "), text.find("\t")) if index >= 0),
        default=len(text),
    )
    name = text[:separator]
    description = text[separator:].strip(" \t")
    return name, description


def _invalid_character(text: str, allowed: frozenset[str]) -> str | None:
    for character in text:
        if character not in allowed:
            return character
    return None


def parse_fastq(lines: list[str]) -> list[Record]:
    records: list[Record] = []
    index = 0
    while index < len(lines):
        record_line = index + 1
        header = lines[index]
        if not header.startswith("@"):
            raise FormatError(record_line, "expected @ at the start of a record")
        if len(lines) - index < 4:
            raise FormatError(record_line, "incomplete record")

        name, description = _header_fields(header)
        if not name:
            raise FormatError(record_line, "empty name")

        sequence = lines[index + 1]
        invalid_base = _invalid_character(sequence, _VALID_BASES)
        if invalid_base is not None:
            raise FormatError(
                record_line + 1,
                f"invalid character {show(invalid_base)} in sequence",
            )

        if not lines[index + 2].startswith("+"):
            raise FormatError(record_line + 2, "expected + line")

        quality = lines[index + 3]
        invalid_quality = _invalid_character(
            quality,
            frozenset(chr(codepoint) for codepoint in range(33, 127)),
        )
        if invalid_quality is not None:
            raise FormatError(
                record_line + 3,
                f"invalid character {show(invalid_quality)} in quality",
            )
        normalized_sequence = sequence.upper()
        if len(normalized_sequence) != len(quality):
            raise FormatError(
                record_line + 3,
                "sequence and quality have different lengths",
            )

        records.append(Record(name, description, normalized_sequence, quality))
        index += 4
    return records


def format_fastq(records: list[Record]) -> str:
    output: list[str] = []
    for record in records:
        header = "@" + record.name
        if record.description:
            header += " " + record.description
        output.extend((header, record.seq, "+", record.qual or ""))
    if not output:
        return ""
    return "\n".join(output) + "\n"

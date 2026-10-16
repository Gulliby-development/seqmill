from dataclasses import dataclass
from typing import BinaryIO, TextIO

from seqmill.core.seq import FormatError, Record, SeqError
from seqmill.formats.fasta import parse_fasta
from seqmill.formats.fastq import parse_fastq


@dataclass
class Io:
    stdin: BinaryIO
    stdout: TextIO
    stderr: TextIO


@dataclass
class Loaded:
    display: str
    kind: str
    records: list[Record]


def display_path(path: str) -> str:
    if path == "-":
        return "<stdin>"
    return path.replace("\\", "/")


def split_lines(text: str) -> list[str]:
    text = text.removeprefix("\ufeff")
    if not text:
        return []

    ends_with_lf = text.endswith("\n")
    lines = text.split("\n")
    if ends_with_lf:
        lines.pop()

    result: list[str] = []
    for index, line in enumerate(lines):
        has_line_ending = index < len(lines) - 1 or ends_with_lf
        if has_line_ending and line.endswith("\r"):
            line = line[:-1]
        result.append(line)
    return result


def read_text(io: Io, path: str) -> tuple[str, str]:
    display = display_path(path)
    try:
        if path == "-":
            contents = io.stdin.read()
        else:
            with open(path, "rb") as source:
                contents = source.read()
    except (OSError, ValueError, UnicodeError) as error:
        raise SeqError(f"{display}: cannot read the file") from error

    try:
        text = contents.decode("utf-8")
    except UnicodeError as error:
        raise SeqError(f"{display}: input is not valid UTF-8") from error
    return display, text


def check_stdin_once(paths: list[str]) -> None:
    if paths.count("-") > 1:
        raise SeqError("standard input given twice")


def load_records(io: Io, path: str) -> Loaded:
    display, text = read_text(io, path)
    lines = split_lines(text)
    is_fastq = bool(lines and lines[0].startswith("@"))
    kind = "fastq" if is_fastq else "fasta"

    try:
        records = parse_fastq(lines) if is_fastq else parse_fasta(lines)
    except FormatError as error:
        raise SeqError(f"{display}:{error.line}: {error.message}") from error
    return Loaded(display, kind, records)

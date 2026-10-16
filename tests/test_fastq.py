from dataclasses import fields, is_dataclass
from inspect import signature

import pytest

from seqmill.core.seq import FormatError, Record, SeqError, show
from seqmill.formats import fastq, files
from seqmill.formats.files import Io, Loaded, split_lines


def test_fastq_api_reuses_core_types_and_exposes_contract_signatures() -> None:
    assert fastq.Record is Record
    assert fastq.FormatError is FormatError
    assert files.SeqError is SeqError
    assert fastq.show is show
    assert list(signature(fastq.parse_fastq).parameters) == ["lines"]
    assert list(signature(fastq.format_fastq).parameters) == ["records"]
    assert signature(fastq.parse_fastq).return_annotation == list[Record]
    assert signature(fastq.format_fastq).return_annotation is str

    assert is_dataclass(Record)
    assert [field.name for field in fields(Record)] == ["name", "description", "seq", "qual"]
    assert is_dataclass(Io)
    assert [field.name for field in fields(Io)] == ["stdin", "stdout", "stderr"]
    assert is_dataclass(Loaded)
    assert [field.name for field in fields(Loaded)] == ["display", "kind", "records"]


@pytest.mark.parametrize(
    ("lines", "expected"),
    [
        ([], []),
        (
            ["@a desc", "acgtn", "+ignored", "!~ABC"],
            [Record("a", "desc", "ACGTN", "!~ABC")],
        ),
        (
            ["@a\t desc\t", "", "+", ""],
            [Record("a", "desc", "", "")],
        ),
        (
            ["@a", "A", "+", "!", "@a", "C", "+", "~"],
            [Record("a", "", "A", "!"), Record("a", "", "C", "~")],
        ),
        (["@é", "", "+anything", ""], [Record("é", "", "", "")]),
        (["@a", "A", "+", "@"], [Record("a", "", "A", "@")]),
        (
            ["@a\t  lane 1 \t", "a", "+ignored text", "~"],
            [Record("a", "lane 1", "A", "~")],
        ),
    ],
)
def test_parse_fastq_valid_records(lines: list[str], expected: list[Record]) -> None:
    assert fastq.parse_fastq(lines) == expected


@pytest.mark.parametrize(
    ("lines", "line", "message"),
    [
        ([""], 1, "expected @ at the start of a record"),
        (["x", "R", "x", " "], 1, "expected @ at the start of a record"),
        (["@"], 1, "incomplete record"),
        (["@a", "R", "x"], 1, "incomplete record"),
        (["@", "R", "x", " "], 1, "empty name"),
        (["@ a", "A", "+", "!"], 1, "empty name"),
        (["@a", "rZ", "x", " "], 2, 'invalid character "r" in sequence'),
        (["@a", "A ", "+", "!"], 2, 'invalid character "?" in sequence'),
        (["@a", "ſ", "+", "!"], 2, 'invalid character "?" in sequence'),
        (["@a", "A", "x", " "], 3, "expected + line"),
        (["@a", "A", "x", "\x7f"], 3, "expected + line"),
        (["@a", "A", "+", " "], 4, 'invalid character "?" in quality'),
        (["@a", "AC", "+", " "], 4, 'invalid character "?" in quality'),
        (["@a", "A", "+", "\x7f"], 4, 'invalid character "?" in quality'),
        (["@a", "A", "+", "é"], 4, 'invalid character "?" in quality'),
        (["@a", "A", "+", "\té"], 4, 'invalid character "?" in quality'),
        (["@a", "AC", "+", "!"], 4, "sequence and quality have different lengths"),
        (["@a", "A", "+", "!!"], 4, "sequence and quality have different lengths"),
        (["@a", "AC", "+", "é"], 4, 'invalid character "?" in quality'),
        (["@a", "A", "+", "!", ""], 5, "expected @ at the start of a record"),
        (
            ["@ok", "A", "+", "!", "@bad", "AC", "+", "é"],
            8,
            'invalid character "?" in quality',
        ),
    ],
)
def test_parse_fastq_reports_first_error_with_exact_line_and_message(
    lines: list[str], line: int, message: str
) -> None:
    with pytest.raises(FormatError) as caught:
        fastq.parse_fastq(lines)

    assert caught.value.line == line
    assert caught.value.message == message
    assert str(caught.value) == message


@pytest.mark.parametrize(
    ("records", "expected"),
    [
        ([], ""),
        ([Record("a", "", "A", "!")], "@a\nA\n+\n!\n"),
        ([Record("a", "desc", "AC", "!~")], "@a desc\nAC\n+\n!~\n"),
        ([Record("a", "", "", "")], "@a\n\n+\n\n"),
        (
            [Record("a", "lane 1", "AC", "!~"), Record("b", "", "GT", "AB")],
            "@a lane 1\nAC\n+\n!~\n@b\nGT\n+\nAB\n",
        ),
    ],
)
def test_format_fastq_exact_lf_output(records: list[Record], expected: str) -> None:
    assert fastq.format_fastq(records) == expected


def test_fastq_records_round_trip_through_writer_and_line_splitter() -> None:
    records = [
        Record("r1", "lane 1", "ACGTN", "!~ABC"),
        Record("r2", "", "", ""),
    ]

    serialized = fastq.format_fastq(records)

    assert serialized == "@r1 lane 1\nACGTN\n+\n!~ABC\n@r2\n\n+\n\n"
    assert fastq.parse_fastq(split_lines(serialized)) == records

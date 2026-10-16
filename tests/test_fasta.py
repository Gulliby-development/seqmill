from dataclasses import fields, is_dataclass
from inspect import signature

import pytest

from seqmill.core.seq import FormatError, Record, SeqError, show
from seqmill.formats import fasta, files
from seqmill.formats.files import Io, Loaded, split_lines


def test_fasta_api_reuses_core_types_and_exposes_contract_signatures() -> None:
    assert fasta.Record is Record
    assert fasta.FormatError is FormatError
    assert files.SeqError is SeqError
    assert fasta.show is show
    assert list(signature(fasta.parse_fasta).parameters) == ["lines"]
    assert list(signature(fasta.format_fasta).parameters) == ["records", "width"]
    assert signature(fasta.parse_fasta).return_annotation == list[Record]
    assert signature(fasta.format_fasta).return_annotation is str

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
        ([">a"], [Record("a", "", "")]),
        ([">a desc", "ac", "gTn"], [Record("a", "desc", "ACGTN")]),
        ([">a\t  two words\t", "A"], [Record("a", "two words", "A")]),
        ([">a", ">a", "C"], [Record("a", "", ""), Record("a", "", "C")]),
        ([" \t", ">a", "\t", "AC"], [Record("a", "", "AC")]),
        ([">é info", "N"], [Record("é", "info", "N")]),
    ],
)
def test_parse_fasta_valid_records(lines: list[str], expected: list[Record]) -> None:
    assert fasta.parse_fasta(lines) == expected


@pytest.mark.parametrize(
    ("lines", "line", "message"),
    [
        ([">"], 1, "empty name"),
        (["> a"], 1, "empty name"),
        ([">\ta"], 1, "empty name"),
        (["AC"], 1, "sequence data before the first header"),
        (["R", ">a"], 1, "sequence data before the first header"),
        (["", ">a", " \t", "ACrZ"], 4, 'invalid character "r" in sequence'),
        ([">a", "AC G"], 2, 'invalid character "?" in sequence'),
        ([">a", "ſ"], 2, 'invalid character "?" in sequence'),
    ],
)
def test_parse_fasta_reports_first_error_with_exact_line_and_message(
    lines: list[str], line: int, message: str
) -> None:
    with pytest.raises(FormatError) as caught:
        fasta.parse_fasta(lines)

    assert caught.value.line == line
    assert caught.value.message == message
    assert str(caught.value) == message


@pytest.mark.parametrize(
    ("records", "width", "expected"),
    [
        ([], 0, ""),
        ([Record("a", "", "")], 60, ">a\n"),
        ([Record("a", "desc", "ACGTN")], 2, ">a desc\nAC\nGT\nN\n"),
        ([Record("a", "", "ACG")], 0, ">a\nACG\n"),
        ([Record("a", "", "AC")], 1, ">a\nA\nC\n"),
        ([Record("a", "", "EFMAK*X")], 3, ">a\nEFM\nAK*\nX\n"),
        ([Record("a", "", "ACG")], 60, ">a\nACG\n"),
        ([Record("a", "", "ACG")], 1000, ">a\nACG\n"),
        (
            [Record("a", "first", "AC"), Record("b", "", "GT")],
            2,
            ">a first\nAC\n>b\nGT\n",
        ),
    ],
)
def test_format_fasta_exact_output(
    records: list[Record], width: int, expected: str
) -> None:
    assert fasta.format_fasta(records, width) == expected


def test_fasta_normalized_records_round_trip_through_writer_and_line_splitter() -> None:
    records = fasta.parse_fasta([">a lane 1", "ac", "GT", ">a", "n"])

    serialized = fasta.format_fasta(records, 2)

    assert serialized == ">a lane 1\nAC\nGT\n>a\nN\n"
    assert fasta.parse_fasta(split_lines(serialized)) == records

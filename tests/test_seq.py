from dataclasses import fields, is_dataclass

import pytest

from seqmill.core.seq import (
    CODONS,
    FormatError,
    Record,
    SeqError,
    reverse_complement,
    show,
    translate,
)


def test_seq_error_exposes_message_and_inherits_exception() -> None:
    error = SeqError("bad")

    assert isinstance(error, Exception)
    assert error.message == "bad"
    assert str(error) == "bad"


def test_format_error_exposes_line_and_message_and_inherits_exception() -> None:
    error = FormatError(7, "bad")

    assert isinstance(error, Exception)
    assert error.line == 7
    assert error.message == "bad"
    assert str(error) == "bad"


def test_record_is_a_mutable_dataclass_with_the_contract_field_order() -> None:
    assert is_dataclass(Record)
    assert [field.name for field in fields(Record)] == ["name", "description", "seq", "qual"]

    record = Record("r", "", "ACGN")
    assert record.name == "r"
    assert record.description == ""
    assert record.seq == "ACGN"
    assert record.qual is None

    record.seq = "TGCA"
    assert record.seq == "TGCA"


def test_record_accepts_fastq_quality() -> None:
    record = Record("q", "lane 1", "AC", "!~")

    assert record == Record("q", "lane 1", "AC", "!~")
    assert record.qual == "!~"


def test_codon_table_is_the_exact_t_c_a_g_ordered_literal() -> None:
    assert CODONS == "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    assert len(CODONS) == 64


@pytest.mark.parametrize(
    ("sequence", "expected"),
    [
        ("", ""),
        ("A", "T"),
        ("C", "G"),
        ("G", "C"),
        ("T", "A"),
        ("N", "N"),
        ("ACGTN", "NACGT"),
        ("AAACCCGGGTTT", "AAACCCGGGTTT"),
        ("ACGTTGCA", "TGCAACGT"),
    ],
)
def test_reverse_complement(sequence: str, expected: str) -> None:
    assert reverse_complement(sequence) == expected


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("", '""'),
        ("!", '"!"'),
        ("~", '"~"'),
        ('"', '"""'),
        ("\\", '"' + "\\" + '"'),
        ("A C", '"A?C"'),
        ("A\tB", '"A?B"'),
        ("A\rB", '"A?B"'),
        ("A\nB", '"A?B"'),
        ("A\x7fB", '"A?B"'),
        ("AéB", '"A?B"'),
        (" \t\r\n\x7fé", '"??????"'),
    ],
)
def test_show_quotes_printable_text_and_replaces_each_unsafe_character(
    text: str, expected: str
) -> None:
    assert show(text) == expected


@pytest.mark.parametrize(
    ("sequence", "expected"),
    [
        ("", ""),
        ("A", ""),
        ("AT", ""),
        ("ATG", "M"),
        ("ATGA", "M"),
        ("ATGAA", "M"),
        ("TTT", "F"),
        ("TTC", "F"),
        ("TTA", "L"),
        ("TTG", "L"),
        ("CCC", "P"),
        ("AAA", "K"),
        ("GGG", "G"),
        ("ATGNNNTAA", "MX*"),
        ("NAA", "X"),
        ("ANA", "X"),
        ("AAN", "X"),
    ],
)
def test_translate_complete_codons_and_n_containing_codons(sequence: str, expected: str) -> None:
    assert translate(sequence) == expected


@pytest.mark.parametrize(
    ("stop_codon", "expected"),
    [
        ("TAA", "M*"),
        ("TAG", "M*"),
        ("TGA", "M*"),
    ],
)
def test_translate_includes_each_stop_codon_by_default(stop_codon: str, expected: str) -> None:
    assert translate("ATG" + stop_codon) == expected


@pytest.mark.parametrize("stop_codon", ["TAA", "TAG", "TGA"])
def test_translate_to_stop_omits_each_stop_codon(stop_codon: str) -> None:
    assert translate("ATG" + stop_codon + "TGG", to_stop=True) == "M"


def test_translate_to_stop_returns_empty_when_the_first_codon_is_a_stop() -> None:
    assert translate("TAAATG", to_stop=True) == ""


def test_translate_without_to_stop_preserves_codons_after_a_stop() -> None:
    assert translate("ATGTAATGG") == "M*W"

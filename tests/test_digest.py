from inspect import signature
from io import BytesIO, StringIO
from pathlib import Path
from types import ModuleType

import pytest

from seqmill.core.options import OptionSpec, parse_options
from seqmill.core.seq import SeqError, show
from seqmill.core.text import render_table
from seqmill.formats.files import Io, check_stdin_once, load_records
from seqmill.tools import digest


class CountingStringIO(StringIO):
    def __init__(self) -> None:
        super().__init__()
        self.write_count = 0

    def write(self, value: str) -> int:
        self.write_count += 1
        return super().write(value)


def make_io(data: bytes = b"") -> Io:
    return Io(BytesIO(data), CountingStringIO(), StringIO())


def write_fasta(tmp_path: Path, name: str, text: str) -> str:
    path = tmp_path / name
    path.write_bytes(text.encode("utf-8"))
    return str(path)


def expected_table(rows: list[list[str]]) -> str:
    lines = render_table(["name", "fragment", "start", "end", "length"], rows, {1, 2, 3, 4})
    return "\n".join(lines) + "\n"


def test_digest_api_reuses_upstream_contracts_without_sibling_tool_imports() -> None:
    run_signature = signature(digest.run)
    assert list(run_signature.parameters) == ["argv", "io"]
    assert run_signature.parameters["argv"].annotation == list[str]
    assert run_signature.parameters["io"].annotation is Io
    assert run_signature.return_annotation is int

    assert digest.OptionSpec is OptionSpec
    assert digest.parse_options is parse_options
    assert digest.Io is Io
    assert digest.check_stdin_once is check_stdin_once
    assert digest.load_records is load_records
    assert digest.SeqError is SeqError
    assert digest.show is show
    assert digest.render_table is render_table

    for value in vars(digest).values():
        if isinstance(value, ModuleType):
            assert not value.__name__.startswith("seqmill.tools.")
        else:
            module_name = getattr(value, "__module__", "")
            assert not (
                module_name.startswith("seqmill.tools.") and module_name != digest.__name__
            )


def test_supplied_gene_example_renders_exact_table_with_one_lf_write(tmp_path: Path) -> None:
    path = write_fasta(
        tmp_path,
        "gene.fa",
        ">gene1 test construct\n"
        "GAATTCATGGCCAAAGGATCCTTTGGGAAATAGCCCGGGATGAAACCCTGAGTTTTCATTTGGGAAACAT\n",
    )
    streams = make_io()

    assert digest.run(["--enzyme", "EcoRI,BamHI", path], streams) == 0
    assert streams.stdout.getvalue() == (
        "name   fragment  start  end  length\n"
        "gene1         1      1    1       1\n"
        "gene1         2      2   16      15\n"
        "gene1         3     17   70      54\n"
    )
    assert streams.stdout.write_count == 1
    assert streams.stderr.getvalue() == ""


def test_supplied_circular_example_renders_exact_wrapped_fragment_table(tmp_path: Path) -> None:
    path = write_fasta(
        tmp_path,
        "ring.fa",
        ">ring\nGGATCCAAAAAAAAAAGAATTCTTTTTTTTTTGGATCCA\n",
    )
    streams = make_io()

    assert digest.run(["--enzyme=bamhi,ecori", "--circular", path], streams) == 0
    assert streams.stdout.getvalue() == (
        "name  fragment  start  end  length\n"
        "ring         1      2   17      16\n"
        "ring         2     18   33      16\n"
        "ring         3     34    1       7\n"
    )
    assert streams.stdout.write_count == 1
    assert streams.stderr.getvalue() == ""


@pytest.mark.parametrize(
    ("enzyme", "site", "offset"),
    [
        ("ApaI", "GGGCCC", 5),
        ("BamHI", "GGATCC", 1),
        ("EcoRI", "GAATTC", 1),
        ("EcoRV", "GATATC", 3),
        ("HindIII", "AAGCTT", 1),
        ("KpnI", "GGTACC", 5),
        ("NcoI", "CCATGG", 1),
        ("NotI", "GCGGCCGC", 2),
        ("PstI", "CTGCAG", 5),
        ("SacI", "GAGCTC", 5),
        ("SalI", "GTCGAC", 1),
        ("SmaI", "CCCGGG", 3),
        ("XbaI", "TCTAGA", 1),
        ("XhoI", "CTCGAG", 1),
    ],
)
def test_each_enzyme_recognizes_its_literal_site_and_exact_cut_offset(
    tmp_path: Path, enzyme: str, site: str, offset: int
) -> None:
    path = write_fasta(tmp_path, "site.fa", f">s\n{site}\n")
    streams = make_io()

    assert digest.run(["--enzyme", enzyme, path], streams) == 0
    assert streams.stdout.getvalue() == expected_table(
        [
            ["s", "1", "1", str(offset), str(offset)],
            ["s", "2", str(offset + 1), str(len(site)), str(len(site) - offset)],
        ]
    )
    assert streams.stdout.write_count == 1
    assert streams.stderr.getvalue() == ""


def test_missing_enzyme_precedes_missing_file() -> None:
    streams = make_io()

    with pytest.raises(SeqError) as caught:
        digest.run([], streams)

    assert caught.value.message == "missing --enzyme"
    assert streams.stdout.getvalue() == ""
    assert streams.stdout.write_count == 0
    assert streams.stderr.getvalue() == ""


def test_unknown_enzyme_is_reported_before_missing_file() -> None:
    streams = make_io()

    with pytest.raises(SeqError) as caught:
        digest.run(["--enzyme=bad"], streams)

    assert caught.value.message == (
        'unknown enzyme "bad" (expected ApaI, BamHI, EcoRI, EcoRV, HindIII, KpnI, '
        "NcoI, NotI, PstI, SacI, SalI, SmaI, XbaI or XhoI)"
    )
    assert streams.stdout.getvalue() == ""
    assert streams.stderr.getvalue() == ""


def test_empty_enzyme_item_has_exact_unknown_enzyme_diagnostic() -> None:
    streams = make_io()

    with pytest.raises(SeqError) as caught:
        digest.run(["--enzyme=EcoRI,,BamHI"], streams)

    assert caught.value.message == (
        'unknown enzyme "" (expected ApaI, BamHI, EcoRI, EcoRV, HindIII, KpnI, '
        "NcoI, NotI, PstI, SacI, SalI, SmaI, XbaI or XhoI)"
    )


def test_first_unknown_enzyme_wins_and_user_text_is_sanitized() -> None:
    streams = make_io()

    with pytest.raises(SeqError) as caught:
        digest.run(["--enzyme=bad,wrong"], streams)

    assert caught.value.message == (
        'unknown enzyme "bad" (expected ApaI, BamHI, EcoRI, EcoRV, HindIII, KpnI, '
        "NcoI, NotI, PstI, SacI, SalI, SmaI, XbaI or XhoI)"
    )

    with pytest.raises(SeqError) as unsafe:
        digest.run(["--enzyme=Bad\nName"], streams)
    assert unsafe.value.message.startswith('unknown enzyme "Bad?Name" ')


def test_unknown_option_syntax_precedes_enzyme_validation() -> None:
    streams = make_io()

    with pytest.raises(SeqError) as caught:
        digest.run(["--enzyme=bad", "--unknown"], streams)

    assert caught.value.message == 'unknown option "--unknown"'
    assert streams.stdout.getvalue() == ""
    assert streams.stderr.getvalue() == ""


def test_case_insensitive_repetition_is_deduplicated_and_last_value_wins(
    tmp_path: Path,
) -> None:
    path = write_fasta(tmp_path, "repeat.fa", ">s\nGAATTC\n")
    streams = make_io()

    assert digest.run(
        ["--enzyme=bad", "--enzyme=eCoRi,EcoRI,ECORI", path],
        streams,
    ) == 0
    assert streams.stdout.getvalue() == expected_table(
        [["s", "1", "1", "1", "1"], ["s", "2", "2", "6", "5"]]
    )
    assert streams.stdout.write_count == 1
    assert streams.stderr.getvalue() == ""


def test_repeated_enzyme_option_uses_only_the_last_value(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "last.fa", ">s\nGAATTC\n")
    streams = make_io()

    assert digest.run(["--enzyme=BamHI", "--enzyme=EcoRI", path], streams) == 0
    assert streams.stdout.getvalue() == expected_table(
        [["s", "1", "1", "1", "1"], ["s", "2", "2", "6", "5"]]
    )


def test_missing_file_follows_valid_own_checks() -> None:
    streams = make_io()

    with pytest.raises(SeqError) as caught:
        digest.run(["--enzyme=EcoRI"], streams)

    assert caught.value.message == "missing FILE"
    assert streams.stdout.getvalue() == ""
    assert streams.stderr.getvalue() == ""


def test_duplicate_stdin_is_rejected_before_any_file_is_read() -> None:
    streams = make_io(b">s\nGAATTC\n")

    with pytest.raises(SeqError) as caught:
        digest.run(["--enzyme=EcoRI", "-", "-", "missing.fa"], streams)

    assert caught.value.message == "standard input given twice"
    assert streams.stdin.tell() == 0
    assert streams.stdout.getvalue() == ""
    assert streams.stdout.write_count == 0
    assert streams.stderr.getvalue() == ""


def test_linear_ecori_fragments_have_inclusive_boundaries(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "linear.fa", ">s\nGAATTC\n")
    streams = make_io()

    assert digest.run(["--enzyme=EcoRI", path], streams) == 0
    assert streams.stdout.getvalue() == expected_table(
        [["s", "1", "1", "1", "1"], ["s", "2", "2", "6", "5"]]
    )


def test_circular_ecori_fragment_wraps_back_to_its_start(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "circle.fa", ">s\nGAATTC\n")
    streams = make_io()

    assert digest.run(["--enzyme=EcoRI", "--circular", path], streams) == 0
    assert streams.stdout.getvalue() == expected_table([["s", "1", "2", "1", "6"]])


def test_circular_cut_equal_to_sequence_length_is_retained(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "cut-n.fa", ">s\nAATTCG\n")
    streams = make_io()

    assert digest.run(["--enzyme=EcoRI", "--circular", path], streams) == 0
    assert streams.stdout.getvalue() == expected_table([["s", "1", "1", "6", "6"]])


def test_circular_cut_above_sequence_length_is_normalized(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "wrapped.fa", ">s\nCCGGGC\n")
    streams = make_io()

    assert digest.run(["--enzyme=ApaI", "--circular", path], streams) == 0
    assert streams.stdout.getvalue() == expected_table([["s", "1", "2", "1", "6"]])


def test_overlapping_sites_from_different_enzymes_both_contribute_cuts(
    tmp_path: Path,
) -> None:
    path = write_fasta(tmp_path, "overlap.fa", ">s\nGGGCCCGGG\n")
    streams = make_io()

    assert digest.run(["--enzyme=ApaI,SmaI", path], streams) == 0
    assert streams.stdout.getvalue() == expected_table(
        [["s", "1", "1", "5", "5"], ["s", "2", "6", "6", "1"], ["s", "3", "7", "9", "3"]]
    )


def test_shared_cut_coordinates_are_deduplicated(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "shared.fa", ">s\nGGGCCCATGG\n")
    streams = make_io()

    assert digest.run(["--enzyme=ApaI,NcoI", path], streams) == 0
    assert streams.stdout.getvalue() == expected_table(
        [["s", "1", "1", "5", "5"], ["s", "2", "6", "10", "5"]]
    )


def test_n_in_a_site_prevents_a_restriction_match(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "n.fa", ">s\nGAATNC\n")
    streams = make_io()

    assert digest.run(["--enzyme=EcoRI", path], streams) == 0
    assert streams.stdout.getvalue() == expected_table([["s", "1", "1", "6", "6"]])


@pytest.mark.parametrize("circular", [False, True])
def test_record_shorter_than_site_has_no_cut_on_a_line_or_ring(
    tmp_path: Path, circular: bool
) -> None:
    path = write_fasta(tmp_path, "short.fa", ">s\nGAATT\n")
    streams = make_io()
    args = ["--enzyme=EcoRI"]
    if circular:
        args.append("--circular")
    args.append(path)

    assert digest.run(args, streams) == 0
    assert streams.stdout.getvalue() == expected_table([["s", "1", "1", "5", "5"]])


def test_record_without_cuts_is_one_fragment(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "uncut.fa", ">s\nACGTN\n")
    streams = make_io()

    assert digest.run(["--enzyme=EcoRI", path], streams) == 0
    assert streams.stdout.getvalue() == expected_table([["s", "1", "1", "5", "5"]])


def test_empty_record_uses_the_specified_zero_length_fragment(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "empty-record.fa", ">empty\n")
    streams = make_io()

    assert digest.run(["--enzyme=EcoRI", path], streams) == 0
    assert streams.stdout.getvalue() == expected_table([["empty", "1", "1", "0", "0"]])


def test_empty_file_has_no_rows_and_returns_one(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "empty.fa", "")
    streams = make_io()

    assert digest.run(["--enzyme=EcoRI", path], streams) == 1
    assert streams.stdout.getvalue() == ""
    assert streams.stdout.write_count == 0
    assert streams.stderr.getvalue() == ""


def test_fragment_numbering_restarts_for_each_record(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "records.fa", ">first\nGAATTC\n>second\nGAATTC\n")
    streams = make_io()

    assert digest.run(["--enzyme=EcoRI", path], streams) == 0
    assert streams.stdout.getvalue() == expected_table(
        [
            ["first", "1", "1", "1", "1"],
            ["first", "2", "2", "6", "5"],
            ["second", "1", "1", "1", "1"],
            ["second", "2", "2", "6", "5"],
        ]
    )


def test_file_and_record_order_are_preserved(tmp_path: Path) -> None:
    first = write_fasta(tmp_path, "first.fa", ">first\nGAATTC\n>second\nGAATTC\n")
    second = write_fasta(tmp_path, "second.fa", ">third\nGAATTC\n")
    streams = make_io()

    assert digest.run(["--enzyme=EcoRI", first, second], streams) == 0
    assert streams.stdout.getvalue() == expected_table(
        [
            ["first", "1", "1", "1", "1"],
            ["first", "2", "2", "6", "5"],
            ["second", "1", "1", "1", "1"],
            ["second", "2", "2", "6", "5"],
            ["third", "1", "1", "1", "1"],
            ["third", "2", "2", "6", "5"],
        ]
    )
    assert streams.stdout.write_count == 1


@pytest.mark.parametrize("failure", ["unreadable", "malformed", "invalid_utf8"])
def test_late_file_failure_leaves_stdout_empty(
    tmp_path: Path, failure: str
) -> None:
    first = write_fasta(tmp_path, "first.fa", ">first\nGAATTC\n")
    bad = tmp_path / "bad.fa"
    if failure == "malformed":
        bad.write_bytes(b">bad\nACGR\n")
    elif failure == "invalid_utf8":
        bad.write_bytes(b">bad\nAC\xff\n")
    else:
        bad = tmp_path / "missing.fa"
    streams = make_io()

    with pytest.raises(SeqError) as caught:
        digest.run(["--enzyme=EcoRI", first, str(bad)], streams)

    display = bad.as_posix()
    if failure == "unreadable":
        expected = f"{display}: cannot read the file"
    elif failure == "malformed":
        expected = f'{display}:2: invalid character "R" in sequence'
    else:
        expected = f"{display}: input is not valid UTF-8"
    assert caught.value.message == expected
    assert streams.stdout.getvalue() == ""
    assert streams.stdout.write_count == 0
    assert streams.stderr.getvalue() == ""


def test_invalid_utf8_from_stdin_bytes_uses_the_binary_input_reader() -> None:
    streams = make_io(b">s\nAC\xff\n")

    with pytest.raises(SeqError) as caught:
        digest.run(["--enzyme=EcoRI", "-"], streams)

    assert caught.value.message == "<stdin>: input is not valid UTF-8"
    assert streams.stdout.getvalue() == ""
    assert streams.stdout.write_count == 0
    assert streams.stderr.getvalue() == ""


def test_valid_stdin_is_digested_through_the_binary_input_reader() -> None:
    streams = make_io(b">s\nGAATTC\n")

    assert digest.run(["--enzyme=EcoRI", "-"], streams) == 0
    assert streams.stdout.getvalue() == expected_table(
        [["s", "1", "1", "1", "1"], ["s", "2", "2", "6", "5"]]
    )
    assert streams.stdout.write_count == 1
    assert streams.stderr.getvalue() == ""

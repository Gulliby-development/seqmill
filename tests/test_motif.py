import ast
from inspect import signature
from io import BytesIO, StringIO
from pathlib import Path
from types import ModuleType

import pytest

from seqmill.core.options import OptionSpec, parse_options
from seqmill.core.seq import SeqError, show
from seqmill.core.text import render_table
from seqmill.formats.files import Io, check_stdin_once, load_records
from seqmill.tools import motif


class CountingStringIO(StringIO):
    def __init__(self) -> None:
        super().__init__()
        self.write_count = 0

    def write(self, value: str) -> int:
        self.write_count += 1
        return super().write(value)


def make_io(data: bytes = b"") -> Io:
    return Io(BytesIO(data), CountingStringIO(), StringIO())


def write_fasta(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_bytes(text.encode("utf-8"))
    return path


def expected_table(rows: list[list[str]]) -> str:
    lines = render_table(
        ["name", "strand", "start", "end", "mismatches", "match"],
        rows,
        {2, 3, 4},
    )
    return "\n".join(lines) + "\n"


def test_motif_api_reuses_upstream_contracts_without_sibling_tool_imports() -> None:
    run_signature = signature(motif.run)
    assert list(run_signature.parameters) == ["argv", "io"]
    assert run_signature.parameters["argv"].annotation == list[str]
    assert run_signature.parameters["io"].annotation is Io
    assert run_signature.return_annotation is int

    assert motif.OptionSpec is OptionSpec
    assert motif.parse_options is parse_options
    assert motif.Io is Io
    assert motif.check_stdin_once is check_stdin_once
    assert motif.load_records is load_records
    assert motif.SeqError is SeqError
    assert motif.show is show
    assert motif.render_table is render_table

    for value in vars(motif).values():
        if isinstance(value, ModuleType):
            assert not value.__name__.startswith("seqmill.tools.")
        else:
            module_name = getattr(value, "__module__", "")
            assert not (
                module_name.startswith("seqmill.tools.") and module_name != motif.__name__
            )

    source = Path(motif.__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    prohibited_imports = {"locale", "os", "random", "subprocess", "time"}
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            assert all(alias.name.split(".")[0] not in prohibited_imports for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            assert node.module.split(".")[0] not in prohibited_imports


def test_supplied_gene_example_renders_exact_table_and_writes_once(tmp_path: Path) -> None:
    path = write_fasta(
        tmp_path,
        "gene.fa",
        ">gene1 test construct\n"
        "GAATTCATGGCCAAAGGATCCTTTGGGAAATAGCCCGGGATGAAACCCTGAGTTTTCATTTGGGAAACAT\n",
    )
    streams = make_io()

    assert motif.run(["--both-strands", "--mismatches", "1", "ATGAAR", str(path)], streams) == 0
    assert streams.stdout.getvalue() == (
        "name   strand  start  end  mismatches  match\n"
        "gene1  -           3    8           1  ATTCAT\n"
        "gene1  +          40   45           0  ATGAAA\n"
        "gene1  -          54   59           0  TTTCAT\n"
    )
    assert streams.stdout.write_count == 1
    assert streams.stderr.getvalue() == ""


@pytest.mark.parametrize(
    ("pattern", "rows"),
    [
        ("A", [["s", "+", "1", "1", "0", "A"]]),
        ("C", [["s", "+", "2", "2", "0", "C"]]),
        ("G", [["s", "+", "3", "3", "0", "G"]]),
        ("T", [["s", "+", "4", "4", "0", "T"]]),
        ("R", [["s", "+", "1", "1", "0", "A"], ["s", "+", "3", "3", "0", "G"]]),
        ("Y", [["s", "+", "2", "2", "0", "C"], ["s", "+", "4", "4", "0", "T"]]),
        ("S", [["s", "+", "2", "2", "0", "C"], ["s", "+", "3", "3", "0", "G"]]),
        ("W", [["s", "+", "1", "1", "0", "A"], ["s", "+", "4", "4", "0", "T"]]),
        ("K", [["s", "+", "3", "3", "0", "G"], ["s", "+", "4", "4", "0", "T"]]),
        ("M", [["s", "+", "1", "1", "0", "A"], ["s", "+", "2", "2", "0", "C"]]),
        (
            "B",
            [["s", "+", "2", "2", "0", "C"], ["s", "+", "3", "3", "0", "G"], ["s", "+", "4", "4", "0", "T"]],
        ),
        (
            "D",
            [["s", "+", "1", "1", "0", "A"], ["s", "+", "3", "3", "0", "G"], ["s", "+", "4", "4", "0", "T"]],
        ),
        (
            "H",
            [["s", "+", "1", "1", "0", "A"], ["s", "+", "2", "2", "0", "C"], ["s", "+", "4", "4", "0", "T"]],
        ),
        (
            "V",
            [["s", "+", "1", "1", "0", "A"], ["s", "+", "2", "2", "0", "C"], ["s", "+", "3", "3", "0", "G"]],
        ),
        (
            "N",
            [
                ["s", "+", "1", "1", "0", "A"],
                ["s", "+", "2", "2", "0", "C"],
                ["s", "+", "3", "3", "0", "G"],
                ["s", "+", "4", "4", "0", "T"],
                ["s", "+", "5", "5", "0", "N"],
            ],
        ),
    ],
)
def test_each_iupac_symbol_matches_its_declared_set_and_n_only_matches_n(
    tmp_path: Path, pattern: str, rows: list[list[str]]
) -> None:
    path = write_fasta(tmp_path, "bases.fa", ">s\nACGTN\n")
    streams = make_io()

    assert motif.run([pattern, str(path)], streams) == 0
    assert streams.stdout.getvalue() == expected_table(rows)
    assert streams.stdout.write_count == 1
    assert streams.stderr.getvalue() == ""


@pytest.mark.parametrize(
    ("pattern", "sequence", "rows"),
    [
        ("A", "T", [["s", "-", "1", "1", "0", "T"]]),
        ("C", "G", [["s", "-", "1", "1", "0", "G"]]),
        ("G", "C", [["s", "-", "1", "1", "0", "C"]]),
        ("T", "A", [["s", "-", "1", "1", "0", "A"]]),
        ("R", "C", [["s", "-", "1", "1", "0", "C"]]),
        ("Y", "A", [["s", "-", "1", "1", "0", "A"]]),
        ("S", "C", [["s", "+", "1", "1", "0", "C"], ["s", "-", "1", "1", "0", "C"]]),
        ("W", "A", [["s", "+", "1", "1", "0", "A"], ["s", "-", "1", "1", "0", "A"]]),
        ("K", "A", [["s", "-", "1", "1", "0", "A"]]),
        ("M", "G", [["s", "-", "1", "1", "0", "G"]]),
        ("B", "A", [["s", "-", "1", "1", "0", "A"]]),
        ("V", "T", [["s", "-", "1", "1", "0", "T"]]),
        ("D", "C", [["s", "-", "1", "1", "0", "C"]]),
        ("H", "G", [["s", "-", "1", "1", "0", "G"]]),
        ("N", "N", [["s", "+", "1", "1", "0", "N"], ["s", "-", "1", "1", "0", "N"]]),
    ],
)
def test_reverse_strand_uses_every_iupac_complement_pair(
    tmp_path: Path, pattern: str, sequence: str, rows: list[list[str]]
) -> None:
    path = write_fasta(tmp_path, "complement.fa", f">s\n{sequence}\n")
    streams = make_io()

    assert motif.run(["--both-strands", pattern, str(path)], streams) == 0
    assert streams.stdout.getvalue() == expected_table(rows)
    assert streams.stderr.getvalue() == ""


def test_overlapping_windows_are_all_reported(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "overlap.fa", ">s\nAAA\n")
    streams = make_io()

    assert motif.run(["AA", str(path)], streams) == 0
    assert streams.stdout.getvalue() == expected_table(
        [["s", "+", "1", "2", "0", "AA"], ["s", "+", "2", "3", "0", "AA"]]
    )
    assert streams.stdout.write_count == 1


@pytest.mark.parametrize(
    ("mismatches", "sequence", "expected_rows", "expected_status"),
    [
        (0, "ATG", [["s", "+", "1", "3", "0", "ATG"]], 0),
        (1, "ACG", [["s", "+", "1", "3", "1", "ACG"]], 0),
        (1, "CCG", [], 1),
        (2, "CCG", [["s", "+", "1", "3", "2", "CCG"]], 0),
    ],
)
def test_mismatch_limit_includes_the_boundary_and_rejects_one_above(
    tmp_path: Path,
    mismatches: int,
    sequence: str,
    expected_rows: list[list[str]],
    expected_status: int,
) -> None:
    path = write_fasta(tmp_path, "mismatches.fa", f">s\n{sequence}\n")
    streams = make_io()

    assert motif.run([f"--mismatches={mismatches}", "ATG", str(path)], streams) == expected_status
    assert streams.stdout.getvalue() == (
        expected_table(expected_rows) if expected_rows else ""
    )
    assert streams.stdout.write_count == (1 if expected_rows else 0)
    assert streams.stderr.getvalue() == ""


def test_mismatch_eight_is_accepted_at_its_maximum_boundary(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "eight.fa", ">s\nCCCCCCCCA\n")
    streams = make_io()

    assert motif.run(["--mismatches=8", "AAAAAAAAA", str(path)], streams) == 0
    assert streams.stdout.getvalue() == expected_table(
        [["s", "+", "1", "9", "8", "CCCCCCCCA"]]
    )


def test_lowercase_pattern_is_normalized_for_search(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "lowercase.fa", ">s\natgaaa\n")
    streams = make_io()

    assert motif.run(["atg", str(path)], streams) == 0
    assert streams.stdout.getvalue() == expected_table(
        [["s", "+", "1", "3", "0", "ATG"]]
    )


def test_pattern_length_sixty_four_is_valid(tmp_path: Path) -> None:
    pattern = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    path = write_fasta(
        tmp_path,
        "sixty-four.fa",
        ">s\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n",
    )
    streams = make_io()

    assert motif.run([pattern, str(path)], streams) == 0
    assert streams.stdout.getvalue() == expected_table(
        [["s", "+", "1", "64", "0", pattern]]
    )


@pytest.mark.parametrize(
    ("pattern", "expected"),
    [
        ("", 'invalid pattern "" (expected 1 to 64 letters from ACGTNRYSWKMBDHV)'),
        (
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            (
                'invalid pattern "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" '
                "(expected 1 to 64 letters from ACGTNRYSWKMBDHV)"
            ),
        ),
        ("ATGX", 'invalid pattern "ATGX" (expected 1 to 64 letters from ACGTNRYSWKMBDHV)'),
        ("ATGé", 'invalid pattern "ATG?" (expected 1 to 64 letters from ACGTNRYSWKMBDHV)'),
        ("ATG\n", 'invalid pattern "ATG?" (expected 1 to 64 letters from ACGTNRYSWKMBDHV)'),
    ],
)
def test_invalid_pattern_has_exact_diagnostic_before_file_checks(
    pattern: str, expected: str
) -> None:
    streams = make_io()

    with pytest.raises(SeqError) as caught:
        motif.run([pattern], streams)

    assert caught.value.message == expected
    assert streams.stdout.getvalue() == ""
    assert streams.stdout.write_count == 0
    assert streams.stderr.getvalue() == ""


def test_pattern_longer_than_sixty_four_is_rejected(tmp_path: Path) -> None:
    streams = make_io()
    path = tmp_path / "not-read.fa"

    with pytest.raises(SeqError) as caught:
        motif.run(
            [
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                str(path),
            ],
            streams,
        )

    assert caught.value.message == (
        'invalid pattern "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" '
        "(expected 1 to 64 letters from ACGTNRYSWKMBDHV)"
    )
    assert not path.exists()
    assert streams.stdout.getvalue() == ""


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("-1", 'invalid --mismatches "-1" (expected an integer from 0 to 8)'),
        ("9", 'invalid --mismatches "9" (expected an integer from 0 to 8)'),
    ],
)
def test_invalid_mismatch_bounds_precede_pattern_and_file_validation(
    value: str, expected: str
) -> None:
    streams = make_io()

    with pytest.raises(SeqError) as caught:
        motif.run([f"--mismatches={value}", "XYZ"], streams)

    assert caught.value.message == expected
    assert streams.stdout.getvalue() == ""
    assert streams.stderr.getvalue() == ""


def test_missing_pattern_is_reported_before_missing_file() -> None:
    streams = make_io()

    with pytest.raises(SeqError) as caught:
        motif.run([], streams)

    assert caught.value.message == "missing PATTERN"
    assert streams.stdout.getvalue() == ""
    assert streams.stderr.getvalue() == ""


def test_syntax_error_is_left_to_right_and_precedes_deferred_value_validation() -> None:
    streams = make_io()

    with pytest.raises(SeqError) as caught:
        motif.run(["--mismatches=9", "--unknown"], streams)

    assert caught.value.message == 'unknown option "--unknown"'
    assert streams.stdout.getvalue() == ""
    assert streams.stdout.write_count == 0
    assert streams.stderr.getvalue() == ""


def test_missing_option_value_has_exact_syntax_diagnostic() -> None:
    streams = make_io()

    with pytest.raises(SeqError) as caught:
        motif.run(["--mismatches"], streams)

    assert caught.value.message == "--mismatches needs a value"
    assert streams.stdout.getvalue() == ""
    assert streams.stderr.getvalue() == ""


def test_last_mismatch_option_value_wins(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "last-wins.fa", ">s\nACG\n")
    streams = make_io()

    assert motif.run(
        ["--mismatches=9", "ATG", "--mismatches", "1", str(path)],
        streams,
    ) == 0
    assert streams.stdout.getvalue() == expected_table(
        [["s", "+", "1", "3", "1", "ACG"]]
    )


def test_missing_file_follows_valid_pattern_checks() -> None:
    streams = make_io()

    with pytest.raises(SeqError) as caught:
        motif.run(["ATG"], streams)

    assert caught.value.message == "missing FILE"
    assert streams.stdout.getvalue() == ""
    assert streams.stderr.getvalue() == ""


def test_invalid_pattern_precedes_duplicate_stdin_check() -> None:
    streams = make_io(b">s\nATG\n")

    with pytest.raises(SeqError) as caught:
        motif.run(["XYZ", "-", "-"], streams)

    assert caught.value.message == (
        'invalid pattern "XYZ" (expected 1 to 64 letters from ACGTNRYSWKMBDHV)'
    )
    assert streams.stdin.tell() == 0
    assert streams.stdout.getvalue() == ""
    assert streams.stderr.getvalue() == ""


def test_duplicate_stdin_is_rejected_before_any_input_is_read() -> None:
    streams = make_io(b">s\nATG\n")

    with pytest.raises(SeqError) as caught:
        motif.run(["ATG", "-", "-", "missing.fa"], streams)

    assert caught.value.message == "standard input given twice"
    assert streams.stdin.tell() == 0
    assert streams.stdout.getvalue() == ""
    assert streams.stdout.write_count == 0
    assert streams.stderr.getvalue() == ""


def test_empty_file_has_no_rows_and_returns_one(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "empty.fa", "")
    streams = make_io()

    assert motif.run(["ATG", str(path)], streams) == 1
    assert streams.stdout.getvalue() == ""
    assert streams.stdout.write_count == 0
    assert streams.stderr.getvalue() == ""


def test_empty_sequence_and_shorter_record_have_no_rows(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "empty-record.fa", ">empty\n")
    streams = make_io()

    assert motif.run(["ATG", str(path)], streams) == 1
    assert streams.stdout.getvalue() == ""
    assert streams.stdout.write_count == 0
    assert streams.stderr.getvalue() == ""


def test_match_coordinates_and_window_are_on_the_original_reverse_strand_record(
    tmp_path: Path,
) -> None:
    path = write_fasta(tmp_path, "reverse.fa", ">s\nCCAT\n")
    streams = make_io()

    assert motif.run(["--both-strands", "ATG", str(path)], streams) == 0
    assert streams.stdout.getvalue() == expected_table(
        [["s", "-", "2", "4", "0", "CAT"]]
    )


def test_palindrome_reports_both_strands_at_the_same_coordinates(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "palindrome.fa", ">s\nAT\n")
    streams = make_io()

    assert motif.run(["--both-strands", "AT", str(path)], streams) == 0
    assert streams.stdout.getvalue() == expected_table(
        [
            ["s", "+", "1", "2", "0", "AT"],
            ["s", "-", "1", "2", "0", "AT"],
        ]
    )


def test_same_start_sorts_plus_before_minus_even_with_more_mismatches(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "strand-order.fa", ">s\nATT\n")
    streams = make_io()

    assert motif.run(
        ["--both-strands", "--mismatches=2", "AAG", str(path)],
        streams,
    ) == 0
    assert streams.stdout.getvalue() == expected_table(
        [
            ["s", "+", "1", "3", "2", "ATT"],
            ["s", "-", "1", "3", "1", "ATT"],
        ]
    )


def test_file_and_record_order_are_preserved(tmp_path: Path) -> None:
    first = write_fasta(tmp_path, "first.fa", ">first\nAAA\n>second\nAAA\n")
    second = write_fasta(tmp_path, "second.fa", ">third\nAAA\n")
    streams = make_io()

    assert motif.run(["AA", str(first), str(second)], streams) == 0
    assert streams.stdout.getvalue() == expected_table(
        [
            ["first", "+", "1", "2", "0", "AA"],
            ["first", "+", "2", "3", "0", "AA"],
            ["second", "+", "1", "2", "0", "AA"],
            ["second", "+", "2", "3", "0", "AA"],
            ["third", "+", "1", "2", "0", "AA"],
            ["third", "+", "2", "3", "0", "AA"],
        ]
    )
    assert streams.stdout.write_count == 1
    assert streams.stderr.getvalue() == ""


@pytest.mark.parametrize("failure", ["unreadable", "malformed", "invalid_utf8"])
def test_late_file_failure_leaves_stdout_empty(tmp_path: Path, failure: str) -> None:
    first = write_fasta(tmp_path, "first.fa", ">first\nATG\n")
    bad = tmp_path / "bad.fa"
    if failure == "malformed":
        bad.write_bytes(b">bad\nATGR\n")
    elif failure == "invalid_utf8":
        bad.write_bytes(b">bad\nAC\xff\n")
    else:
        bad = tmp_path / "missing.fa"
    streams = make_io()

    with pytest.raises(SeqError) as caught:
        motif.run(["ATG", str(first), str(bad)], streams)

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


@pytest.mark.parametrize("source", ["stdin", "file"])
def test_invalid_utf8_bytes_are_decoded_from_the_real_input_channel(
    tmp_path: Path, source: str
) -> None:
    raw = b">s\nAC\xff\n"
    streams = make_io(raw if source == "stdin" else b"")
    path = "-"
    if source == "file":
        file_path = tmp_path / "invalid.fa"
        file_path.write_bytes(raw)
        path = str(file_path)
        expected_display = file_path.as_posix()
    else:
        expected_display = "<stdin>"

    with pytest.raises(SeqError) as caught:
        motif.run(["AC", path], streams)

    assert caught.value.message == f"{expected_display}: input is not valid UTF-8"
    assert streams.stdout.getvalue() == ""
    assert streams.stdout.write_count == 0
    assert streams.stderr.getvalue() == ""


def test_valid_stdin_is_searched_through_the_binary_input_reader() -> None:
    streams = make_io(b">s\nATG\n")

    assert motif.run(["ATG", "-"], streams) == 0
    assert streams.stdout.getvalue() == expected_table(
        [["s", "+", "1", "3", "0", "ATG"]]
    )
    assert streams.stdout.write_count == 1
    assert streams.stderr.getvalue() == ""

from inspect import signature
from io import BytesIO, StringIO
from pathlib import Path
from types import ModuleType

import pytest

from seqmill.core.options import OptionSpec, parse_options
from seqmill.core.seq import SeqError, reverse_complement
from seqmill.core.text import render_table
from seqmill.formats.files import Io, check_stdin_once, load_records
from seqmill.tools import orfs


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
    lines = render_table(["name", "frame", "start", "end", "aa"], rows, {1, 2, 3, 4})
    return "\n".join(lines) + "\n"


def test_orfs_api_reuses_upstream_contracts_without_sibling_tool_imports() -> None:
    run_signature = signature(orfs.run)
    assert list(run_signature.parameters) == ["argv", "io"]
    assert run_signature.parameters["argv"].annotation == list[str]
    assert run_signature.parameters["io"].annotation is Io
    assert run_signature.return_annotation is int

    assert orfs.OptionSpec is OptionSpec
    assert orfs.parse_options is parse_options
    assert orfs.Io is Io
    assert orfs.check_stdin_once is check_stdin_once
    assert orfs.load_records is load_records
    assert orfs.SeqError is SeqError
    assert orfs.reverse_complement is reverse_complement
    assert orfs.render_table is render_table

    for value in vars(orfs).values():
        if isinstance(value, ModuleType):
            assert not value.__name__.startswith("seqmill.tools.")
        else:
            module_name = getattr(value, "__module__", "")
            assert not (
                module_name.startswith("seqmill.tools.") and module_name != orfs.__name__
            )


def test_supplied_gene_example_renders_exact_table_with_one_lf_write(tmp_path: Path) -> None:
    path = write_fasta(
        tmp_path,
        "gene.fa",
        ">gene1 test construct\n"
        "GAATTCATGGCCAAAGGATCCTTTGGGAAATAGCCCGGGATGAAACCCTGAGTTTTCATTTGGGAAACAT\n",
    )
    streams = make_io()

    assert orfs.run(["--min", "2", path], streams) == 0
    assert streams.stdout.getvalue() == (
        "name   frame  start  end  aa\n"
        "gene1     +1      7   33   8\n"
        "gene1     +1     40   51   3\n"
        "gene1     -1     56   70   4\n"
    )
    assert streams.stdout.write_count == 1
    assert streams.stderr.getvalue() == ""


def test_supplied_gene_example_uses_default_minimum_of_thirty(tmp_path: Path) -> None:
    path = write_fasta(
        tmp_path,
        "gene.fa",
        ">gene1 test construct\n"
        "GAATTCATGGCCAAAGGATCCTTTGGGAAATAGCCCGGGATGAAACCCTGAGTTTTCATTTGGGAAACAT\n",
    )
    streams = make_io()

    assert orfs.run([path], streams) == 1
    assert streams.stdout.getvalue() == ""
    assert streams.stdout.write_count == 0
    assert streams.stderr.getvalue() == ""


@pytest.mark.parametrize(
    ("sequence", "frame", "start", "end"),
    [
        ("ATGAAATAA", 1, 1, 9),
        ("CATGAAATAA", 2, 2, 10),
        ("CCATGAAATAA", 3, 3, 11),
        ("TTATTTCAT", -1, 1, 9),
        ("TTATTTCATG", -2, 1, 9),
        ("TTATTTCATGG", -3, 1, 9),
    ],
)
def test_scans_each_of_the_six_frames_and_maps_reverse_coordinates(
    tmp_path: Path, sequence: str, frame: int, start: int, end: int
) -> None:
    path = write_fasta(tmp_path, "frame.fa", f">s\n{sequence}\n")
    streams = make_io()

    assert orfs.run(["--min=1", path], streams) == 0
    assert streams.stdout.getvalue() == expected_table(
        [["s", f"{frame:+d}", str(start), str(end), "2"]]
    )
    assert streams.stdout.write_count == 1
    assert streams.stderr.getvalue() == ""


@pytest.mark.parametrize("stop", ["TAA", "TAG", "TGA"])
def test_each_stop_codon_closes_an_orf_inclusively(tmp_path: Path, stop: str) -> None:
    path = write_fasta(tmp_path, "stop.fa", f">s\nATG{stop}\n")
    streams = make_io()

    assert orfs.run(["--min", "1", path], streams) == 0
    assert streams.stdout.getvalue() == expected_table([["s", "+1", "1", "6", "1"]])


def test_nested_start_does_not_replace_first_start(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "nested.fa", ">s\nATGATGAAATAA\n")
    streams = make_io()

    assert orfs.run(["--min", "1", path], streams) == 0
    assert streams.stdout.getvalue() == expected_table([["s", "+1", "1", "12", "3"]])


def test_scanning_resumes_after_a_stop_and_keeps_nonoverlapping_orfs(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "resumed.fa", ">s\nATGAAATAAATGCCCTAG\n")
    streams = make_io()

    assert orfs.run(["--min", "2", path], streams) == 0
    assert streams.stdout.getvalue() == expected_table(
        [["s", "+1", "1", "9", "2"], ["s", "+1", "10", "18", "2"]]
    )


def test_orfs_from_multiple_frames_sort_by_original_start_then_end(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "sorted.fa", ">s\nATGAAATAAAATGAAATAA\n")
    streams = make_io()

    assert orfs.run(["--min", "1", path], streams) == 0
    assert streams.stdout.getvalue() == expected_table(
        [["s", "+1", "1", "9", "2"], ["s", "+2", "11", "19", "2"]]
    )


def test_n_codon_counts_toward_amino_acid_length_but_is_not_a_start_or_stop(
    tmp_path: Path,
) -> None:
    path = write_fasta(tmp_path, "n.fa", ">s\nATGNNNTAA\n")
    streams = make_io()

    assert orfs.run(["--min", "2", path], streams) == 0
    assert streams.stdout.getvalue() == expected_table([["s", "+1", "1", "9", "2"]])


@pytest.mark.parametrize("sequence", ["ATGAAAT", "ATGAAA", ""])
def test_partial_codon_unterminated_start_and_empty_sequence_have_no_orf(
    tmp_path: Path, sequence: str
) -> None:
    path = write_fasta(tmp_path, "none.fa", f">s\n{sequence}\n")
    streams = make_io()

    assert orfs.run(["--min", "1", path], streams) == 1
    assert streams.stdout.getvalue() == ""
    assert streams.stdout.write_count == 0
    assert streams.stderr.getvalue() == ""


def test_empty_fasta_file_has_no_rows(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "empty.fa", "")
    streams = make_io()

    assert orfs.run([path], streams) == 1
    assert streams.stdout.getvalue() == ""
    assert streams.stdout.write_count == 0


def test_minimum_is_inclusive_at_the_orf_amino_acid_threshold(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "threshold.fa", ">s\nATGAAATAA\n")
    streams = make_io()

    assert orfs.run(["--min=2", path], streams) == 0
    assert streams.stdout.getvalue() == expected_table([["s", "+1", "1", "9", "2"]])


def test_minimum_above_threshold_filters_the_orf(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "threshold.fa", ">s\nATGAAATAA\n")
    streams = make_io()

    assert orfs.run(["--min=3", path], streams) == 1
    assert streams.stdout.getvalue() == ""


@pytest.mark.parametrize(
    ("raw_minimum", "expected"),
    [
        ("0", 'invalid --min "0" (expected an integer from 1 to 100000)'),
        ("100001", 'invalid --min "100001" (expected an integer from 1 to 100000)'),
    ],
)
def test_minimum_rejects_values_beyond_inclusive_bounds(
    raw_minimum: str, expected: str
) -> None:
    streams = make_io()

    with pytest.raises(SeqError) as caught:
        orfs.run([f"--min={raw_minimum}"], streams)

    assert caught.value.message == expected
    assert streams.stdout.getvalue() == ""
    assert streams.stdout.write_count == 0
    assert streams.stderr.getvalue() == ""


def test_minimum_accepts_upper_bound_even_when_no_orf_is_long_enough(
    tmp_path: Path,
) -> None:
    path = write_fasta(tmp_path, "short.fa", ">s\nATGAAATAA\n")
    streams = make_io()

    assert orfs.run(["--min=100000", path], streams) == 1
    assert streams.stdout.getvalue() == ""


def test_option_syntax_error_precedes_deferred_minimum_validation() -> None:
    streams = make_io()

    with pytest.raises(SeqError) as caught:
        orfs.run(["--min=0", "--unknown"], streams)

    assert caught.value.message == 'unknown option "--unknown"'
    assert streams.stdout.getvalue() == ""
    assert streams.stderr.getvalue() == ""


def test_last_minimum_value_wins_even_when_earlier_value_is_invalid(
    tmp_path: Path,
) -> None:
    path = write_fasta(tmp_path, "last-value.fa", ">s\nATGAAATAA\n")
    streams = make_io()

    assert orfs.run(["--min=0", "--min=2", path], streams) == 0
    assert streams.stdout.getvalue() == expected_table([["s", "+1", "1", "9", "2"]])


def test_missing_file_is_reported_after_options_are_validated() -> None:
    streams = make_io()

    with pytest.raises(SeqError) as caught:
        orfs.run([], streams)

    assert caught.value.message == "missing FILE"
    assert streams.stdout.getvalue() == ""
    assert streams.stderr.getvalue() == ""


def test_minimum_error_precedes_missing_file_check() -> None:
    streams = make_io()

    with pytest.raises(SeqError) as caught:
        orfs.run(["--min=0"], streams)

    assert caught.value.message == 'invalid --min "0" (expected an integer from 1 to 100000)'


def test_duplicate_stdin_is_rejected_before_any_file_is_read() -> None:
    streams = make_io(b">s\nATGAAATAA\n")

    with pytest.raises(SeqError) as caught:
        orfs.run(["--min=1", "-", "-"], streams)

    assert caught.value.message == "standard input given twice"
    assert streams.stdin.tell() == 0
    assert streams.stdout.getvalue() == ""
    assert streams.stdout.write_count == 0
    assert streams.stderr.getvalue() == ""


def test_files_and_records_keep_input_order(tmp_path: Path) -> None:
    first = write_fasta(tmp_path, "first.fa", ">first\nATGAAATAA\n>empty\n\n")
    second = write_fasta(tmp_path, "second.fa", ">last\nATGCCCTAG\n")
    streams = make_io()

    assert orfs.run(["--min=1", first, second], streams) == 0
    assert streams.stdout.getvalue() == expected_table(
        [["first", "+1", "1", "9", "2"], ["last", "+1", "1", "9", "2"]]
    )
    assert streams.stdout.write_count == 1


@pytest.mark.parametrize("failure", ["unreadable", "malformed", "invalid_utf8"])
def test_late_file_failure_leaves_stdout_empty(
    tmp_path: Path, failure: str
) -> None:
    first = write_fasta(tmp_path, "first.fa", ">first\nATGAAATAA\n")
    bad = tmp_path / "bad.fa"
    if failure == "malformed":
        bad.write_bytes(b">bad\nACGR\n")
    elif failure == "invalid_utf8":
        bad.write_bytes(b">bad\nAC\xff\n")
    else:
        bad = tmp_path / "missing.fa"
    streams = make_io()

    with pytest.raises(SeqError) as caught:
        orfs.run(["--min=1", first, str(bad)], streams)

    display = str(bad).replace("\\", "/")
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


def test_invalid_utf8_from_stdin_bytes_reaches_the_real_binary_reader() -> None:
    streams = make_io(b">s\nAC\xff\n")

    with pytest.raises(SeqError) as caught:
        orfs.run(["--min=1", "-"], streams)

    assert caught.value.message == "<stdin>: input is not valid UTF-8"
    assert streams.stdout.getvalue() == ""
    assert streams.stdout.write_count == 0
    assert streams.stderr.getvalue() == ""


def test_valid_stdin_is_scanned_through_the_binary_reader() -> None:
    streams = make_io(b">s\nATGAAATAA\n")

    assert orfs.run(["--min=1", "-"], streams) == 0
    assert streams.stdout.getvalue() == expected_table([["s", "+1", "1", "9", "2"]])
    assert streams.stdout.write_count == 1
    assert streams.stderr.getvalue() == ""

from inspect import signature
from io import BytesIO, StringIO
from pathlib import Path

import pytest

from seqmill.core.seq import SeqError
from seqmill.formats.files import Io
from seqmill.tools import convert, stats


class RecordingStream(StringIO):
    def __init__(self) -> None:
        super().__init__()
        self.writes: list[str] = []

    def write(self, value: str) -> int:
        self.writes.append(value)
        return super().write(value)


def make_io(data: bytes = b"") -> tuple[Io, RecordingStream, RecordingStream]:
    stdout = RecordingStream()
    stderr = RecordingStream()
    return Io(BytesIO(data), stdout, stderr), stdout, stderr


def write_fixture(tmp_path: Path, name: str, data: bytes) -> str:
    (tmp_path / name).write_bytes(data)
    return name


def assert_failure(io: Io, argv: list[str], message: str) -> None:
    with pytest.raises(SeqError) as caught:
        stats.run(argv, io)

    assert caught.value.message == message
    assert io.stdout.getvalue() == ""
    assert io.stderr.getvalue() == ""
    assert io.stdout.writes == []
    assert io.stderr.writes == []


def test_stats_api_and_tools_package_marker_are_exact() -> None:
    run_signature = signature(stats.run)
    assert list(run_signature.parameters) == ["argv", "io"]
    assert run_signature.parameters["argv"].annotation == list[str]
    assert run_signature.parameters["io"].annotation is Io
    assert run_signature.return_annotation is int

    marker = Path(__file__).resolve().parents[1] / "seqmill" / "tools" / "__init__.py"
    assert marker.read_bytes() == b""


def test_stats_renders_the_supplied_reads_example_exactly(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(
        tmp_path,
        "reads.fa",
        b">seq1 first read\nACGTACGTAC\nGTNNAC\n>seq2\nGGGCCC\n>seq3 empty\n",
    )
    io, stdout, stderr = make_io()

    assert stats.run(["reads.fa"], io) == 0

    assert stdout.getvalue() == (
        "file      records  bases  min  max  mean  n50   gc%  n\n"
        "reads.fa        3     22    0   16   7.3   16  65.0  2\n"
    )
    assert stdout.writes == [stdout.getvalue()]
    assert stderr.getvalue() == ""
    assert stderr.writes == []


def test_stats_empty_file_has_header_and_zero_row(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, "empty.fa", b"")
    io, stdout, stderr = make_io()

    assert stats.run(["empty.fa"], io) == 0

    assert stdout.getvalue() == (
        "file      records  bases  min  max  mean  n50  gc%  n\n"
        "empty.fa        0      0    -    -     -    -    -  0\n"
    )
    assert stdout.writes == [stdout.getvalue()]
    assert stderr.getvalue() == ""


@pytest.mark.parametrize(
    ("data", "expected_cells"),
    [
        (b">empty\n", ["empty.fa", "1", "0", "0", "0", "0.0", "-", "-", "0"]),
        (b">n-only\nNNN\n", ["empty.fa", "1", "3", "3", "3", "3.0", "3", "-", "3"]),
        (
            b">canonical-and-n\nACGTNN\n",
            ["empty.fa", "1", "6", "6", "6", "6.0", "6", "50.0", "2"],
        ),
    ],
)
def test_stats_empty_record_and_base_composition_boundaries(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    data: bytes,
    expected_cells: list[str],
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, "empty.fa", data)
    io, stdout, stderr = make_io()

    assert stats.run(["empty.fa"], io) == 0

    assert stdout.getvalue().splitlines()[1].split() == expected_cells
    assert stdout.getvalue().endswith("\n")
    assert stdout.writes == [stdout.getvalue()]
    assert stderr.getvalue() == ""


@pytest.mark.parametrize(
    ("data", "expected_max", "expected_n50", "expected_gc"),
    [
        (b">a\nAAAA\n>b\nCC\n>c\nGG\n", "4", "4", "50.0"),
        (b">a\nAAA\n>b\nCCC\n>c\nGG\n", "3", "3", "62.5"),
    ],
)
def test_stats_n50_uses_the_first_length_reaching_half_of_total_bases(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    data: bytes,
    expected_max: str,
    expected_n50: str,
    expected_gc: str,
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, "lengths.fa", data)
    io, stdout, _ = make_io()

    assert stats.run(["lengths.fa"], io) == 0

    assert stdout.getvalue().splitlines()[1].split() == [
        "lengths.fa",
        "3",
        "8",
        "2",
        expected_max,
        "2.7",
        expected_n50,
        expected_gc,
        "0",
    ]


def test_stats_keeps_mixed_format_files_in_argument_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, "first.fa", b">a\nAC\n>b empty\n")
    write_fixture(tmp_path, "second.fq", b"@b desc\nGGN\n+ignored\n!!!\n")
    io, stdout, stderr = make_io()

    assert stats.run(["second.fq", "first.fa"], io) == 0

    assert stdout.getvalue() == (
        "file       records  bases  min  max  mean  n50    gc%  n\n"
        "second.fq        1      3    3    3   3.0    3  100.0  1\n"
        "first.fa         2      2    0    2   1.0    2   50.0  0\n"
    )
    assert stdout.writes == [stdout.getvalue()]
    assert stderr.getvalue() == ""


def test_stats_reads_standard_input_as_a_single_successful_write() -> None:
    io, stdout, stderr = make_io(b">s\nacgt\n")

    assert stats.run(["-"], io) == 0

    assert stdout.getvalue() == (
        "file     records  bases  min  max  mean  n50   gc%  n\n"
        "<stdin>        1      4    4    4   4.0    4  50.0  0\n"
    )
    assert stdout.writes == [stdout.getvalue()]
    assert stderr.getvalue() == ""


def test_stats_reports_missing_file_before_reading_any_input() -> None:
    io, _, _ = make_io()

    assert_failure(io, [], "missing FILE")


def test_stats_rejects_duplicate_stdin_before_attempting_file_reads(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    io, _, _ = make_io()

    assert_failure(io, ["missing.fa", "-", "-"], "standard input given twice")
    assert io.stdin.tell() == 0


@pytest.mark.parametrize(
    ("argument", "message"),
    [
        ("-x", 'unknown option "-x"'),
        ("--", 'unknown option "--"'),
        ("--help", 'unknown option "--help"'),
    ],
)
def test_stats_safely_quotes_unknown_options(argument: str, message: str) -> None:
    io, _, _ = make_io()

    assert_failure(io, [argument], message)


def test_conversion_width_1001_has_the_exact_deferred_option_error() -> None:
    io, stdout, stderr = make_io()

    with pytest.raises(SeqError) as caught:
        convert.run_revcomp(["--width", "1001", "missing.fa"], io)

    assert caught.value.message == (
        'invalid --width "1001" (expected an integer from 0 to 1000)'
    )
    assert stdout.getvalue() == ""
    assert stderr.getvalue() == ""
    assert stdout.writes == []
    assert stderr.writes == []


def test_translate_flag_assignment_has_the_exact_unknown_option_error() -> None:
    io, stdout, stderr = make_io()

    with pytest.raises(SeqError) as caught:
        convert.run_translate(["--to-stop=yes"], io)

    assert caught.value.message == 'unknown option "--to-stop=yes"'
    assert stdout.getvalue() == ""
    assert stderr.getvalue() == ""
    assert stdout.writes == []
    assert stderr.writes == []


@pytest.mark.parametrize(
    ("failure", "message"),
    [
        ("missing", "missing.fa: cannot read the file"),
        ("format", 'format.fa:2: invalid character "R" in sequence'),
        ("utf8", "utf8.fa: input is not valid UTF-8"),
    ],
)
def test_stats_reports_first_file_failure_and_keeps_output_atomic(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    failure: str,
    message: str,
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, "good.fa", b">good\nACGT\n")
    if failure == "format":
        write_fixture(tmp_path, "format.fa", b">bad\nACGR\n")
        args = ["good.fa", "format.fa", "missing-after-error.fa"]
    elif failure == "utf8":
        write_fixture(tmp_path, "utf8.fa", b">s\nAC\xff\n")
        args = ["good.fa", "utf8.fa", "missing-after-error.fa"]
    else:
        args = ["good.fa", "missing.fa"]
    io, _, _ = make_io()

    assert_failure(io, args, message)


def test_stats_passes_exact_malformed_stdin_bytes_through_existing_decoder() -> None:
    io, _, _ = make_io(b">s\nAC\xff\n")

    assert_failure(io, ["-"], "<stdin>: input is not valid UTF-8")


def test_stats_passes_exact_malformed_file_bytes_through_existing_decoder(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, "bad.fa", b">s\nAC\xff\n")
    io, _, _ = make_io()

    assert_failure(io, ["bad.fa"], "bad.fa: input is not valid UTF-8")


def test_stats_normalizes_bom_crlf_lowercase_and_missing_final_lf(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, "normalized.fa", b"\xef\xbb\xbf>s desc\r\nacgt\r\nnn")
    io, stdout, stderr = make_io()

    assert stats.run(["normalized.fa"], io) == 0

    assert stdout.getvalue().splitlines()[1].split() == [
        "normalized.fa",
        "1",
        "6",
        "6",
        "6",
        "6.0",
        "6",
        "50.0",
        "2",
    ]
    assert stdout.writes == [stdout.getvalue()]
    assert stderr.getvalue() == ""

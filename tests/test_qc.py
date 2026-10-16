import re
from inspect import signature
from io import BytesIO, StringIO
from pathlib import Path
from typing import BinaryIO

import pytest

from seqmill.core.seq import SeqError
from seqmill.formats.files import Io
from seqmill.tools.qc import run


class WriteSpy(StringIO):
    def __init__(self, name: str, events: list[tuple[str, str]]) -> None:
        super().__init__()
        self.name = name
        self.events = events

    def write(self, text: str) -> int:
        self.events.append((self.name, text))
        return super().write(text)


class ReadSpy(BytesIO):
    def __init__(self, initial_bytes: bytes = b"") -> None:
        super().__init__(initial_bytes)
        self.read_calls = 0

    def read(self, size: int = -1) -> bytes:
        self.read_calls += 1
        return super().read(size)


def make_io(data: bytes = b"", *, stdin: BinaryIO | None = None) -> tuple[Io, list[tuple[str, str]]]:
    events: list[tuple[str, str]] = []
    io = Io(
        stdin if stdin is not None else BytesIO(data),
        WriteSpy("stdout", events),
        WriteSpy("stderr", events),
    )
    return io, events


def assert_no_writes(io: Io) -> None:
    assert io.stdout.getvalue() == ""
    assert io.stderr.getvalue() == ""
    assert io.stdout.events == []
    assert io.stderr.events == []


def assert_error(argv: list[str], io: Io, message: str) -> None:
    with pytest.raises(SeqError) as caught:
        run(argv, io)
    assert caught.value.message == message
    assert_no_writes(io)


def test_run_has_the_exact_public_signature() -> None:
    run_signature = signature(run)
    assert list(run_signature.parameters) == ["argv", "io"]
    assert run_signature.parameters["argv"].annotation == list[str]
    assert run_signature.parameters["io"].annotation is Io
    assert run_signature.return_annotation is int


def test_default_options_keep_a_low_quality_read_at_minimum_length_one(
    tmp_path: Path,
) -> None:
    path = tmp_path / "default.fq"
    path.write_bytes(b"@r desc\nA\n+\n!\n")
    io, events = make_io()

    assert run([str(path)], io) == 0

    assert io.stdout.getvalue() == "@r desc\nA\n+\n!\n"
    assert io.stderr.getvalue() == "kept 1 of 1 reads: 0 too short, 0 low quality, 0 bases trimmed\n"
    assert events == [
        ("stdout", "@r desc\nA\n+\n!\n"),
        ("stderr", "kept 1 of 1 reads: 0 too short, 0 low quality, 0 bases trimmed\n"),
    ]


@pytest.mark.parametrize(
    ("option", "value", "message"),
    [
        ("--trim-end", "-1", 'invalid --trim-end "-1" (expected an integer from 0 to 93)'),
        ("--trim-end", "94", 'invalid --trim-end "94" (expected an integer from 0 to 93)'),
        ("--min-length", "-1", 'invalid --min-length "-1" (expected an integer from 0 to 1000000)'),
        ("--min-length", "1000001", 'invalid --min-length "1000001" (expected an integer from 0 to 1000000)'),
        ("--min-quality", "-1", 'invalid --min-quality "-1" (expected an integer from 0 to 93)'),
        ("--min-quality", "94", 'invalid --min-quality "94" (expected an integer from 0 to 93)'),
        ("--trim-end", "-0", 'invalid --trim-end "-0" (expected an integer from 0 to 93)'),
        ("--trim-end", "+3", 'invalid --trim-end "+3" (expected an integer from 0 to 93)'),
        ("--trim-end", "007", 'invalid --trim-end "007" (expected an integer from 0 to 93)'),
        ("--min-length", "1.5", 'invalid --min-length "1.5" (expected an integer from 0 to 1000000)'),
        ("--min-quality", "1000000000", 'invalid --min-quality "1000000000" (expected an integer from 0 to 93)'),
    ],
)
def test_option_values_reject_invalid_integer_forms_and_bounds(
    tmp_path: Path, option: str, value: str, message: str
) -> None:
    io, _ = make_io()
    assert_error([f"{option}={value}", str(tmp_path / "unused.fq")], io, message)


def test_inclusive_option_boundaries_are_accepted(tmp_path: Path) -> None:
    empty_read = tmp_path / "empty.fq"
    empty_read.write_bytes(b"@empty\n\n+\n\n")
    io, _ = make_io()

    assert run(
        ["--trim-end=0", "--min-length=0", "--min-quality=93", str(empty_read)],
        io,
    ) == 0

    assert io.stdout.getvalue() == "@empty\n\n+\n\n"
    assert io.stderr.getvalue() == "kept 1 of 1 reads: 0 too short, 0 low quality, 0 bases trimmed\n"


def test_upper_trim_and_length_option_boundaries_are_accepted(tmp_path: Path) -> None:
    path = tmp_path / "one.fq"
    path.write_bytes(b"@r\nA\n+\n~\n")
    io, _ = make_io()

    assert run(["--trim-end=93", "--min-length=1000000", str(path)], io) == 0

    assert io.stdout.getvalue() == ""
    assert io.stderr.getvalue() == "kept 0 of 1 reads: 1 too short, 0 low quality, 0 bases trimmed\n"


def test_separate_equals_and_interspersed_options_are_equivalent(tmp_path: Path) -> None:
    path = tmp_path / "mixed.fq"
    path.write_bytes(b"@r\nAC\n+\nI5\n")
    first_io, _ = make_io()
    second_io, _ = make_io()

    assert run(["--trim-end", "20", str(path), "--min-length=2", "--min-quality", "20"], first_io) == 0
    assert run(["--trim-end=20", "--min-length", "2", "--min-quality=20", str(path)], second_io) == 0

    assert first_io.stdout.getvalue() == second_io.stdout.getvalue() == "@r\nAC\n+\nI5\n"
    assert first_io.stderr.getvalue() == second_io.stderr.getvalue() == (
        "kept 1 of 1 reads: 0 too short, 0 low quality, 0 bases trimmed\n"
    )


def test_last_repeated_value_wins_even_when_earlier_value_is_invalid(tmp_path: Path) -> None:
    path = tmp_path / "repeat.fq"
    path.write_bytes(b"@r\nA\n+\nI\n")
    io, _ = make_io()

    assert run(["--trim-end=94", "--trim-end=0", str(path)], io) == 0

    assert io.stdout.getvalue() == "@r\nA\n+\nI\n"
    assert io.stderr.getvalue() == "kept 1 of 1 reads: 0 too short, 0 low quality, 0 bases trimmed\n"


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        (["--help"], 'unknown option "--help"'),
        (["--"], 'unknown option "--"'),
        (["-x"], 'unknown option "-x"'),
        (["--mystery=1"], 'unknown option "--mystery=1"'),
        (["--mystery\n=é"], 'unknown option "--mystery?=?"'),
        (["--min-quality=1\n2"], 'invalid --min-quality "1?2" (expected an integer from 0 to 93)'),
    ],
)
def test_syntax_errors_and_safe_quoting_leave_both_streams_untouched(
    tmp_path: Path, argv: list[str], message: str
) -> None:
    io, _ = make_io()
    assert_error([*argv, str(tmp_path / "unused.fq")], io, message)


def test_missing_option_value_is_reported_when_the_option_is_last() -> None:
    io, _ = make_io()
    assert_error(["--trim-end"], io, "--trim-end needs a value")


def test_syntax_error_precedes_deferred_invalid_value(tmp_path: Path) -> None:
    io, _ = make_io()
    assert_error(
        ["--trim-end=bad", "--unknown", str(tmp_path / "unused.fq")],
        io,
        'unknown option "--unknown"',
    )


def test_option_validation_precedes_missing_file() -> None:
    io, _ = make_io()
    assert_error(
        ["--min-quality=94"],
        io,
        'invalid --min-quality "94" (expected an integer from 0 to 93)',
    )


def test_invalid_options_are_reported_in_usage_order() -> None:
    io, _ = make_io()
    assert_error(
        ["--min-quality=94", "--min-length=-1", "--trim-end=94", "reads.fq"],
        io,
        'invalid --trim-end "94" (expected an integer from 0 to 93)',
    )


def test_min_length_validation_precedes_min_quality_validation() -> None:
    io, _ = make_io()
    assert_error(
        ["--min-quality=94", "--min-length=-1", "reads.fq"],
        io,
        'invalid --min-length "-1" (expected an integer from 0 to 1000000)',
    )


def test_missing_file_and_duplicate_stdin_are_checked_before_reads() -> None:
    missing_io, _ = make_io()
    assert_error([], missing_io, "missing FILE")

    stdin = ReadSpy(b"@r\nA\n+\nI\n")
    duplicate_io, _ = make_io(stdin=stdin)
    assert_error(["-", "-"], duplicate_io, "standard input given twice")
    assert stdin.read_calls == 0


@pytest.mark.parametrize("contents", [b"", b" \t\r\n\t\n"])
def test_empty_and_blank_only_fasta_are_accepted_as_no_reads(
    tmp_path: Path, contents: bytes
) -> None:
    path = tmp_path / "empty.fa"
    path.write_bytes(contents)
    io, events = make_io()

    assert run([str(path)], io) == 0

    assert io.stdout.getvalue() == ""
    assert io.stderr.getvalue() == "kept 0 of 0 reads: 0 too short, 0 low quality, 0 bases trimmed\n"
    assert events == [
        ("stdout", ""),
        ("stderr", "kept 0 of 0 reads: 0 too short, 0 low quality, 0 bases trimmed\n"),
    ]


@pytest.mark.parametrize(
    ("contents", "line", "message"),
    [
        (b">header only\n", None, "not a FASTQ file"),
        (b">populated\nACGT\n", None, "not a FASTQ file"),
        (b">bad\nACR\n", 2, 'invalid character "R" in sequence'),
    ],
)
def test_fasta_records_are_rejected_after_successful_parsing(
    tmp_path: Path, contents: bytes, line: int | None, message: str
) -> None:
    path = tmp_path / "input.fa"
    path.write_bytes(contents)
    io, _ = make_io()
    full_message = (
        f"{path.as_posix()}:{line}: {message}" if line is not None else f"{path.as_posix()}: {message}"
    )

    assert_error([str(path)], io, full_message)


def test_malformed_fastq_error_is_preserved_exactly(tmp_path: Path) -> None:
    path = tmp_path / "bad.fq"
    path.write_bytes(b"@r\nAC\nx\nII\n")
    io, _ = make_io()

    assert_error([str(path)], io, f"{path.as_posix()}:3: expected + line")


def test_unreadable_file_error_uses_display_path(tmp_path: Path) -> None:
    path = tmp_path / "does-not-exist.fq"
    io, _ = make_io()

    assert_error([str(path)], io, f"{path.as_posix()}: cannot read the file")


def test_invalid_utf8_is_decoded_from_binary_stdin_before_format_parsing() -> None:
    io, _ = make_io(data=b">s\nAC\xff\n")

    assert_error(["-"], io, "<stdin>: input is not valid UTF-8")


def test_invalid_utf8_file_is_reported_before_format_validation(tmp_path: Path) -> None:
    path = tmp_path / "bad-encoding.fq"
    path.write_bytes(b"@r\nAC\xff\n+\n!!\n")
    io, _ = make_io()

    assert_error([str(path)], io, f"{path.as_posix()}: input is not valid UTF-8")


def test_bom_crlf_lowercase_and_final_line_without_lf_are_normalized(tmp_path: Path) -> None:
    path = tmp_path / "normalized.fq"
    path.write_bytes(b"\xef\xbb\xbf@r lane 1\r\nacgt\r\n+ignored\r\nIIII")
    io, _ = make_io()

    assert run([str(path)], io) == 0

    assert io.stdout.getvalue() == "@r lane 1\nACGT\n+\nIIII\n"
    assert io.stderr.getvalue() == "kept 1 of 1 reads: 0 too short, 0 low quality, 0 bases trimmed\n"


def test_single_stdin_file_uses_the_binary_input_channel() -> None:
    io, _ = make_io(data=b"@r\r\nac\r\n+\r\nI5")

    assert run(["-"], io) == 0

    assert io.stdout.getvalue() == "@r\nAC\n+\nI5\n"
    assert io.stderr.getvalue() == "kept 1 of 1 reads: 0 too short, 0 low quality, 0 bases trimmed\n"


def test_later_file_failure_keeps_both_streams_atomic(tmp_path: Path) -> None:
    first = tmp_path / "first.fq"
    later = tmp_path / "later.fq"
    first.write_bytes(b"@first\nAC\n+\nII\n")
    later.write_bytes(b"@bad\nA\nx\n!\n")
    io, _ = make_io()

    assert_error(
        [str(first), str(later)],
        io,
        f"{later.as_posix()}:3: expected + line",
    )


def test_non_fastq_first_file_wins_before_a_later_malformed_file(tmp_path: Path) -> None:
    first = tmp_path / "first.fa"
    later = tmp_path / "later.fq"
    first.write_bytes(b">dna\nAC\n")
    later.write_bytes(b"@bad\nA\nx\n!\n")
    io, _ = make_io()

    assert_error(
        [str(first), str(later)],
        io,
        f"{first.as_posix()}: not a FASTQ file",
    )


@pytest.mark.parametrize(
    ("options", "contents", "expected_stdout", "expected_summary"),
    [
        (
            ["--trim-end=20", "--min-length=3"],
            b"@r\nACGT\n+\nI!I#\n",
            "@r\nACG\n+\nI!I\n",
            "kept 1 of 1 reads: 0 too short, 0 low quality, 1 bases trimmed\n",
        ),
        (
            ["--trim-end=20"],
            b"@r\nAC\n+\nI5\n",
            "@r\nAC\n+\nI5\n",
            "kept 1 of 1 reads: 0 too short, 0 low quality, 0 bases trimmed\n",
        ),
        (
            ["--trim-end=0"],
            b"@r\nAC\n+\n!!\n",
            "@r\nAC\n+\n!!\n",
            "kept 1 of 1 reads: 0 too short, 0 low quality, 0 bases trimmed\n",
        ),
        (
            ["--trim-end=20", "--min-length=0"],
            b"@r\nAC\n+\n##\n",
            "@r\n\n+\n\n",
            "kept 1 of 1 reads: 0 too short, 0 low quality, 2 bases trimmed\n",
        ),
        (
            ["--min-quality=20"],
            b"@r\nAC\n+\n55\n",
            "@r\nAC\n+\n55\n",
            "kept 1 of 1 reads: 0 too short, 0 low quality, 0 bases trimmed\n",
        ),
        (
            ["--min-quality=30"],
            b"@r\nAC\n+\n5I\n",
            "@r\nAC\n+\n5I\n",
            "kept 1 of 1 reads: 0 too short, 0 low quality, 0 bases trimmed\n",
        ),
        (
            ["--min-quality=27"],
            b"@r\nACG\n+\n!II\n",
            "",
            "kept 0 of 1 reads: 0 too short, 1 low quality, 0 bases trimmed\n",
        ),
        (
            ["--trim-end=20", "--min-length=1"],
            b"@r\nAC\n+\n##\n",
            "",
            "kept 0 of 1 reads: 1 too short, 0 low quality, 2 bases trimmed\n",
        ),
        (
            ["--trim-end=20", "--min-length=1", "--min-quality=93"],
            b"@r\nA\n+\n!\n",
            "",
            "kept 0 of 1 reads: 1 too short, 0 low quality, 1 bases trimmed\n",
        ),
    ],
)
def test_trimming_quality_and_rejection_boundaries(
    tmp_path: Path,
    options: list[str],
    contents: bytes,
    expected_stdout: str,
    expected_summary: str,
) -> None:
    path = tmp_path / "boundary.fq"
    path.write_bytes(contents)
    io, events = make_io()

    assert run([*options, str(path)], io) == 0

    assert io.stdout.getvalue() == expected_stdout
    assert io.stderr.getvalue() == expected_summary
    assert events == [("stdout", expected_stdout), ("stderr", expected_summary)]


def test_empty_read_passes_zero_length_and_quality_thresholds(tmp_path: Path) -> None:
    path = tmp_path / "empty-read.fq"
    path.write_bytes(b"@empty\n\n+\n\n")
    io, _ = make_io()

    assert run(["--min-length=0", "--min-quality=93", str(path)], io) == 0

    assert io.stdout.getvalue() == "@empty\n\n+\n\n"
    assert io.stderr.getvalue() == "kept 1 of 1 reads: 0 too short, 0 low quality, 0 bases trimmed\n"


def test_empty_read_is_too_short_with_the_default_minimum_length(tmp_path: Path) -> None:
    path = tmp_path / "empty-read.fq"
    path.write_bytes(b"@empty\n\n+\n\n")
    io, _ = make_io()

    assert run([str(path)], io) == 0

    assert io.stdout.getvalue() == ""
    assert io.stderr.getvalue() == "kept 0 of 1 reads: 1 too short, 0 low quality, 0 bases trimmed\n"


def test_supplied_five_read_example_matches_exact_output_and_summary(tmp_path: Path) -> None:
    path = tmp_path / "reads.fq"
    path.write_bytes(
        b"@r1 lane1\nACGTACGT\n+\nIIIIII##\n"
        b"@r2\nGGGGCCCC\n+\n5555!!!!\n"
        b"@r3\nACGT\n+\nII##\n"
        b"@r4\nTTTTAAAA\n+\n++++++++\n"
        b"@r5\nCCCCGGGG\n+\nIIIIIIII\n"
    )
    io, events = make_io()

    assert run(
        ["--trim-end", "20", "--min-length", "4", "--min-quality", "25", str(path)],
        io,
    ) == 0

    expected_stdout = "@r1 lane1\nACGTAC\n+\nIIIIII\n@r5\nCCCCGGGG\n+\nIIIIIIII\n"
    expected_stderr = "kept 2 of 5 reads: 2 too short, 1 low quality, 16 bases trimmed\n"
    assert io.stdout.getvalue() == expected_stdout
    assert io.stderr.getvalue() == expected_stderr
    assert events == [("stdout", expected_stdout), ("stderr", expected_stderr)]


def test_multifile_accounting_preserves_metadata_order_and_dropped_trims(
    tmp_path: Path,
) -> None:
    first = tmp_path / "first.fq"
    second = tmp_path / "second.fq"
    first.write_bytes(b"@dup left\nACGT\n+\nIIII\n@short lane\nA\n+\n!\n")
    second.write_bytes(
        b"@dup right\nGG\n+\n##\n@low\nACG\n+\n!II\n@keep\nTT\n+\nII\n"
    )
    io, events = make_io()

    assert run(
        ["--trim-end=20", "--min-length=2", "--min-quality=27", str(first), str(second)],
        io,
    ) == 0

    expected_stdout = "@dup left\nACGT\n+\nIIII\n@keep\nTT\n+\nII\n"
    expected_stderr = "kept 2 of 5 reads: 2 too short, 1 low quality, 3 bases trimmed\n"
    assert io.stdout.getvalue() == expected_stdout
    assert io.stderr.getvalue() == expected_stderr
    counts = re.fullmatch(
        r"kept (?P<kept>\d+) of (?P<total>\d+) reads: "
        r"(?P<short>\d+) too short, (?P<low>\d+) low quality, "
        r"\d+ bases trimmed\n",
        io.stderr.getvalue(),
    )
    assert counts is not None
    assert (
        int(counts["kept"]) + int(counts["short"]) + int(counts["low"])
        == int(counts["total"])
    )
    assert events == [("stdout", expected_stdout), ("stderr", expected_stderr)]


def test_all_rejected_reads_still_write_empty_stdout_then_one_summary(tmp_path: Path) -> None:
    path = tmp_path / "rejected.fq"
    path.write_bytes(b"@short\nA\n+\n!\n@low\nAC\n+\n!I\n")
    io, events = make_io()

    assert run(["--trim-end=20", "--min-length=1", "--min-quality=30", str(path)], io) == 0

    expected_summary = "kept 0 of 2 reads: 1 too short, 1 low quality, 1 bases trimmed\n"
    assert io.stdout.getvalue() == ""
    assert io.stderr.getvalue() == expected_summary
    assert events == [("stdout", ""), ("stderr", expected_summary)]

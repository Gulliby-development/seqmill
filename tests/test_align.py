from io import BytesIO, StringIO
from pathlib import Path

import pytest

from seqmill.core.seq import SeqError
from seqmill.formats.files import Io
from seqmill.tools.align import run


class CountingText:
    def __init__(self) -> None:
        self.writes: list[str] = []

    def write(self, text: str) -> int:
        self.writes.append(text)
        return len(text)

    def getvalue(self) -> str:
        return "".join(self.writes)


def make_io(stdin: bytes = b"") -> tuple[Io, CountingText, StringIO]:
    stdout = CountingText()
    stderr = StringIO()
    return Io(BytesIO(stdin), stdout, stderr), stdout, stderr


def write_fasta(
    tmp_path: Path, filename: str, name: str, sequence: str, description: str = ""
) -> Path:
    path = tmp_path / filename
    suffix = f" {description}" if description else ""
    path.write_text(f">{name}{suffix}\n{sequence}\n", encoding="utf-8")
    return path


def run_pair(
    tmp_path: Path,
    a_sequence: str,
    b_sequence: str,
    options: tuple[str, ...] = (),
    a_description: str = "",
    b_description: str = "",
) -> tuple[int, CountingText, StringIO]:
    a_path = write_fasta(tmp_path, "a.fa", "alpha", a_sequence, a_description)
    b_path = write_fasta(tmp_path, "b.fa", "beta", b_sequence, b_description)
    io, stdout, stderr = make_io()
    status = run([*options, a_path.as_posix(), b_path.as_posix()], io)
    assert len(stdout.writes) == 1
    assert stdout.getvalue().endswith("\n")
    assert "\r" not in stdout.getvalue()
    assert stderr.getvalue() == ""
    return status, stdout, stderr


def assert_error(
    argv: list[str],
    io: Io,
    stdout: CountingText,
    stderr: StringIO,
    message: str,
) -> None:
    with pytest.raises(SeqError) as caught:
        run(argv, io)
    assert caught.value.message == message
    assert stdout.writes == []
    assert stderr.getvalue() == ""


def test_global_example_matches_the_supplied_output_exactly(tmp_path: Path) -> None:
    status, stdout, stderr = run_pair(
        tmp_path,
        "GATTACAGATTACA",
        "GCATGCAGATTCA",
        ("--width", "10"),
        a_description="first sequence",
        b_description="second sequence",
    )

    assert status == 0
    assert stdout.getvalue() == (
        "score: 5\n"
        "length: 14\n"
        "identity: 10/14 (71.4%)\n"
        "A: alpha 1-14\n"
        "B: beta 1-13\n"
        "\n"
        "GATTACAGAT\n"
        "|xx|x|||||\n"
        "GCATGCAGAT\n"
        "\n"
        "TACA\n"
        "|-||\n"
        "T-CA\n"
    )
    assert len(stdout.writes) == 1
    assert stderr.getvalue() == ""


def test_local_example_matches_the_supplied_output_exactly(tmp_path: Path) -> None:
    status, stdout, stderr = run_pair(
        tmp_path,
        "GATTACAGATTACA",
        "GCATGCAGATTCA",
        ("--mode", "local"),
        a_description="first sequence",
        b_description="second sequence",
    )

    assert status == 0
    assert stdout.getvalue() == (
        "score: 6\n"
        "length: 6\n"
        "identity: 6/6 (100.0%)\n"
        "A: alpha 6-11\n"
        "B: beta 6-11\n"
        "\n"
        "CAGATT\n"
        "||||||\n"
        "CAGATT\n"
    )
    assert len(stdout.writes) == 1
    assert stderr.getvalue() == ""


def test_defaults_align_all_bases_and_ignore_descriptions(tmp_path: Path) -> None:
    status, stdout, _ = run_pair(
        tmp_path, "ACG", "ACG", a_description="alpha details", b_description="beta details"
    )

    assert status == 0
    assert stdout.getvalue() == (
        "score: 3\n"
        "length: 3\n"
        "identity: 3/3 (100.0%)\n"
        "A: alpha 1-3\n"
        "B: beta 1-3\n"
        "\n"
        "ACG\n"
        "|||\n"
        "ACG\n"
    )


@pytest.mark.parametrize(
    ("option", "expected"),
    [
        (
            "--width=1",
            (
                "score: 3\nlength: 3\nidentity: 3/3 (100.0%)\n"
                "A: alpha 1-3\nB: beta 1-3\n\nA\n|\nA\n\nC\n|\nC\n\nG\n|\nG\n"
            ),
        ),
        (
            "--width=2",
            (
                "score: 3\nlength: 3\nidentity: 3/3 (100.0%)\n"
                "A: alpha 1-3\nB: beta 1-3\n\nAC\n||\nAC\n\nG\n|\nG\n"
            ),
        ),
        (
            "--width=10",
            (
                "score: 3\nlength: 3\nidentity: 3/3 (100.0%)\n"
                "A: alpha 1-3\nB: beta 1-3\n\nACG\n|||\nACG\n"
            ),
        ),
        (
            "--width=60",
            (
                "score: 3\nlength: 3\nidentity: 3/3 (100.0%)\n"
                "A: alpha 1-3\nB: beta 1-3\n\nACG\n|||\nACG\n"
            ),
        ),
        (
            "--width=1000",
            (
                "score: 3\nlength: 3\nidentity: 3/3 (100.0%)\n"
                "A: alpha 1-3\nB: beta 1-3\n\nACG\n|||\nACG\n"
            ),
        ),
    ],
)
def test_width_boundaries_make_exact_lf_blocks(
    tmp_path: Path, option: str, expected: str
) -> None:
    status, stdout, _ = run_pair(tmp_path, "ACG", "ACG", (option,))

    assert status == 0
    assert stdout.getvalue() == expected
    assert stdout.getvalue().endswith("\n")
    assert "\r" not in stdout.getvalue()


def test_option_forms_repetition_and_interspersed_operands(tmp_path: Path) -> None:
    a_path = write_fasta(tmp_path, "a.fa", "alpha", "AC")
    b_path = write_fasta(tmp_path, "b.fa", "beta", "AC")
    io, stdout, stderr = make_io()

    status = run(
        ["--match", "2", a_path.as_posix(), "--match=1", "--width", "1000", b_path.as_posix()],
        io,
    )

    assert status == 0
    assert stdout.getvalue().startswith("score: 2\nlength: 2\n")
    assert len(stdout.writes) == 1
    assert stderr.getvalue() == ""


def test_valid_scoring_endpoints_are_accepted(tmp_path: Path) -> None:
    status, stdout, _ = run_pair(
        tmp_path,
        "A",
        "A",
        ("--match=100", "--mismatch=-100", "--gap=-100", "--width=1000"),
    )

    assert status == 0
    assert stdout.getvalue().startswith(
        "score: 100\nlength: 1\nidentity: 1/1 (100.0%)\n"
    )


@pytest.mark.parametrize(
    "options",
    [
        ("--mode=global",),
        ("--mode", "local"),
        ("--match=1",),
        ("--match", "1"),
        ("--match=100",),
        ("--match", "100"),
        ("--mismatch=0",),
        ("--mismatch", "0"),
        ("--mismatch=-100",),
        ("--mismatch", "-100"),
        ("--gap=-1",),
        ("--gap", "-1"),
        ("--gap=-100",),
        ("--gap", "-100"),
        ("--width=1",),
        ("--width", "1"),
        ("--width=1000",),
        ("--width", "1000"),
    ],
)
def test_option_endpoints_are_accepted_with_both_value_syntaxes(
    tmp_path: Path, options: tuple[str, ...]
) -> None:
    status, _, _ = run_pair(tmp_path, "A", "A", options)

    assert status == 0


@pytest.mark.parametrize(
    ("option", "message"),
    [
        ("--mode=LOCAL", 'invalid --mode "LOCAL" (expected global or local)'),
        ("--match=-0", 'invalid --match "-0" (expected an integer from 1 to 100)'),
        ("--match=+3", 'invalid --match "+3" (expected an integer from 1 to 100)'),
        ("--match=007", 'invalid --match "007" (expected an integer from 1 to 100)'),
        (
            "--match=1234567890",
            'invalid --match "1234567890" (expected an integer from 1 to 100)',
        ),
        ("--match=0", 'invalid --match "0" (expected an integer from 1 to 100)'),
        ("--mismatch=1", 'invalid --mismatch "1" (expected an integer from -100 to 0)'),
        ("--gap=0", 'invalid --gap "0" (expected an integer from -100 to -1)'),
        ("--width=0", 'invalid --width "0" (expected an integer from 1 to 1000)'),
        ("--width=1001", 'invalid --width "1001" (expected an integer from 1 to 1000)'),
    ],
)
def test_option_values_use_exact_bounds_and_integer_grammar(
    option: str, message: str
) -> None:
    io, stdout, stderr = make_io()

    assert_error([option], io, stdout, stderr, message)


@pytest.mark.parametrize(
    ("option", "message"),
    [
        ("--help", 'unknown option "--help"'),
        ("--", 'unknown option "--"'),
        ("-x", 'unknown option "-x"'),
        ("--mode=other", 'invalid --mode "other" (expected global or local)'),
    ],
)
def test_unknown_and_invalid_choice_options_are_exact(
    option: str, message: str
) -> None:
    io, stdout, stderr = make_io()

    assert_error([option], io, stdout, stderr, message)


def test_missing_option_value_precedes_missing_a() -> None:
    io, stdout, stderr = make_io()

    assert_error(["--gap"], io, stdout, stderr, "--gap needs a value")


def test_syntax_error_after_bad_value_precedes_deferred_validation() -> None:
    io, stdout, stderr = make_io()

    assert_error(
        ["--match=bad", "--unknown"],
        io,
        stdout,
        stderr,
        'unknown option "--unknown"',
    )


def test_option_value_errors_follow_spec_order() -> None:
    io, stdout, stderr = make_io()

    assert_error(
        [
            "--width=0",
            "--gap=0",
            "--mismatch=1",
            "--match=0",
            "--mode=other",
        ],
        io,
        stdout,
        stderr,
        'invalid --mode "other" (expected global or local)',
    )


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        ([], "missing A"),
        (["a.fa"], "missing B"),
        (["a.fa", "b.fa", "extra.fa"], "too many arguments"),
    ],
)
def test_operand_count_errors_precede_file_reads(argv: list[str], message: str) -> None:
    io, stdout, stderr = make_io()

    assert_error(argv, io, stdout, stderr, message)


def test_duplicate_stdin_is_rejected_before_reading() -> None:
    io, stdout, stderr = make_io()

    assert_error(["-", "-"], io, stdout, stderr, "standard input given twice")


@pytest.mark.parametrize(
    "bad_operand",
    [
        "A",
        "B",
    ],
)
def test_no_records_is_reported_for_either_operand(
    tmp_path: Path, bad_operand: str
) -> None:
    a_path = tmp_path / "a.fa"
    b_path = tmp_path / "b.fa"
    a_path.write_text(">alpha\nA\n", encoding="utf-8")
    b_path.write_text(">beta\nA\n", encoding="utf-8")
    (a_path if bad_operand == "A" else b_path).write_bytes(b"")
    io, stdout, stderr = make_io()

    assert_error(
        [a_path.as_posix(), b_path.as_posix()],
        io,
        stdout,
        stderr,
        f"{(a_path if bad_operand == 'A' else b_path).as_posix()}: no sequences",
    )


@pytest.mark.parametrize(
    ("bad_operand", "name"),
    [
        ("A", "alpha"),
        ("B", "beta"),
    ],
)
def test_empty_first_sequence_is_rejected_in_operand_order(
    tmp_path: Path, bad_operand: str, name: str
) -> None:
    a_path = tmp_path / "a.fa"
    b_path = tmp_path / "b.fa"
    a_path.write_text(">alpha\n\n", encoding="utf-8")
    b_path.write_text(">beta\nA\n", encoding="utf-8")
    if bad_operand == "B":
        a_path.write_text(">alpha\nA\n", encoding="utf-8")
        b_path.write_text(">beta\n\n", encoding="utf-8")
    io, stdout, stderr = make_io()
    path = a_path if bad_operand == "A" else b_path

    assert_error(
        [a_path.as_posix(), b_path.as_posix()],
        io,
        stdout,
        stderr,
        f'{path.as_posix()}: sequence "{name}" is empty',
    )


SEQUENCE_2000 = (
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
    "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA"
)
SEQUENCE_2001 = SEQUENCE_2000 + "A"


def test_2000_by_2000_alignment_is_supported(tmp_path: Path) -> None:
    assert len(SEQUENCE_2000) == 2000
    status, stdout, _ = run_pair(
        tmp_path,
        SEQUENCE_2000,
        SEQUENCE_2000,
        ("--width=1000",),
    )

    lines = stdout.getvalue().splitlines()
    assert status == 0
    assert lines[:6] == [
        "score: 2000",
        "length: 2000",
        "identity: 2000/2000 (100.0%)",
        "A: alpha 1-2000",
        "B: beta 1-2000",
        "",
    ]
    assert len(lines) == 13
    assert lines[6] == SEQUENCE_2000[:1000]
    assert len(lines[7]) == 1000
    assert set(lines[7]) == {"|"}
    assert lines[8] == SEQUENCE_2000[:1000]
    assert lines[9] == ""
    assert lines[10] == SEQUENCE_2000[1000:]
    assert len(lines[11]) == 1000
    assert set(lines[11]) == {"|"}
    assert lines[12] == SEQUENCE_2000[1000:]


@pytest.mark.parametrize(
    ("bad_operand", "expected_name"),
    [
        ("A", "alpha"),
        ("B", "beta"),
    ],
)
def test_2001_base_first_sequence_is_rejected_in_either_operand(
    tmp_path: Path, bad_operand: str, expected_name: str
) -> None:
    assert len(SEQUENCE_2001) == 2001
    a_path = write_fasta(
        tmp_path,
        "a.fa",
        "alpha",
        SEQUENCE_2001 if bad_operand == "A" else "A",
    )
    b_path = write_fasta(
        tmp_path,
        "b.fa",
        "beta",
        SEQUENCE_2001 if bad_operand == "B" else "A",
    )
    io, stdout, stderr = make_io()
    path = a_path if bad_operand == "A" else b_path

    assert_error(
        [a_path.as_posix(), b_path.as_posix()],
        io,
        stdout,
        stderr,
        f'{path.as_posix()}: sequence "{expected_name}" is longer than 2000 bases',
    )


def test_first_record_only_is_aligned_but_all_records_are_parsed(tmp_path: Path) -> None:
    a_path = tmp_path / "a.fa"
    b_path = write_fasta(tmp_path, "b.fa", "beta", "AC")
    a_path.write_text(">first details\nAC\n>second\nG\n", encoding="utf-8")
    io, stdout, stderr = make_io()

    status = run([a_path.as_posix(), b_path.as_posix()], io)

    assert status == 0
    assert stdout.getvalue().startswith(
        "score: 2\nlength: 2\nidentity: 2/2 (100.0%)\nA: first 1-2\nB: beta 1-2\n"
    )
    assert "second" not in stdout.getvalue()
    assert len(stdout.writes) == 1
    assert stderr.getvalue() == ""


def test_malformed_later_record_still_fails_before_output(tmp_path: Path) -> None:
    a_path = tmp_path / "a.fa"
    b_path = write_fasta(tmp_path, "b.fa", "beta", "AC")
    a_path.write_text(">alpha\nAC\n>later\nR\n", encoding="utf-8")
    io, stdout, stderr = make_io()

    assert_error(
        [a_path.as_posix(), b_path.as_posix()],
        io,
        stdout,
        stderr,
        f'{a_path.as_posix()}:4: invalid character "R" in sequence',
    )


def test_fasta_and_fastq_first_records_align_identically(tmp_path: Path) -> None:
    fasta_a = tmp_path / "a.fa"
    fasta_b = tmp_path / "b.fa"
    fastq_a = tmp_path / "a.fq"
    fastq_b = tmp_path / "b.fq"
    fasta_a.write_text(">alpha details\nAC\n", encoding="utf-8")
    fasta_b.write_text(">beta details\nAT\n", encoding="utf-8")
    fastq_a.write_text("@alpha details\nAC\n+\n!!\n", encoding="utf-8")
    fastq_b.write_text("@beta details\nAT\n+\n##\n", encoding="utf-8")
    fasta_io, fasta_stdout, fasta_stderr = make_io()
    fastq_io, fastq_stdout, fastq_stderr = make_io()

    fasta_status = run([fasta_a.as_posix(), fasta_b.as_posix()], fasta_io)
    fastq_status = run([fastq_a.as_posix(), fastq_b.as_posix()], fastq_io)

    assert fasta_status == 0
    assert fastq_status == 0
    assert fasta_stdout.getvalue() == fastq_stdout.getvalue()
    assert len(fasta_stdout.writes) == 1
    assert len(fastq_stdout.writes) == 1
    assert fasta_stderr.getvalue() == ""
    assert fastq_stderr.getvalue() == ""


def test_a_validation_error_precedes_an_unreadable_b(tmp_path: Path) -> None:
    a_path = tmp_path / "a.fa"
    a_path.write_text(">alpha\nACR\n", encoding="utf-8")
    b_path = tmp_path / "missing.fa"
    io, stdout, stderr = make_io()

    assert_error(
        [a_path.as_posix(), b_path.as_posix()],
        io,
        stdout,
        stderr,
        f'{a_path.as_posix()}:2: invalid character "R" in sequence',
    )


def test_unreadable_b_is_reported_after_validating_a(tmp_path: Path) -> None:
    a_path = write_fasta(tmp_path, "a.fa", "alpha", "AC")
    b_path = tmp_path / "missing.fa"
    io, stdout, stderr = make_io()

    assert_error(
        [a_path.as_posix(), b_path.as_posix()],
        io,
        stdout,
        stderr,
        f"{b_path.as_posix()}: cannot read the file",
    )


def test_invalid_utf8_file_is_reported_without_output(tmp_path: Path) -> None:
    a_path = write_fasta(tmp_path, "a.fa", "alpha", "AC")
    b_path = tmp_path / "bad.fa"
    b_path.write_bytes(b">beta\nAC\xff\n")
    io, stdout, stderr = make_io()

    assert_error(
        [a_path.as_posix(), b_path.as_posix()],
        io,
        stdout,
        stderr,
        f"{b_path.as_posix()}: input is not valid UTF-8",
    )


def test_invalid_utf8_stdin_uses_the_real_binary_reader(tmp_path: Path) -> None:
    b_path = write_fasta(tmp_path, "b.fa", "beta", "AC")
    io, stdout, stderr = make_io(b">alpha\nAC\xff\n")

    assert_error(
        ["-", b_path.as_posix()],
        io,
        stdout,
        stderr,
        "<stdin>: input is not valid UTF-8",
    )


def test_fastq_format_error_keeps_upstream_line_number(tmp_path: Path) -> None:
    a_path = tmp_path / "a.fq"
    b_path = write_fasta(tmp_path, "b.fa", "beta", "AC")
    a_path.write_text("@alpha\nAC\nnot-plus\n!!\n", encoding="utf-8")
    io, stdout, stderr = make_io()

    assert_error(
        [a_path.as_posix(), b_path.as_posix()],
        io,
        stdout,
        stderr,
        f"{a_path.as_posix()}:3: expected + line",
    )


def test_single_stdin_operand_is_loaded_through_io(tmp_path: Path) -> None:
    b_path = write_fasta(tmp_path, "b.fa", "beta", "AC")
    io, stdout, stderr = make_io(b">alpha\nAC\n")

    status = run(["-", b_path.as_posix()], io)

    assert status == 0
    assert stdout.getvalue().startswith(
        "score: 2\nlength: 2\nidentity: 2/2 (100.0%)\nA: alpha 1-2\nB: beta 1-2\n"
    )
    assert len(stdout.writes) == 1
    assert stderr.getvalue() == ""


@pytest.mark.parametrize(
    ("options", "expected"),
    [
        (("--mode", "local"), "score: 0\nno alignment\n"),
        (("--mode=local", "--mismatch=0"), "score: 0\nno alignment\n"),
    ],
)
def test_zero_score_local_alignment_writes_once_and_returns_one(
    tmp_path: Path, options: tuple[str, ...], expected: str
) -> None:
    status, stdout, stderr = run_pair(tmp_path, "A", "C", options)

    assert status == 1
    assert stdout.getvalue() == expected
    assert len(stdout.writes) == 1
    assert stderr.getvalue() == ""


def test_local_traceback_stops_at_a_zero_score_diagonal(tmp_path: Path) -> None:
    status, stdout, _ = run_pair(
        tmp_path, "CA", "TA", ("--mode=local", "--mismatch=0")
    )

    assert status == 0
    assert stdout.getvalue() == (
        "score: 1\n"
        "length: 1\n"
        "identity: 1/1 (100.0%)\n"
        "A: alpha 2-2\n"
        "B: beta 2-2\n"
        "\n"
        "A\n"
        "|\n"
        "A\n"
    )


@pytest.mark.parametrize(
    ("a_sequence", "b_sequence", "expected_start"),
    [
        ("AAA", "A", "A: alpha 1-1\nB: beta 1-1\n"),
        ("A", "AA", "A: alpha 1-1\nB: beta 1-1\n"),
    ],
)
def test_local_equal_maxima_choose_smallest_i_then_smallest_j(
    tmp_path: Path, a_sequence: str, b_sequence: str, expected_start: str
) -> None:
    status, stdout, _ = run_pair(
        tmp_path, a_sequence, b_sequence, ("--mode=local",)
    )

    assert status == 0
    assert stdout.getvalue().startswith(
        "score: 1\nlength: 1\nidentity: 1/1 (100.0%)\n"
    )
    assert expected_start in stdout.getvalue()
    assert stdout.getvalue().endswith("\n")


def test_local_alignment_uses_substring_coordinates(tmp_path: Path) -> None:
    status, stdout, _ = run_pair(
        tmp_path, "TTACGTT", "GGACGCC", ("--mode", "local")
    )

    assert status == 0
    assert stdout.getvalue() == (
        "score: 3\n"
        "length: 3\n"
        "identity: 3/3 (100.0%)\n"
        "A: alpha 3-5\n"
        "B: beta 3-5\n"
        "\n"
        "ACG\n"
        "|||\n"
        "ACG\n"
    )


@pytest.mark.parametrize(
    ("a_sequence", "b_sequence", "expected_a", "expected_markers", "expected_b"),
    [
        ("A", "CA", "-A", "-|", "CA"),
        ("CA", "A", "CA", "-|", "-A"),
        ("A", "AC", "A-", "|-", "AC"),
        ("AC", "A", "AC", "|-", "A-"),
        ("AC", "ATC", "A-C", "|-|", "ATC"),
        ("ATC", "AC", "ATC", "|-|", "A-C"),
    ],
)
def test_global_traceback_renders_leading_internal_and_trailing_gaps(
    tmp_path: Path,
    a_sequence: str,
    b_sequence: str,
    expected_a: str,
    expected_markers: str,
    expected_b: str,
) -> None:
    status, stdout, _ = run_pair(
        tmp_path,
        a_sequence,
        b_sequence,
        ("--match=100", "--mismatch=-100", "--gap=-1"),
    )

    assert status == 0
    assert stdout.getvalue().splitlines()[-3:] == [
        expected_a,
        expected_markers,
        expected_b,
    ]


def test_diagonal_wins_a_tie_with_up(tmp_path: Path) -> None:
    status, stdout, _ = run_pair(tmp_path, "AA", "A")

    assert status == 0
    assert stdout.getvalue().splitlines()[-3:] == ["AA", "-|", "-A"]


def test_up_wins_a_tie_with_left_when_diagonal_is_worse(tmp_path: Path) -> None:
    status, stdout, _ = run_pair(
        tmp_path,
        "AC",
        "CA",
        ("--match=1", "--mismatch=-100", "--gap=-1"),
    )

    assert status == 0
    assert stdout.getvalue().splitlines()[-3:] == ["-AC", "-|-", "CA-"]


@pytest.mark.parametrize(
    ("a_sequence", "b_sequence", "score", "identity", "marker"),
    [
        ("N", "N", 1, "1/1 (100.0%)", "|"),
        ("N", "A", -1, "0/1 (0.0%)", "x"),
    ],
)
def test_n_matches_only_n_and_markers_show_literal_equality(
    tmp_path: Path,
    a_sequence: str,
    b_sequence: str,
    score: int,
    identity: str,
    marker: str,
) -> None:
    status, stdout, _ = run_pair(tmp_path, a_sequence, b_sequence)

    assert status == 0
    assert stdout.getvalue().splitlines()[:3] == [
        f"score: {score}",
        "length: 1",
        f"identity: {identity}",
    ]
    assert stdout.getvalue().splitlines()[-2] == marker


def test_identity_uses_half_up_rounding_for_one_of_sixteen(tmp_path: Path) -> None:
    b_sequence = "ACCCCCCCCCCCCCCC"
    assert len(b_sequence) == 16
    status, stdout, _ = run_pair(
        tmp_path,
        "A",
        b_sequence,
        ("--match=100", "--mismatch=-100", "--gap=-1"),
    )

    assert status == 0
    assert stdout.getvalue() == (
        "score: 85\n"
        "length: 16\n"
        "identity: 1/16 (6.3%)\n"
        "A: alpha 1-1\n"
        "B: beta 1-16\n"
        "\n"
        "A---------------\n"
        "|---------------\n"
        "ACCCCCCCCCCCCCCC\n"
    )


def test_nonequal_global_alignment_reports_zero_identity(tmp_path: Path) -> None:
    status, stdout, _ = run_pair(tmp_path, "A", "C")

    assert status == 0
    assert "identity: 0/1 (0.0%)\n" in stdout.getvalue()

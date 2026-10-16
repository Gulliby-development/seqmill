from collections.abc import Callable
from inspect import signature
from io import BytesIO, StringIO
from pathlib import Path

import pytest

from seqmill.core.seq import SeqError
from seqmill.formats.files import Io
from seqmill.tools import convert


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


def assert_failure(
    runner: Callable[[list[str], Io], int], argv: list[str], io: Io, message: str
) -> None:
    with pytest.raises(SeqError) as caught:
        runner(argv, io)

    assert caught.value.message == message
    assert io.stdout.getvalue() == ""
    assert io.stderr.getvalue() == ""
    assert io.stdout.writes == []
    assert io.stderr.writes == []


def test_conversion_apis_have_the_exact_argv_io_signature() -> None:
    for runner in (convert.run_revcomp, convert.run_translate):
        run_signature = signature(runner)
        assert list(run_signature.parameters) == ["argv", "io"]
        assert run_signature.parameters["argv"].annotation == list[str]
        assert run_signature.parameters["io"].annotation is Io
        assert run_signature.return_annotation is int


def test_revcomp_reverses_the_complement_and_preserves_metadata() -> None:
    io, stdout, stderr = make_io(b">read1 lane one\nACGTN\n>read2\nA\n")
    assert convert.run_revcomp(["-"], io) == 0

    assert stdout.getvalue() == ">read1 lane one\nNACGT\n>read2\nT\n"
    assert stdout.writes == [stdout.getvalue()]
    assert stderr.getvalue() == ""
    assert stderr.writes == []


def test_revcomp_reverses_fastq_sequence_and_quality_without_wrapping() -> None:
    io, stdout, stderr = make_io(b"@r description\nACGTN\n+ignored\n!#5AI\n")

    assert convert.run_revcomp(["--width=1", "-"], io) == 0

    assert stdout.getvalue() == "@r description\nNACGT\n+\nIA5#!\n"
    assert stdout.writes == [stdout.getvalue()]
    assert stderr.getvalue() == ""


def test_revcomp_keeps_an_empty_fasta_record_header_only() -> None:
    io, stdout, stderr = make_io(b">empty description\n")

    assert convert.run_revcomp(["-"], io) == 0

    assert stdout.getvalue() == ">empty description\n"
    assert stdout.writes == [stdout.getvalue()]
    assert stderr.getvalue() == ""


def test_revcomp_keeps_an_empty_fastq_record_with_four_lines() -> None:
    io, stdout, stderr = make_io(b"@empty description\n\n+\n\n")

    assert convert.run_revcomp(["-"], io) == 0

    assert stdout.getvalue() == "@empty description\n\n+\n\n"
    assert stdout.writes == [stdout.getvalue()]
    assert stderr.getvalue() == ""


def test_revcomp_preserves_mixed_file_order_and_each_file_format(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, "first.fa", b">first desc\nAC\n")
    write_fixture(tmp_path, "middle.fq", b"@middle\nGT\n+\n!#\n")
    write_fixture(tmp_path, "last.fa", b">last\nN\n")
    io, stdout, stderr = make_io()

    assert convert.run_revcomp(["first.fa", "middle.fq", "last.fa"], io) == 0

    assert stdout.getvalue() == ">first desc\nGT\n@middle\nAC\n+\n#!\n>last\nN\n"
    assert stdout.writes == [stdout.getvalue()]
    assert stderr.getvalue() == ""


@pytest.mark.parametrize(
    ("option", "expected"),
    [
        ([], ">r\nNACGT\n"),
        (["--width=0"], ">r\nNACGT\n"),
        (["--width", "1"], ">r\nN\nA\nC\nG\nT\n"),
        (["--width=2"], ">r\nNA\nCG\nT\n"),
        (["--width", "3"], ">r\nNAC\nGT\n"),
        (["--width=60"], ">r\nNACGT\n"),
        (["--width=1000"], ">r\nNACGT\n"),
    ],
)
def test_revcomp_supports_width_forms_and_boundaries(
    option: list[str], expected: str
) -> None:
    io, stdout, stderr = make_io(b">r\nACGTN\n")

    assert convert.run_revcomp([*option, "-"], io) == 0

    assert stdout.getvalue() == expected
    assert stdout.writes == [expected]
    assert stderr.getvalue() == ""


@pytest.mark.parametrize(
    ("sequence", "width", "expected"),
    [
        ("A", "1", ">r\nT\n"),
        ("AC", "2", ">r\nGT\n"),
        ("ACG", "3", ">r\nCGT\n"),
        ("ACGT", "3", ">r\nACG\nT\n"),
        (
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            "60",
            ">r\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n",
        ),
        (
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            "60",
            ">r\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\nT\n",
        ),
        (
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            "60",
            ">r\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\nTTT\n",
        ),
    ],
)
def test_revcomp_wraps_short_and_boundary_length_sequences(
    sequence: str, width: str, expected: str
) -> None:
    io, stdout, _ = make_io(f">r\n{sequence}\n".encode("ascii"))

    assert convert.run_revcomp([f"--width={width}", "-"], io) == 0

    assert stdout.getvalue() == expected
    assert stdout.writes == [expected]


@pytest.mark.parametrize(
    ("frame", "expected"),
    [
        ("1", "MKP"),
        ("2", "*N"),
        ("3", "ET"),
        ("-1", "GFH"),
        ("-2", "GF"),
        ("-3", "VS"),
    ],
)
def test_translate_handles_all_six_frames(frame: str, expected: str) -> None:
    io, stdout, stderr = make_io(b">gene description\nATGAAACCC\n")

    assert convert.run_translate(["--frame", frame, "-"], io) == 0

    assert stdout.getvalue() == f">gene description\n{expected}\n"
    assert stdout.writes == [stdout.getvalue()]
    assert stderr.getvalue() == ""


def test_translate_renders_the_supplied_gene_example_exactly() -> None:
    io, stdout, stderr = make_io(
        b">gene1 test construct\n"
        b"GAATTCATGGCCAAAGGATCCTTTGGGAAATAGCCCGGGATGAAACCCTGAGTTTTCATTTGGGAAACAT\n"
    )

    assert convert.run_translate(["-"], io) == 0

    assert stdout.getvalue() == ">gene1 test construct\nEFMAKGSFGK*PGMKP*VFIWET\n"
    assert stdout.writes == [stdout.getvalue()]
    assert stderr.getvalue() == ""


@pytest.mark.parametrize(
    ("sequence", "option", "expected"),
    [
        ("ATGTAAGGG", [], "M*G"),
        ("ATGTAAGGG", ["--to-stop"], "M"),
        ("TAAATG", ["--to-stop"], ""),
        ("ATGNAA", [], "MX"),
        ("ATGA", [], "M"),
        ("A", [], ""),
        ("", [], ""),
    ],
)
def test_translate_stop_ambiguous_partial_and_empty_sequences(
    sequence: str, option: list[str], expected: str
) -> None:
    io, stdout, stderr = make_io(f">read desc\n{sequence}\n".encode("ascii"))

    assert convert.run_translate([*option, "-"], io) == 0

    expected_output = f">read desc\n{expected}\n" if expected else ">read desc\n"
    assert stdout.getvalue() == expected_output
    assert stdout.writes == [expected_output]
    assert stderr.getvalue() == ""


def test_translate_mixed_fasta_and_fastq_files_become_ordered_fasta(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, "first.fa", b">first alpha\nATG\n")
    write_fixture(tmp_path, "middle.fq", b"@middle beta\nATGAAA\n+\n!!!!!!\n")
    write_fixture(tmp_path, "last.fa", b">last gamma\nATGAAACCC\n")
    io, stdout, stderr = make_io()

    assert convert.run_translate(["first.fa", "middle.fq", "last.fa"], io) == 0

    assert stdout.getvalue() == (
        ">first alpha\nM\n>middle beta\nMK\n>last gamma\nMKP\n"
    )
    assert stdout.writes == [stdout.getvalue()]
    assert stderr.getvalue() == ""


@pytest.mark.parametrize(
    ("option", "expected"),
    [
        ([], ">r\nMKP\n"),
        (["--width=0"], ">r\nMKP\n"),
        (["--width=1"], ">r\nM\nK\nP\n"),
        (["--width", "2"], ">r\nMK\nP\n"),
        (["--width=3"], ">r\nMKP\n"),
        (["--width=60"], ">r\nMKP\n"),
        (["--width=1000"], ">r\nMKP\n"),
    ],
)
def test_translate_supports_width_forms_and_boundaries(option: list[str], expected: str) -> None:
    io, stdout, stderr = make_io(b">r\nATGAAACCC\n")

    assert convert.run_translate([*option, "-"], io) == 0

    assert stdout.getvalue() == expected
    assert stdout.writes == [expected]
    assert stderr.getvalue() == ""


@pytest.mark.parametrize(
    ("sequence", "expected_lines"),
    [
        (
            (
                "GCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT"
                "GCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT"
                "GCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT"
            ),
            [
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            ],
        ),
        (
            (
                "GCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT"
                "GCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT"
                "GCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT"
                "GCT"
            ),
            [
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                "A",
            ],
        ),
        (
            (
                "GCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT"
                "GCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT"
                "GCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT"
                "GCTGCTGCT"
            ),
            [
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                "AAA",
            ],
        ),
    ],
)
def test_translate_wraps_protein_lengths_sixty_sixty_one_and_sixty_three(
    sequence: str, expected_lines: list[str]
) -> None:
    io, stdout, _ = make_io(f">protein\n{sequence}\n".encode("ascii"))

    assert convert.run_translate(["--width=60", "-"], io) == 0

    assert stdout.getvalue() == ">protein\n" + "\n".join(expected_lines) + "\n"
    assert stdout.writes == [stdout.getvalue()]


@pytest.mark.parametrize(
    ("runner", "argv"),
    [
        (convert.run_revcomp, []),
        (convert.run_translate, []),
    ],
)
def test_conversions_report_missing_file(runner: Callable[[list[str], Io], int], argv: list[str]) -> None:
    io, _, _ = make_io()

    assert_failure(runner, argv, io, "missing FILE")


@pytest.mark.parametrize(
    "runner",
    [convert.run_revcomp, convert.run_translate],
)
def test_conversions_reject_duplicate_stdin_before_reading(
    runner: Callable[[list[str], Io], int],
) -> None:
    io, _, _ = make_io(b">s\nAC\n")

    assert_failure(runner, ["missing.fa", "-", "-"], io, "standard input given twice")
    assert io.stdin.tell() == 0


@pytest.mark.parametrize(
    ("runner", "option"),
    [
        (convert.run_revcomp, "-x"),
        (convert.run_revcomp, "--"),
        (convert.run_revcomp, "--help"),
        (convert.run_translate, "-x"),
        (convert.run_translate, "--"),
        (convert.run_translate, "--help"),
    ],
)
def test_conversions_safely_quote_unknown_options(
    runner: Callable[[list[str], Io], int], option: str
) -> None:
    io, _, _ = make_io()

    assert_failure(runner, [option], io, f'unknown option "{option}"')


def test_translate_flag_assignment_is_an_unknown_option() -> None:
    io, _, _ = make_io()

    assert_failure(
        convert.run_translate,
        ["--to-stop=yes"],
        io,
        'unknown option "--to-stop=yes"',
    )


@pytest.mark.parametrize(
    ("runner", "argv", "message"),
    [
        (convert.run_revcomp, ["--width"], "--width needs a value"),
        (convert.run_translate, ["--width"], "--width needs a value"),
        (convert.run_translate, ["--frame"], "--frame needs a value"),
    ],
)
def test_conversions_report_missing_option_values(
    runner: Callable[[list[str], Io], int], argv: list[str], message: str
) -> None:
    io, _, _ = make_io()

    assert_failure(runner, argv, io, message)


@pytest.mark.parametrize(
    ("raw_width", "shown_width"),
    [
        ("-0", "-0"),
        ("+3", "+3"),
        ("007", "007"),
        ("1001", "1001"),
        ("1000000000", "1000000000"),
    ],
)
@pytest.mark.parametrize("runner", [convert.run_revcomp, convert.run_translate])
def test_conversions_reject_invalid_width_lexemes_and_bounds(
    runner: Callable[[list[str], Io], int], raw_width: str, shown_width: str
) -> None:
    io, _, _ = make_io()

    assert_failure(
        runner,
        [f"--width={raw_width}", "missing.fa"],
        io,
        f'invalid --width "{shown_width}" (expected an integer from 0 to 1000)',
    )


def test_translate_rejects_invalid_frame_with_exact_choice_diagnostic() -> None:
    io, _, _ = make_io()

    assert_failure(
        convert.run_translate,
        ["--frame", "4", "missing.fa"],
        io,
        'invalid --frame "4" (expected 1, 2, 3, -1, -2 or -3)',
    )


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        (
            ["--width=1001", "--frame=4", "missing.fa"],
            'invalid --frame "4" (expected 1, 2, 3, -1, -2 or -3)',
        ),
        (
            ["--frame=4", "--width=1001", "missing.fa"],
            'invalid --frame "4" (expected 1, 2, 3, -1, -2 or -3)',
        ),
    ],
)
def test_translate_validates_frame_before_width_independent_of_argv_order(
    argv: list[str], message: str
) -> None:
    io, _, _ = make_io()

    assert_failure(convert.run_translate, argv, io, message)


def test_translate_finds_syntax_errors_before_deferred_choice_errors() -> None:
    io, _, _ = make_io()

    assert_failure(
        convert.run_translate,
        ["--frame=4", "--unknown"],
        io,
        'unknown option "--unknown"',
    )


@pytest.mark.parametrize(
    ("argv", "expected"),
    [
        (["--width", "2", "-"], ">r\nNA\nCG\nT\n"),
        (["--width=2", "-"], ">r\nNA\nCG\nT\n"),
        (["--width=3", "--width", "2", "-"], ">r\nNA\nCG\nT\n"),
    ],
)
def test_revcomp_accepts_both_value_syntaxes_and_last_value_wins(
    argv: list[str], expected: str
) -> None:
    io, stdout, _ = make_io(b">r\nACGTN\n")

    assert convert.run_revcomp(argv, io) == 0

    assert stdout.getvalue() == expected
    assert stdout.writes == [expected]


def test_translate_accepts_interspersed_files_negative_frame_and_last_value_wins(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, "first.fa", b">first\nATGAAACCC\n")
    write_fixture(tmp_path, "second.fa", b">second\nATGAAACCC\n")
    io, stdout, stderr = make_io()

    assert convert.run_translate(
        ["first.fa", "--frame=-1", "--width=3", "second.fa", "--frame", "-2"], io
    ) == 0

    assert stdout.getvalue() == ">first\nGF\n>second\nGF\n"
    assert stdout.writes == [stdout.getvalue()]
    assert stderr.getvalue() == ""


@pytest.mark.parametrize(
    ("runner", "failure", "message"),
    [
        (convert.run_revcomp, "missing", "missing.fa: cannot read the file"),
        (convert.run_revcomp, "format", 'format.fa:2: invalid character "R" in sequence'),
        (convert.run_revcomp, "utf8", "utf8.fa: input is not valid UTF-8"),
        (convert.run_translate, "missing", "missing.fa: cannot read the file"),
        (convert.run_translate, "format", 'format.fa:2: invalid character "R" in sequence'),
        (convert.run_translate, "utf8", "utf8.fa: input is not valid UTF-8"),
    ],
)
def test_conversions_leave_output_untouched_when_a_later_file_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    runner: Callable[[list[str], Io], int],
    failure: str,
    message: str,
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, "good.fa", b">good\nATG\n")
    if failure == "format":
        write_fixture(tmp_path, "format.fa", b">bad\nACGR\n")
        argv = ["good.fa", "format.fa", "missing-after-error.fa"]
    elif failure == "utf8":
        write_fixture(tmp_path, "utf8.fa", b">bad\nAC\xff\n")
        argv = ["good.fa", "utf8.fa", "missing-after-error.fa"]
    else:
        argv = ["good.fa", "missing.fa"]
    io, _, _ = make_io()

    assert_failure(runner, argv, io, message)


@pytest.mark.parametrize("runner", [convert.run_revcomp, convert.run_translate])
def test_conversions_decode_exact_malformed_stdin_bytes_atomically(
    runner: Callable[[list[str], Io], int],
) -> None:
    io, _, _ = make_io(b">s\nAC\xff\n")

    assert_failure(runner, ["-"], io, "<stdin>: input is not valid UTF-8")


@pytest.mark.parametrize("runner", [convert.run_revcomp, convert.run_translate])
def test_conversions_decode_exact_malformed_file_bytes_atomically(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, runner: Callable[[list[str], Io], int]
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, "bad.fa", b">s\nAC\xff\n")
    io, _, _ = make_io()

    assert_failure(runner, ["bad.fa"], io, "bad.fa: input is not valid UTF-8")


@pytest.mark.parametrize(
    ("runner", "expected"),
    [
        (convert.run_revcomp, ">r desc\nNNACGT\n"),
        (convert.run_translate, ">r desc\nTX\n"),
    ],
)
def test_conversions_normalize_bom_crlf_lowercase_and_missing_final_lf(
    runner: Callable[[list[str], Io], int], expected: str
) -> None:
    io, stdout, stderr = make_io(b"\xef\xbb\xbf>r desc\r\nacgt\r\nnn")

    assert runner(["-"], io) == 0

    assert stdout.getvalue() == expected
    assert stdout.writes == [expected]
    assert stderr.getvalue() == ""


@pytest.mark.parametrize(
    "runner",
    [convert.run_revcomp, convert.run_translate],
)
def test_empty_conversions_still_make_exactly_one_empty_write(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, runner: Callable[[list[str], Io], int]
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, "empty.fa", b"")
    io, stdout, stderr = make_io()

    assert runner(["empty.fa"], io) == 0

    assert stdout.getvalue() == ""
    assert stdout.writes == [""]
    assert stderr.getvalue() == ""
    assert stderr.writes == []

import ast
from inspect import signature
from io import BytesIO, StringIO
from pathlib import Path
from types import SimpleNamespace

import pytest

from seqmill import cli
from seqmill.tools import align, convert, digest, motif, orfs, qc, stats

EXPECTED_USAGE = (
    "usage: seqmill stats FILE...\n"
    "       seqmill revcomp [--width N] FILE...\n"
    "       seqmill translate [--frame F] [--to-stop] [--width N] FILE...\n"
    "       seqmill orfs [--min N] FILE...\n"
    "       seqmill digest --enzyme NAME[,NAME]... [--circular] FILE...\n"
    "       seqmill align [--mode global|local] [--match N] [--mismatch N] [--gap N] [--width N] A B\n"
    "       seqmill motif [--mismatches N] [--both-strands] PATTERN FILE...\n"
    "       seqmill qc [--trim-end Q] [--min-length N] [--min-quality Q] FILE...\n"
)


_LONG_FASTA = b""">long
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
A
"""

class RecordingText(StringIO):
    def __init__(self) -> None:
        super().__init__()
        self.writes: list[str] = []

    def write(self, value: str) -> int:
        self.writes.append(value)
        return super().write(value)


class FalseyBytesIO(BytesIO):
    def __bool__(self) -> bool:
        return False


class FalseyText(RecordingText):
    def __bool__(self) -> bool:
        return False


def invoke(
    argv: list[str], stdin_data: bytes = b""
) -> tuple[int, RecordingText, RecordingText]:
    stdout = RecordingText()
    stderr = RecordingText()
    result = cli.main(argv, stdin=BytesIO(stdin_data), stdout=stdout, stderr=stderr)
    return result, stdout, stderr


def write_fixture(tmp_path: Path, name: str, data: bytes) -> str:
    (tmp_path / name).write_bytes(data)
    return name


def assert_cli_error(argv: list[str], message: str) -> None:
    result, stdout, stderr = invoke(argv)
    expected_stderr = f"seqmill: {message}\n"
    assert result == 2
    assert stdout.getvalue() == ""
    assert stdout.writes == []
    assert stderr.getvalue() == expected_stderr
    assert stderr.writes == [expected_stderr]


def assert_cli_success(
    argv: list[str], expected_stdout: str, expected_stderr: str = ""
) -> None:
    result, stdout, stderr = invoke(argv)
    assert result == 0
    assert stdout.getvalue() == expected_stdout
    assert stderr.getvalue() == expected_stderr


def test_cli_api_usage_and_single_dispatch_dictionary_are_exact() -> None:
    main_signature = signature(cli.main)
    assert list(main_signature.parameters) == ["argv", "stdin", "stdout", "stderr"]
    assert main_signature.parameters["argv"].annotation == list[str]
    assert all(
        main_signature.parameters[name].default is None
        for name in ("stdin", "stdout", "stderr")
    )
    assert main_signature.return_annotation is int
    assert cli.USAGE == EXPECTED_USAGE

    expected_commands = {
        "stats": stats.run,
        "revcomp": convert.run_revcomp,
        "translate": convert.run_translate,
        "orfs": orfs.run,
        "digest": digest.run,
        "align": align.run,
        "motif": motif.run,
        "qc": qc.run,
    }
    assert type(cli._COMMANDS) is dict
    assert cli._COMMANDS == expected_commands

    cli_source = (Path(__file__).resolve().parents[1] / "seqmill" / "cli.py").read_text(
        encoding="utf-8"
    )
    tree = ast.parse(cli_source)
    assert sum(isinstance(node, ast.Dict) for node in ast.walk(tree)) == 1
    imports = {
        node.module
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom)
    }
    assert {
        "seqmill.tools.align",
        "seqmill.tools.convert",
        "seqmill.tools.digest",
        "seqmill.tools.motif",
        "seqmill.tools.orfs",
        "seqmill.tools.qc",
        "seqmill.tools.stats",
    } <= imports


def test_cli_passes_trailing_arguments_and_the_selected_io_unchanged(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen: list[tuple[list[str], object]] = []

    def capture(argv: list[str], io: object) -> int:
        seen.append((argv, io))
        return 1

    monkeypatch.setitem(cli._COMMANDS, "stats", capture)
    stdin = BytesIO(b"unused")
    stdout = RecordingText()
    stderr = RecordingText()
    arguments = ["stats", "--odd=value", "-", "tail"]

    assert cli.main(arguments, stdin=stdin, stdout=stdout, stderr=stderr) == 1
    assert seen[0][0] == ["--odd=value", "-", "tail"]
    assert seen[0][1].stdin is stdin
    assert seen[0][1].stdout is stdout
    assert seen[0][1].stderr is stderr
    assert stdout.getvalue() == ""
    assert stderr.getvalue() == ""


def test_process_entry_configures_utf8_lf_without_decoding_stdin() -> None:
    entry_source = (
        Path(__file__).resolve().parents[1] / "seqmill" / "__main__.py"
    ).read_text(encoding="utf-8")
    assert entry_source == (
        "import sys\n\n"
        "from seqmill.cli import main\n\n"
        'sys.stdout.reconfigure(encoding="utf-8", newline="\\n")\n'
        'sys.stderr.reconfigure(encoding="utf-8", newline="\\n")\n'
        "sys.exit(main(sys.argv[1:]))\n"
    )


def test_application_sources_keep_the_cli_deterministic_and_dependency_free() -> None:
    source_root = Path(__file__).resolve().parents[1] / "seqmill"
    forbidden_imports = {
        "argparse",
        "getopt",
        "http",
        "locale",
        "random",
        "requests",
        "socket",
        "subprocess",
        "time",
        "urllib",
    }
    for source_path in source_root.rglob("*.py"):
        tree = ast.parse(source_path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported = {alias.name.split(".", 1)[0] for alias in node.names}
                assert imported.isdisjoint(forbidden_imports)
            elif isinstance(node, ast.ImportFrom) and node.module is not None:
                assert node.module.split(".", 1)[0] not in forbidden_imports
            if isinstance(node, ast.Attribute) and node.attr in {"environ", "getenv"}:
                assert not (isinstance(node.value, ast.Name) and node.value.id == "os")


def test_help_is_exact_and_only_the_single_argument_is_help() -> None:
    result, stdout, stderr = invoke(["--help"])
    assert result == 0
    assert stdout.getvalue() == EXPECTED_USAGE
    assert stderr.getvalue() == ""
    assert stdout.writes == [EXPECTED_USAGE]
    assert stderr.writes == []


def test_missing_command_has_exact_prefixed_diagnostic() -> None:
    assert_cli_error(
        [],
        "missing command (stats, revcomp, translate, orfs, digest, align, motif or qc)",
    )


def test_unknown_command_has_exact_prefixed_diagnostic() -> None:
    assert_cli_error(["bogus"], 'unknown command "bogus"')


@pytest.mark.parametrize(
    ("argument", "message"),
    [
        ("two words", 'unknown command "two?words"'),
        ("two\twords", 'unknown command "two?words"'),
        ("two\nwords", 'unknown command "two?words"'),
        ("caf\u00e9", 'unknown command "caf?"'),
    ],
)
def test_unknown_command_text_is_safely_quoted_and_kept_to_one_line(
    argument: str, message: str
) -> None:
    assert_cli_error([argument], message)


def test_help_in_any_other_position_is_not_a_help_request() -> None:
    assert_cli_error(["--help", "stats"], 'unknown command "--help"')
    assert_cli_error(["stats", "--help"], 'unknown option "--help"')


def test_none_streams_use_current_sys_streams(monkeypatch: pytest.MonkeyPatch) -> None:
    stdin = BytesIO(b">s\nACGT\n")
    stdout = RecordingText()
    stderr = RecordingText()
    monkeypatch.setattr(cli.sys, "stdin", SimpleNamespace(buffer=stdin))
    monkeypatch.setattr(cli.sys, "stdout", stdout)
    monkeypatch.setattr(cli.sys, "stderr", stderr)

    assert cli.main(["stats", "-"], stdin=None, stdout=None, stderr=None) == 0
    assert stdout.getvalue() == (
        "file     records  bases  min  max  mean  n50   gc%  n\n"
        "<stdin>        1      4    4    4   4.0    4  50.0  0\n"
    )
    assert stderr.getvalue() == ""


def test_each_stream_can_be_overridden_independently(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    default_stdin = BytesIO(b">s\nACGT\n")
    default_stdout = RecordingText()
    default_stderr = RecordingText()
    monkeypatch.setattr(cli.sys, "stdin", SimpleNamespace(buffer=default_stdin))
    monkeypatch.setattr(cli.sys, "stdout", default_stdout)
    monkeypatch.setattr(cli.sys, "stderr", default_stderr)

    explicit_stdout = RecordingText()
    explicit_stderr = RecordingText()
    assert cli.main(["stats", "-"], stdout=explicit_stdout) == 0
    assert explicit_stdout.getvalue() == (
        "file     records  bases  min  max  mean  n50   gc%  n\n"
        "<stdin>        1      4    4    4   4.0    4  50.0  0\n"
    )
    assert default_stdout.getvalue() == ""
    assert default_stderr.getvalue() == ""

    assert cli.main(["--help"], stderr=explicit_stderr) == 0
    assert default_stdout.getvalue() == EXPECTED_USAGE
    assert explicit_stderr.getvalue() == ""

    final_stderr = RecordingText()
    assert cli.main([], stderr=final_stderr) == 2
    assert default_stderr.getvalue() == ""
    assert final_stderr.getvalue() == (
        "seqmill: missing command (stats, revcomp, translate, orfs, digest, align, motif or qc)\n"
    )


def test_falsey_explicit_stream_objects_are_not_replaced() -> None:
    stdin = FalseyBytesIO(b">s\nACGT\n")
    stdout = FalseyText()
    stderr = FalseyText()

    assert cli.main(["stats", "-"], stdin=stdin, stdout=stdout, stderr=stderr) == 0
    assert stdout.getvalue() == (
        "file     records  bases  min  max  mean  n50   gc%  n\n"
        "<stdin>        1      4    4    4   4.0    4  50.0  0\n"
    )
    assert stderr.getvalue() == ""


def test_stats_runs_through_main_and_matches_the_supplied_example(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(
        tmp_path,
        "reads.fa",
        b">seq1 first read\nACGTACGTAC\nGTNNAC\n>seq2\nGGGCCC\n>seq3 empty\n",
    )

    assert_cli_success(
        ["stats", "reads.fa"],
        "file      records  bases  min  max  mean  n50   gc%  n\n"
        "reads.fa        3     22    0   16   7.3   16  65.0  2\n",
    )


def test_revcomp_runs_through_main_and_formats_the_reverse_complement(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, "read.fa", b">read description\nAACGN\n")

    assert_cli_success(
        ["revcomp", "read.fa"],
        ">read description\nNCGTT\n",
    )


def test_translate_runs_through_main_and_matches_the_supplied_example(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(
        tmp_path,
        "gene.fa",
        b">gene1 test construct\n"
        b"GAATTCATGGCCAAAGGATCCTTTGGGAAATAGCCCGGGATGAAACCCTGAGTTTTCATTTGGGAAACAT\n",
    )

    assert_cli_success(
        ["translate", "gene.fa"],
        ">gene1 test construct\nEFMAKGSFGK*PGMKP*VFIWET\n",
    )


def test_orfs_runs_through_main_and_matches_the_supplied_example(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(
        tmp_path,
        "gene.fa",
        b">gene1 test construct\n"
        b"GAATTCATGGCCAAAGGATCCTTTGGGAAATAGCCCGGGATGAAACCCTGAGTTTTCATTTGGGAAACAT\n",
    )

    assert_cli_success(
        ["orfs", "--min", "2", "gene.fa"],
        "name   frame  start  end  aa\n"
        "gene1     +1      7   33   8\n"
        "gene1     +1     40   51   3\n"
        "gene1     -1     56   70   4\n",
    )


def test_digest_runs_through_main_and_matches_the_supplied_example(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(
        tmp_path,
        "gene.fa",
        b">gene1 test construct\n"
        b"GAATTCATGGCCAAAGGATCCTTTGGGAAATAGCCCGGGATGAAACCCTGAGTTTTCATTTGGGAAACAT\n",
    )

    assert_cli_success(
        ["digest", "--enzyme", "EcoRI,BamHI", "gene.fa"],
        "name   fragment  start  end  length\n"
        "gene1         1      1    1       1\n"
        "gene1         2      2   16      15\n"
        "gene1         3     17   70      54\n",
    )


def test_align_runs_through_main_and_matches_the_supplied_global_example(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, "a.fa", b">a\nGATTACAGATTACA\n")
    write_fixture(tmp_path, "b.fa", b">b\nGCATGCAGATTCA\n")

    assert_cli_success(
        ["align", "--width", "10", "a.fa", "b.fa"],
        "score: 5\n"
        "length: 14\n"
        "identity: 10/14 (71.4%)\n"
        "A: a 1-14\n"
        "B: b 1-13\n"
        "\n"
        "GATTACAGAT\n"
        "|xx|x|||||\n"
        "GCATGCAGAT\n"
        "\n"
        "TACA\n"
        "|-||\n"
        "T-CA\n",
    )


def test_motif_runs_through_main_and_matches_the_supplied_example(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(
        tmp_path,
        "gene.fa",
        b">gene1 test construct\n"
        b"GAATTCATGGCCAAAGGATCCTTTGGGAAATAGCCCGGGATGAAACCCTGAGTTTTCATTTGGGAAACAT\n",
    )

    assert_cli_success(
        ["motif", "--both-strands", "--mismatches", "1", "ATGAAR", "gene.fa"],
        "name   strand  start  end  mismatches  match\n"
        "gene1  -           3    8           1  ATTCAT\n"
        "gene1  +          40   45           0  ATGAAA\n"
        "gene1  -          54   59           0  TTTCAT\n",
    )


def test_qc_runs_through_main_and_matches_the_supplied_output_and_summary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(
        tmp_path,
        "reads.fq",
        b"@r1 lane1\nACGTACGT\n+\nIIIIII##\n"
        b"@r2\nGGGGCCCC\n+\n5555!!!!\n"
        b"@r3\nACGT\n+\nII##\n"
        b"@r4\nTTTTAAAA\n+\n++++++++\n"
        b"@r5\nCCCCGGGG\n+\nIIIIIIII\n",
    )

    assert_cli_success(
        ["qc", "--trim-end", "20", "--min-length", "4", "--min-quality", "25", "reads.fq"],
        "@r1 lane1\nACGTAC\n+\nIIIIII\n@r5\nCCCCGGGG\n+\nIIIIIIII\n",
        "kept 2 of 5 reads: 2 too short, 1 low quality, 16 bases trimmed\n",
    )


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        (
            ["revcomp", "--width", "1001", "--unknown", "missing.fa"],
            'unknown option "--unknown"',
        ),
        (
            ["revcomp", "--width", "1001", "--width"],
            "--width needs a value",
        ),
    ],
)
def test_option_syntax_errors_precede_deferred_value_and_file_errors(
    argv: list[str], message: str
) -> None:
    assert_cli_error(argv, message)


def test_option_values_are_validated_in_usage_order_before_file_reads() -> None:
    assert_cli_error(
        ["translate", "--width", "1001", "--frame", "4", "missing.fa"],
        'invalid --frame "4" (expected 1, 2, 3, -1, -2 or -3)',
    )


def test_orfs_minimum_zero_reports_the_supplied_worked_error() -> None:
    assert_cli_error(
        ["orfs", "--min", "0", "reads.fa"],
        'invalid --min "0" (expected an integer from 1 to 100000)',
    )


def test_repeated_invalid_option_uses_the_last_value() -> None:
    assert_cli_error(
        ["revcomp", "--width", "1001", "--width", "1002", "missing.fa"],
        'invalid --width "1002" (expected an integer from 0 to 1000)',
    )


def test_repeated_width_uses_last_value_and_can_succeed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, "repeat.fa", b">s\nAACGN\n")

    assert_cli_success(
        ["revcomp", "--width", "1001", "--width", "1", "repeat.fa"],
        ">s\nN\nC\nG\nT\nT\n",
    )


def test_invalid_option_is_reported_before_a_missing_file() -> None:
    assert_cli_error(
        ["revcomp", "--width", "1001", "missing.fa"],
        'invalid --width "1001" (expected an integer from 0 to 1000)',
    )


def test_digest_enzyme_checks_precede_file_checks() -> None:
    assert_cli_error(["digest", "missing.fa"], "missing --enzyme")
    assert_cli_error(
        ["digest", "--enzyme", "not-an-enzyme", "missing.fa"],
        'unknown enzyme "not-an-enzyme" (expected ApaI, BamHI, EcoRI, EcoRV, HindIII, KpnI, '
        "NcoI, NotI, PstI, SacI, SalI, SmaI, XbaI or XhoI)",
    )


def test_motif_pattern_checks_precede_file_checks() -> None:
    assert_cli_error(["motif"], "missing PATTERN")
    assert_cli_error(
        ["motif", "Z", "missing.fa"],
        'invalid pattern "Z" (expected 1 to 64 letters from ACGTNRYSWKMBDHV)',
    )
    assert_cli_error(["motif", "ACGT"], "missing FILE")


def test_file_count_and_duplicate_stdin_checks_precede_file_reads() -> None:
    assert_cli_error(["stats"], "missing FILE")
    assert_cli_error(
        ["stats", "missing.fa", "-", "-"],
        "standard input given twice",
    )


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        (["align"], "missing A"),
        (["align", "missing-a.fa"], "missing B"),
        (["align", "missing-a.fa", "missing-b.fa", "extra"], "too many arguments"),
        (["align", "-", "-"], "standard input given twice"),
    ],
)
def test_align_argument_errors_precede_file_reads(argv: list[str], message: str) -> None:
    assert_cli_error(argv, message)


@pytest.mark.parametrize(
    ("failure", "path", "data", "message"),
    [
        ("read", "missing-first.fa", None, "missing-first.fa: cannot read the file"),
        (
            "decode",
            "invalid-first.fa",
            b">bad\nAC\xff\n",
            "invalid-first.fa: input is not valid UTF-8",
        ),
        (
            "format",
            "malformed-first.fa",
            b">bad\nACGR\n",
            'malformed-first.fa:2: invalid character "R" in sequence',
        ),
    ],
)
def test_first_file_read_decode_or_format_error_wins(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    failure: str,
    path: str,
    data: bytes | None,
    message: str,
) -> None:
    monkeypatch.chdir(tmp_path)
    if data is not None:
        write_fixture(tmp_path, path, data)

    assert_cli_error(["stats", path, "later-missing.fa"], message)


@pytest.mark.parametrize(
    ("command", "first_path", "first_data", "later_path", "later_data", "argv", "message"),
    [
        (
            "stats",
            "good.fa",
            b">good\nAC\n",
            "missing-later.fa",
            None,
            ["stats", "good.fa", "missing-later.fa"],
            "missing-later.fa: cannot read the file",
        ),
        (
            "revcomp",
            "good.fa",
            b">good\nAC\n",
            "invalid-later.fa",
            b">bad\nAC\xff\n",
            ["revcomp", "good.fa", "invalid-later.fa"],
            "invalid-later.fa: input is not valid UTF-8",
        ),
        (
            "translate",
            "good.fa",
            b">good\nAC\n",
            "malformed-later.fa",
            b">bad\nACGR\n",
            ["translate", "good.fa", "malformed-later.fa"],
            'malformed-later.fa:2: invalid character "R" in sequence',
        ),
        (
            "orfs",
            "good.fa",
            b">good\nAC\n",
            "missing-later.fa",
            None,
            ["orfs", "good.fa", "missing-later.fa"],
            "missing-later.fa: cannot read the file",
        ),
        (
            "digest",
            "good.fa",
            b">good\nAC\n",
            "invalid-later.fa",
            b">bad\nAC\xff\n",
            ["digest", "--enzyme", "EcoRI", "good.fa", "invalid-later.fa"],
            "invalid-later.fa: input is not valid UTF-8",
        ),
        (
            "align",
            "a.fa",
            b">a\nAC\n",
            "b.fa",
            b">b\nACGR\n",
            ["align", "a.fa", "b.fa"],
            'b.fa:2: invalid character "R" in sequence',
        ),
        (
            "motif",
            "good.fa",
            b">good\nAC\n",
            "malformed-later.fa",
            b">bad\nACGR\n",
            ["motif", "AC", "good.fa", "malformed-later.fa"],
            'malformed-later.fa:2: invalid character "R" in sequence',
        ),
        (
            "qc",
            "good.fq",
            b"@good\nAC\n+\n!!\n",
            "malformed-later.fq",
            b"@bad\nAC\n+\n!\n",
            ["qc", "good.fq", "malformed-later.fq"],
            "malformed-later.fq:4: sequence and quality have different lengths",
        ),
    ],
)
def test_later_file_failure_is_atomic_for_every_command(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    command: str,
    first_path: str,
    first_data: bytes,
    later_path: str,
    later_data: bytes | None,
    argv: list[str],
    message: str,
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, first_path, first_data)
    if later_data is not None:
        write_fixture(tmp_path, later_path, later_data)

    result, stdout, stderr = invoke(argv)
    expected_stderr = f"seqmill: {message}\n"
    assert result == 2
    assert stdout.getvalue() == ""
    assert stdout.writes == []
    assert stderr.getvalue() == expected_stderr
    assert stderr.writes == [expected_stderr]


def test_orfs_without_qualifying_rows_returns_one_without_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, "no-orfs.fa", b">s\nACGTACGT\n")

    result, stdout, stderr = invoke(["orfs", "no-orfs.fa"])
    assert result == 1
    assert stdout.getvalue() == ""
    assert stderr.getvalue() == ""


def test_digest_of_an_empty_file_returns_one_without_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, "empty.fa", b"")

    result, stdout, stderr = invoke(["digest", "--enzyme", "EcoRI", "empty.fa"])
    assert result == 1
    assert stdout.getvalue() == ""
    assert stderr.getvalue() == ""


def test_motif_without_matches_returns_one_without_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, "sequence.fa", b">s\nACGT\n")

    result, stdout, stderr = invoke(["motif", "AAA", "sequence.fa"])
    assert result == 1
    assert stdout.getvalue() == ""
    assert stderr.getvalue() == ""


def test_local_alignment_without_a_positive_score_has_exact_negative_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, "a.fa", b">a\nA\n")
    write_fixture(tmp_path, "b.fa", b">b\nT\n")

    result, stdout, stderr = invoke(["align", "--mode", "local", "a.fa", "b.fa"])
    assert result == 1
    assert stdout.getvalue() == "score: 0\nno alignment\n"
    assert stderr.getvalue() == ""


@pytest.mark.parametrize(
    ("path", "data", "message"),
    [
        ("empty.fa", b"", "empty.fa: no sequences"),
        ("empty-sequence.fa", b">empty\n", 'empty-sequence.fa: sequence "empty" is empty'),
        (
            "long.fa",
            _LONG_FASTA,
            'long.fa: sequence "long" is longer than 2000 bases',
        ),
    ],
)
def test_align_validates_a_before_reading_unreadable_b(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    path: str,
    data: bytes,
    message: str,
) -> None:
    monkeypatch.chdir(tmp_path)
    write_fixture(tmp_path, path, data)

    assert_cli_error(["align", path, "missing-b.fa"], message)

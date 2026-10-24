from io import BytesIO, StringIO

import pytest

from seqmill import cli


def invoke(argv):
    stdout = StringIO()
    stderr = StringIO()
    result = cli.main(argv, stdin=BytesIO(b""), stdout=stdout, stderr=stderr)
    return result, stdout.getvalue(), stderr.getvalue()


def test_version_constant_is_the_project_version():
    assert cli.VERSION == "0.1.0"


def test_version_is_printed_with_the_program_name():
    assert invoke(["--version"]) == (0, "seqmill 0.1.0\n", "")


def test_version_does_not_read_standard_input():
    stdout = StringIO()
    stderr = StringIO()
    result = cli.main(["--version"], stdin=None, stdout=stdout, stderr=stderr)
    assert result == 0
    assert stdout.getvalue() == "seqmill 0.1.0\n"
    assert stderr.getvalue() == ""


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        (["--version", "stats"], 'unknown command "--version"'),
        (["--version", "--help"], 'unknown command "--version"'),
        (["--help", "--version"], 'unknown command "--help"'),
        (["--version", ""], 'unknown command "--version"'),
        (["stats", "--version"], 'unknown option "--version"'),
        (["-version"], 'unknown command "-version"'),
        (["--Version"], 'unknown command "--Version"'),
        (["--version=1"], 'unknown command "--version=1"'),
        (["-V"], 'unknown command "-V"'),
    ],
)
def test_version_is_only_recognized_as_the_single_argument(argv, message):
    result, stdout, stderr = invoke(argv)
    assert result == 2
    assert stdout == ""
    assert stderr == f"seqmill: {message}\n"


def test_help_text_is_unchanged_and_does_not_mention_the_version():
    result, stdout, stderr = invoke(["--help"])
    assert result == 0
    assert stdout == cli.USAGE
    assert "--version" not in cli.USAGE
    assert stderr == ""


def test_missing_command_message_is_unchanged():
    assert invoke([]) == (
        2,
        "",
        "seqmill: missing command (stats, revcomp, translate, orfs, digest, align, motif or qc)\n",
    )

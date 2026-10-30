import re
import subprocess
import sys
from io import BytesIO, StringIO
from pathlib import Path

import pytest

from seqmill.cli import main


def test_main_reports_malformed_utf8_from_real_standard_input() -> None:
    repository_root = Path(__file__).resolve().parents[1]

    result = subprocess.run(
        [sys.executable, "-m", "seqmill", "stats", "-"],
        cwd=repository_root,
        input=b">s\nAC\xff\n",
        capture_output=True,
        timeout=30,
        check=False,
    )

    assert result.returncode == 2
    assert result.stdout == b""
    assert result.stderr == b"seqmill: <stdin>: input is not valid UTF-8\n"


def test_readme_documents_and_runs_the_stats_example(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    readme_path = Path(__file__).resolve().parents[1] / "README.md"
    readme = readme_path.read_text(encoding="utf-8")

    expected_usage = (
        "usage: seqmill stats [--no-header] FILE...\n"
        "       seqmill revcomp [--width N] FILE...\n"
        "       seqmill translate [--frame F] [--to-stop] [--width N] FILE...\n"
        "       seqmill orfs [--min N] [--strand both|forward|reverse] FILE...\n"
        "       seqmill digest --enzyme NAME[,NAME]... [--circular] FILE...\n"
        "       seqmill align [--mode global|local] [--match N] [--mismatch N] [--gap N] [--width N] A B\n"
        "       seqmill motif [--mismatches N] [--both-strands] PATTERN FILE...\n"
        "       seqmill qc [--trim-end Q] [--min-length N] [--min-quality Q] FILE...\n"
    )
    expected_fixture = (
        ">seq1 first read\n"
        "ACGTACGTAC\n"
        "GTNNAC\n"
        ">seq2\n"
        "GGGCCC\n"
        ">seq3 empty\n"
    )
    expected_table = (
        "file      records  bases  min  max  mean  n50   gc%  n\n"
        "reads.fa        3     22    0   16   7.3   16  65.0  2\n"
    )

    assert [line for line in readme.splitlines() if line.startswith("## ")] == [
        "## Install",
        "## Usage",
        "## Commands",
        "## Exit codes",
        "## Development",
    ]
    assert "## Install\n\nThere is nothing to install. Use Python 3.11 or newer and run `python -m seqmill` from the repository root." in readme
    assert "## Usage\n\nRun the commands below as `python -m seqmill COMMAND ...` from the repository root. File arguments are paths or `-` for standard input." in readme
    assert f"```text\n{expected_usage.rstrip()}\n```" in readme
    assert "For example, save this supplied fixture as `reads.fa` in the repository root:" in readme
    assert "Run `python -m seqmill stats reads.fa` from the repository root. Its output is:" in readme

    code_blocks = re.findall(r"```([^\n]*)\n(.*?)\n```", readme, re.DOTALL)
    assert code_blocks == [
        ("text", expected_usage.rstrip()),
        ("text", expected_fixture.rstrip()),
        ("text", expected_table.rstrip()),
        ("sh", "python -m ruff check .\npython -m pytest -q"),
    ]

    expected_commands = (
        "- `stats [--no-header] FILE...` reports each file's record count, total bases, minimum and maximum sequence lengths, mean length, N50, GC percentage and number of `N` bases. With `--no-header` the header line is left out and the rows keep their column widths.\n"
        "- `revcomp [--width N] FILE...` writes reverse-complemented records. FASTA and FASTQ formats are preserved, and FASTQ quality strings are reversed along with their sequences.\n"
        "- `translate [--frame F] [--to-stop] [--width N] FILE...` writes translated proteins as FASTA. Choose any of the six reading frames; negative frames use the reverse complement, and `--to-stop` truncates translation before the first stop codon.\n"
        "- `orfs [--min N] [--strand both|forward|reverse] FILE...` reports complete open reading frames from all six frames, or from only the three forward or the three reverse frames with `--strand`. `--min` sets the minimum amino-acid length, excluding the stop codon.\n"
        "- `digest --enzyme NAME[,NAME]... [--circular] FILE...` reports restriction fragments for one or more named enzymes. Use `--circular` to digest sequences as rings.\n"
        "- `align [--mode global|local] [--match N] [--mismatch N] [--gap N] [--width N] A B` aligns the first record from each input, globally or locally, with configurable scores and output width.\n"
        "- `motif [--mismatches N] [--both-strands] PATTERN FILE...` searches for IUPAC DNA motifs, optionally allowing mismatches and checking both strands.\n"
        "- `qc [--trim-end Q] [--min-length N] [--min-quality Q] FILE...` trims and filters FASTQ reads by end quality, length and mean quality, then reports read and trimming counts on standard error."
    )
    commands_section = readme.split("## Commands\n\n", 1)[1].split("\n## Exit codes\n", 1)[0]
    assert commands_section == f"{expected_commands}\n"
    assert (
        "seqmill exits with 0 on success, 1 for a negative answer (a table command with no rows or a local alignment with no alignment), and 2 for an error. Errors are written to standard error and leave standard output empty."
        in readme
    )
    assert "## Development\n\nFrom the repository root, run the checks with:" in readme

    help_stdout = StringIO()
    help_stderr = StringIO()
    assert main(["--help"], stdin=BytesIO(), stdout=help_stdout, stderr=help_stderr) == 0
    assert help_stdout.getvalue() == expected_usage
    assert help_stderr.getvalue() == ""

    monkeypatch.chdir(tmp_path)
    (tmp_path / "reads.fa").write_bytes(expected_fixture.encode("utf-8"))
    stats_stdout = StringIO()
    stats_stderr = StringIO()
    assert (
        main(
            ["stats", "reads.fa"],
            stdin=BytesIO(),
            stdout=stats_stdout,
            stderr=stats_stderr,
        )
        == 0
    )
    assert stats_stdout.getvalue() == expected_table
    assert stats_stderr.getvalue() == ""

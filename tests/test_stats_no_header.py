from io import BytesIO, StringIO
from pathlib import Path

import pytest

from seqmill.cli import main

READS = ">seq1 first read\nACGTACGTAC\nGTNNAC\n>seq2\nGGGCCC\n>seq3 empty\n"
READS_ROW = "reads.fa        3     22    0   16   7.3   16  65.0  2\n"
SHORT_ROW = "a.fa        1      4    4    4   4.0    4  50.0  0\n"
HEADER = "file      records  bases  min  max  mean  n50   gc%  n\n"


@pytest.mark.parametrize(
    ("argv", "stdin", "status", "stdout", "stderr"),
    [
        (["stats", "reads.fa"], b"", 0, HEADER + READS_ROW, ""),
        (["stats", "--no-header", "reads.fa"], b"", 0, READS_ROW, ""),
        (["stats", "reads.fa", "--no-header"], b"", 0, READS_ROW, ""),
        (["stats", "--no-header", "--no-header", "reads.fa"], b"", 0, READS_ROW, ""),
        (["stats", "--no-header", "a.fa"], b"", 0, SHORT_ROW, ""),
        (
            ["stats", "--no-header", "reads.fa", "a.fa"],
            b"",
            0,
            READS_ROW + "a.fa            1      4    4    4   4.0    4  50.0  0\n",
            "",
        ),
        (
            ["stats", "--no-header", "empty.fa"],
            b"",
            0,
            "empty.fa        0      0    -    -     -    -    -  0\n",
            "",
        ),
        (
            ["stats", "--no-header", "-"],
            b">s\nACGT\n",
            0,
            "<stdin>        1      4    4    4   4.0    4  50.0  0\n",
            "",
        ),
        (["stats", "--no-header"], b"", 2, "", "seqmill: missing FILE\n"),
        (
            ["stats", "--no-header=1", "reads.fa"],
            b"",
            2,
            "",
            'seqmill: unknown option "--no-header=1"\n',
        ),
        (["stats", "--header", "reads.fa"], b"", 2, "", 'seqmill: unknown option "--header"\n'),
        (
            ["revcomp", "--no-header", "reads.fa"],
            b"",
            2,
            "",
            'seqmill: unknown option "--no-header"\n',
        ),
        (
            ["stats", "--no-header", "nope.fa"],
            b"",
            2,
            "",
            "seqmill: nope.fa: cannot read the file\n",
        ),
        (["stats", "--no-header", "--", "reads.fa"], b"", 2, "", 'seqmill: unknown option "--"\n'),
    ],
)
def test_stats_no_header_prints_the_rows_without_the_header_line(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    argv: list[str],
    stdin: bytes,
    status: int,
    stdout: str,
    stderr: str,
) -> None:
    (tmp_path / "reads.fa").write_bytes(READS.encode("utf-8"))
    (tmp_path / "a.fa").write_bytes(b">a\nACGT\n")
    (tmp_path / "empty.fa").write_bytes(b"")
    monkeypatch.chdir(tmp_path)
    out = StringIO()
    err = StringIO()

    assert main(argv, stdin=BytesIO(stdin), stdout=out, stderr=err) == status
    assert out.getvalue() == stdout
    assert err.getvalue() == stderr

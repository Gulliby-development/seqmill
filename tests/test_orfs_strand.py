from io import BytesIO, StringIO

from seqmill.cli import main as cli_main
from seqmill.formats.files import Io
from seqmill.tools import orfs

GENE_FASTA = b">gene1 test construct\nGAATTCATGGCCAAAGGATCCTTTGGGAAATAGCCCGGGATGAAACCCTGAGTTTTCATTTGGGAAACAT\n"


def make_io(data: bytes = b"") -> Io:
    return Io(BytesIO(data), StringIO(), StringIO())


def invoke_cli(argv: list[str], data: bytes = b"") -> tuple[int, Io]:
    streams = make_io(data)
    status = cli_main(argv, streams.stdin, streams.stdout, streams.stderr)
    return status, streams


def assert_lf_only(streams: Io) -> None:
    assert "\r" not in streams.stdout.getvalue()


def test_default_and_both_match_the_full_six_frame_table() -> None:
    expected = "name   frame  start  end  aa\ngene1     +1      7   33   8\ngene1     +1     40   51   3\ngene1     -1     56   70   4\n"
    streams = make_io(GENE_FASTA)

    assert orfs.run(["--min", "2", "-"], streams) == 0
    assert streams.stdout.getvalue() == expected
    assert streams.stderr.getvalue() == ""
    assert_lf_only(streams)

    status, cli_streams = invoke_cli(
        ["orfs", "--min", "2", "--strand", "both", "-"], GENE_FASTA
    )
    assert status == 0
    assert cli_streams.stdout.getvalue() == expected
    assert cli_streams.stderr.getvalue() == ""
    assert_lf_only(cli_streams)


def test_forward_keeps_only_forward_rows_and_uses_selected_widths() -> None:
    gene_streams = make_io(GENE_FASTA)

    assert orfs.run(["--min", "2", "--strand", "forward", "-"], gene_streams) == 0
    assert gene_streams.stdout.getvalue() == (
        "name   frame  start  end  aa\n"
        "gene1     +1      7   33   8\n"
        "gene1     +1     40   51   3\n"
    )
    assert gene_streams.stderr.getvalue() == ""
    assert_lf_only(gene_streams)

    padding_fasta = (
        b">gene1\nATGAAATAA\n"
        b">long_reverse_record\nTTATTTCAT\n"
    )
    both_streams = make_io(padding_fasta)
    assert orfs.run(["--min", "1", "-", "--strand", "both"], both_streams) == 0
    assert both_streams.stdout.getvalue() == (
        "name                 frame  start  end  aa\n"
        "gene1                   +1      1    9   2\n"
        "long_reverse_record     -1      1    9   2\n"
    )
    assert both_streams.stderr.getvalue() == ""
    assert_lf_only(both_streams)

    status, forward_streams = invoke_cli(
        ["orfs", "--min", "1", "-", "--strand", "forward"], padding_fasta
    )
    assert status == 0
    assert forward_streams.stdout.getvalue() == (
        "name   frame  start  end  aa\n"
        "gene1     +1      1    9   2\n"
    )
    assert forward_streams.stderr.getvalue() == ""
    assert_lf_only(forward_streams)


def test_reverse_keeps_only_reverse_rows() -> None:
    status, streams = invoke_cli(
        ["orfs", "--min", "2", "--strand", "reverse", "-"], GENE_FASTA
    )

    assert status == 0
    assert streams.stdout.getvalue() == (
        "name   frame  start  end  aa\n"
        "gene1     -1     56   70   4\n"
    )
    assert streams.stderr.getvalue() == ""
    assert_lf_only(streams)


def test_equals_form_selects_reverse_strand() -> None:
    streams = make_io(GENE_FASTA)

    assert orfs.run(["--min", "2", "--strand=reverse", "-"], streams) == 0
    assert streams.stdout.getvalue() == (
        "name   frame  start  end  aa\n"
        "gene1     -1     56   70   4\n"
    )
    assert streams.stderr.getvalue() == ""
    assert_lf_only(streams)


def test_strand_option_after_operand_selects_reverse_strand() -> None:
    status, streams = invoke_cli(
        ["orfs", "--min", "2", "-", "--strand", "reverse"], GENE_FASTA
    )

    assert status == 0
    assert streams.stdout.getvalue() == (
        "name   frame  start  end  aa\n"
        "gene1     -1     56   70   4\n"
    )
    assert streams.stderr.getvalue() == ""
    assert_lf_only(streams)


def test_last_strand_value_wins_including_after_invalid_value() -> None:
    reverse_streams = make_io(GENE_FASTA)

    assert (
        orfs.run(
            ["--min", "2", "--strand", "forward", "--strand", "reverse", "-"],
            reverse_streams,
        )
        == 0
    )
    assert reverse_streams.stdout.getvalue() == (
        "name   frame  start  end  aa\n"
        "gene1     -1     56   70   4\n"
    )
    assert reverse_streams.stderr.getvalue() == ""
    assert_lf_only(reverse_streams)

    status, forward_streams = invoke_cli(
        ["orfs", "--min", "2", "--strand", "up", "--strand", "forward", "-"],
        GENE_FASTA,
    )
    assert status == 0
    assert forward_streams.stdout.getvalue() == (
        "name   frame  start  end  aa\n"
        "gene1     +1      7   33   8\n"
        "gene1     +1     40   51   3\n"
    )
    assert forward_streams.stderr.getvalue() == ""
    assert_lf_only(forward_streams)


def test_empty_reverse_selection_returns_one_without_output() -> None:
    status, streams = invoke_cli(
        ["orfs", "--min", "5", "--strand", "reverse", "-"], GENE_FASTA
    )

    assert status == 1
    assert streams.stdout.getvalue() == ""
    assert streams.stderr.getvalue() == ""


def test_invalid_strand_up_has_exact_cli_diagnostic() -> None:
    status, streams = invoke_cli(["orfs", "--strand", "up", "-"], GENE_FASTA)

    assert status == 2
    assert streams.stdout.getvalue() == ""
    assert streams.stderr.getvalue() == (
        'seqmill: invalid --strand "up" (expected both, forward or reverse)\n'
    )


def test_invalid_strand_forward_case_has_exact_cli_diagnostic() -> None:
    status, streams = invoke_cli(["orfs", "--strand", "Forward", "-"], GENE_FASTA)

    assert status == 2
    assert streams.stdout.getvalue() == ""
    assert streams.stderr.getvalue() == (
        'seqmill: invalid --strand "Forward" (expected both, forward or reverse)\n'
    )


def test_empty_strand_equals_value_has_exact_cli_diagnostic() -> None:
    status, streams = invoke_cli(["orfs", "--strand=", "-"], GENE_FASTA)

    assert status == 2
    assert streams.stdout.getvalue() == ""
    assert streams.stderr.getvalue() == (
        'seqmill: invalid --strand "" (expected both, forward or reverse)\n'
    )


def test_missing_strand_value_returns_exact_cli_diagnostic() -> None:
    status, streams = invoke_cli(["orfs", "--strand"], GENE_FASTA)

    assert status == 2
    assert streams.stdout.getvalue() == ""
    assert streams.stderr.getvalue() == "seqmill: --strand needs a value\n"


def test_minimum_validation_precedence_and_cli_help() -> None:
    minimum_error = (
        "seqmill: invalid --min \"0\" (expected an integer from 1 to 100000)\n"
    )
    for arguments in (
        ["orfs", "--min", "0", "--strand", "up", "-"],
        ["orfs", "--strand", "up", "--min", "0", "-"],
    ):
        error_status, error_streams = invoke_cli(arguments, GENE_FASTA)
        assert error_status == 2
        assert error_streams.stdout.getvalue() == ""
        assert error_streams.stderr.getvalue() == minimum_error

    help_status, help_streams = invoke_cli(["--help"])
    assert help_status == 0
    assert (
        "       seqmill orfs [--min N] [--strand both|forward|reverse] FILE...\n"
        in help_streams.stdout.getvalue()
    )
    assert help_streams.stderr.getvalue() == ""
    assert_lf_only(help_streams)

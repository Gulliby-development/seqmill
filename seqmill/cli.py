import sys

from seqmill.core.seq import SeqError, show
from seqmill.formats.files import Io
from seqmill.tools.align import run as align_run
from seqmill.tools.convert import run_revcomp, run_translate
from seqmill.tools.digest import run as digest_run
from seqmill.tools.motif import run as motif_run
from seqmill.tools.orfs import run as orfs_run
from seqmill.tools.qc import run as qc_run
from seqmill.tools.stats import run as stats_run

USAGE = """usage: seqmill stats [--no-header] FILE...
       seqmill revcomp [--width N] FILE...
       seqmill translate [--frame F] [--to-stop] [--width N] FILE...
       seqmill orfs [--min N] [--strand both|forward|reverse] FILE...
       seqmill digest --enzyme NAME[,NAME]... [--circular] FILE...
       seqmill align [--mode global|local] [--match N] [--mismatch N] [--gap N] [--width N] A B
       seqmill motif [--mismatches N] [--both-strands] PATTERN FILE...
       seqmill qc [--trim-end Q] [--min-length N] [--min-quality Q] FILE...
"""

VERSION = "0.1.0"

_COMMANDS = {
    "stats": stats_run,
    "revcomp": run_revcomp,
    "translate": run_translate,
    "orfs": orfs_run,
    "digest": digest_run,
    "align": align_run,
    "motif": motif_run,
    "qc": qc_run,
}


def main(argv: list[str], stdin=None, stdout=None, stderr=None) -> int:
    input_stream = sys.stdin.buffer if stdin is None else stdin
    output_stream = sys.stdout if stdout is None else stdout
    error_stream = sys.stderr if stderr is None else stderr

    if argv == ["--help"]:
        output_stream.write(USAGE)
        return 0
    if argv == ["--version"]:
        output_stream.write(f"seqmill {VERSION}\n")
        return 0
    if not argv:
        error_stream.write(
            "seqmill: missing command "
            "(stats, revcomp, translate, orfs, digest, align, motif or qc)\n"
        )
        return 2

    command = _COMMANDS.get(argv[0])
    if command is None:
        error_stream.write(f"seqmill: unknown command {show(argv[0])}\n")
        return 2

    io = Io(input_stream, output_stream, error_stream)
    try:
        return command(argv[1:], io)
    except SeqError as error:
        error_stream.write(f"seqmill: {error.message}\n")
        return 2

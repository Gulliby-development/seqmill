<p align="center"><img src="docs/logo.svg" width="96" alt="" /></p>

# seqmill

seqmill is a deterministic, offline command-line toolbox for DNA sequence files in FASTA and FASTQ formats. It is for students, lab members and scripts that need quick sequence statistics, conversions, searches, alignments or FASTQ quality filtering without a bioinformatics suite.

## Install

There is nothing to install. Use Python 3.11 or newer and run `python -m seqmill` from the repository root. To run the checks, install the development tools first with `python -m pip install pytest ruff`.

## Usage

Run the commands below as `python -m seqmill COMMAND ...` from the repository root. File arguments are paths or `-` for standard input.

```text
usage: seqmill stats FILE...
       seqmill revcomp [--width N] FILE...
       seqmill translate [--frame F] [--to-stop] [--width N] FILE...
       seqmill orfs [--min N] [--strand both|forward|reverse] FILE...
       seqmill digest --enzyme NAME[,NAME]... [--circular] FILE...
       seqmill align [--mode global|local] [--match N] [--mismatch N] [--gap N] [--width N] A B
       seqmill motif [--mismatches N] [--both-strands] PATTERN FILE...
       seqmill qc [--trim-end Q] [--min-length N] [--min-quality Q] FILE...
```

For example, save this supplied fixture as `reads.fa` in the repository root:

```text
>seq1 first read
ACGTACGTAC
GTNNAC
>seq2
GGGCCC
>seq3 empty
```

Run `python -m seqmill stats reads.fa` from the repository root. Its output is:

```text
file      records  bases  min  max  mean  n50   gc%  n
reads.fa        3     22    0   16   7.3   16  65.0  2
```

## Commands

- `stats FILE...` reports each file's record count, total bases, minimum and maximum sequence lengths, mean length, N50, GC percentage and number of `N` bases.
- `revcomp [--width N] FILE...` writes reverse-complemented records. FASTA and FASTQ formats are preserved, and FASTQ quality strings are reversed along with their sequences.
- `translate [--frame F] [--to-stop] [--width N] FILE...` writes translated proteins as FASTA. Choose any of the six reading frames; negative frames use the reverse complement, and `--to-stop` truncates translation before the first stop codon.
- `orfs [--min N] [--strand both|forward|reverse] FILE...` reports complete open reading frames from all six frames, or from only the three forward or the three reverse frames with `--strand`. `--min` sets the minimum amino-acid length, excluding the stop codon.
- `digest --enzyme NAME[,NAME]... [--circular] FILE...` reports restriction fragments for one or more named enzymes. Use `--circular` to digest sequences as rings.
- `align [--mode global|local] [--match N] [--mismatch N] [--gap N] [--width N] A B` aligns the first record from each input, globally or locally, with configurable scores and output width.
- `motif [--mismatches N] [--both-strands] PATTERN FILE...` searches for IUPAC DNA motifs, optionally allowing mismatches and checking both strands.
- `qc [--trim-end Q] [--min-length N] [--min-quality Q] FILE...` trims and filters FASTQ reads by end quality, length and mean quality, then reports read and trimming counts on standard error.

## Exit codes

seqmill exits with 0 on success, 1 for a negative answer (a table command with no rows or a local alignment with no alignment), and 2 for an error. Errors are written to standard error and leave standard output empty.

## Development

From the repository root, run the checks with:

```sh
python -m ruff check .
python -m pytest -q
```

---

Built with BAKED

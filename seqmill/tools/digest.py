from seqmill.core.options import OptionSpec, parse_options
from seqmill.core.seq import SeqError, show
from seqmill.core.text import render_table
from seqmill.formats.files import Io, check_stdin_once, load_records

_ENZYMES: dict[str, tuple[str, int]] = {
    "ApaI": ("GGGCCC", 5),
    "BamHI": ("GGATCC", 1),
    "EcoRI": ("GAATTC", 1),
    "EcoRV": ("GATATC", 3),
    "HindIII": ("AAGCTT", 1),
    "KpnI": ("GGTACC", 5),
    "NcoI": ("CCATGG", 1),
    "NotI": ("GCGGCCGC", 2),
    "PstI": ("CTGCAG", 5),
    "SacI": ("GAGCTC", 5),
    "SalI": ("GTCGAC", 1),
    "SmaI": ("CCCGGG", 3),
    "XbaI": ("TCTAGA", 1),
    "XhoI": ("CTCGAG", 1),
}
_ENZYME_NAMES = {name.casefold(): name for name in _ENZYMES}
_EXPECTED_ENZYMES = (
    "ApaI, BamHI, EcoRI, EcoRV, HindIII, KpnI, NcoI, NotI, "
    "PstI, SacI, SalI, SmaI, XbaI or XhoI"
)
_HEADER = ["name", "fragment", "start", "end", "length"]


def _selected_enzymes(value: str) -> list[str]:
    selected: list[str] = []
    seen: set[str] = set()
    for item in value.split(","):
        name = _ENZYME_NAMES.get(item.casefold())
        if name is None:
            raise SeqError(
                f'unknown enzyme {show(item)} (expected {_EXPECTED_ENZYMES})'
            )
        if name not in seen:
            selected.append(name)
            seen.add(name)
    return selected


def _find_cuts(sequence: str, enzymes: list[str], circular: bool) -> list[int]:
    length = len(sequence)
    cuts: set[int] = set()

    for name in enzymes:
        site, offset = _ENZYMES[name]
        if length < len(site):
            continue

        if circular:
            starts = range(length)
        else:
            starts = range(length - len(site) + 1)

        for position in starts:
            end = position + len(site)
            if end <= length:
                window = sequence[position:end]
            else:
                window = sequence[position:] + sequence[: end - length]
            if window != site:
                continue

            cut = position + offset
            if circular and cut > length:
                cut -= length
            if circular and 1 <= cut <= length or not circular and 1 <= cut < length:
                cuts.add(cut)

    return sorted(cuts)


def _fragment_rows(
    name: str, length: int, cuts: list[int], circular: bool
) -> list[list[str]]:
    if not cuts:
        return [[name, "1", "1", str(length), str(length)]]

    rows: list[list[str]] = []
    if circular and length:
        for index, cut in enumerate(cuts):
            next_cut = cuts[(index + 1) % len(cuts)]
            fragment_length = (next_cut - cut) % length or length
            start = cut % length + 1
            rows.append(
                [
                    name,
                    str(index + 1),
                    str(start),
                    str(next_cut),
                    str(fragment_length),
                ]
            )
        return rows

    previous_cut = 0
    for index, cut in enumerate(cuts, start=1):
        rows.append(
            [name, str(index), str(previous_cut + 1), str(cut), str(cut - previous_cut)]
        )
        previous_cut = cut
    rows.append(
        [
            name,
            str(len(cuts) + 1),
            str(previous_cut + 1),
            str(length),
            str(length - previous_cut),
        ]
    )
    return rows


def run(argv: list[str], io: Io) -> int:
    parsed = parse_options(
        argv,
        [
            OptionSpec("--enzyme", "text"),
            OptionSpec("--circular", "flag", default=False),
        ],
    )

    enzyme_value = parsed.values["--enzyme"]
    if enzyme_value is None:
        raise SeqError("missing --enzyme")
    enzymes = _selected_enzymes(enzyme_value)

    paths = parsed.args
    if not paths:
        raise SeqError("missing FILE")
    check_stdin_once(paths)

    circular = parsed.values["--circular"]
    rows: list[list[str]] = []
    for path in paths:
        loaded = load_records(io, path)
        for record in loaded.records:
            cuts = _find_cuts(record.seq, enzymes, circular)
            rows.extend(_fragment_rows(record.name, len(record.seq), cuts, circular))

    if not rows:
        return 1

    lines = render_table(_HEADER, rows, {1, 2, 3, 4})
    io.stdout.write("\n".join(lines) + "\n")
    return 0

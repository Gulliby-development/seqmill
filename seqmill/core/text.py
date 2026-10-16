def _int_to_decimal(value: int) -> str:
    """Convert an integer to decimal using bounded-size string conversions."""
    chunks: list[int] = []
    while value:
        value, chunk = divmod(value, 1_000_000_000)
        chunks.append(chunk)
    if not chunks:
        return "0"
    return str(chunks[-1]) + "".join(f"{chunk:09d}" for chunk in reversed(chunks[:-1]))


def fixed1(num: int, den: int) -> str:
    """Format ``num / den`` to one decimal place, rounding half up."""
    tenths, remainder = divmod(num * 10, den)
    if remainder * 2 >= den:
        tenths += 1
    return f"{_int_to_decimal(tenths // 10)}.{tenths % 10}"


def render_table(header: list[str], rows: list[list[str]], right: set[int]) -> list[str]:
    """Render a character-aligned table without line endings."""
    widths = [len(cell) for cell in header]
    for row in rows:
        for index, cell in enumerate(row):
            widths[index] = max(widths[index], len(cell))

    def render(cells: list[str]) -> str:
        aligned = (
            cells[index].rjust(widths[index])
            if index in right
            else cells[index].ljust(widths[index])
            for index in range(len(header))
        )
        return "  ".join(aligned).rstrip(" ")

    return [render(header), *(render(row) for row in rows)]

import pytest

from seqmill.core.text import fixed1, render_table


@pytest.mark.parametrize(
    ("num", "den", "expected"),
    [
        (0, 1, "0.0"),
        (1, 3, "0.3"),
        (5, 2, "2.5"),
        (1, 20, "0.1"),
        (1, 40, "0.0"),
        (199, 20, "10.0"),
        (200000000000000000001, 20, "10000000000000000000.1"),
    ],
)
def test_fixed1_contract_examples(num: int, den: int, expected: str) -> None:
    assert fixed1(num, den) == expected


@pytest.mark.parametrize(
    ("num", "den", "expected"),
    [
        (1, 40, "0.0"),
        (2, 40, "0.1"),
        (3, 40, "0.1"),
        (6, 40, "0.2"),
        (19, 20, "1.0"),
        (39, 20, "2.0"),
    ],
)
def test_fixed1_rounds_half_up_at_decimal_boundaries(num: int, den: int, expected: str) -> None:
    assert fixed1(num, den) == expected


def test_fixed1_handles_integer_parts_beyond_python_decimal_conversion_limit() -> None:
    assert fixed1(10**4300, 1) == "1" + "0" * 4300 + ".0"


def test_render_table_matches_the_exact_contract_example() -> None:
    assert render_table(["name", "n"], [["a", "2"], ["long", "10"]], {1}) == [
        "name   n",
        "a      2",
        "long  10",
    ]


def test_render_table_returns_a_header_when_there_are_no_rows() -> None:
    assert render_table(["name", "count"], [], {1}) == ["name  count"]


def test_render_table_uses_header_width_when_rows_are_narrower() -> None:
    assert render_table(["identifier", "n"], [["x", "1"]], {1}) == [
        "identifier  n",
        "x           1",
    ]


def test_render_table_right_aligns_headers_as_well_as_cells() -> None:
    assert render_table(["label", "n"], [["x", "123"]], {0, 1}) == [
        "label    n",
        "    x  123",
    ]


def test_render_table_handles_empty_cells_and_removes_trailing_spaces() -> None:
    assert render_table(["id", "value"], [["", "1"], ["long", ""]], {1}) == [
        "id    value",
        "          1",
        "long",
    ]


def test_render_table_measures_unicode_by_characters() -> None:
    assert render_table(["x", "n"], [["ééé", "10"], ["猫猫", "2"]], {1}) == [
        "x     n",
        "ééé  10",
        "猫猫    2",
    ]

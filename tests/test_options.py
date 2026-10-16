from dataclasses import FrozenInstanceError

import pytest

from seqmill.core.options import OptionSpec, Parsed, parse_options
from seqmill.core.seq import SeqError


def test_option_spec_fields_and_defaults_match_the_contract() -> None:
    spec = OptionSpec("--mode", "choice")
    assert (spec.name, spec.kind, spec.default, spec.lo, spec.hi, spec.choices) == (
        "--mode",
        "choice",
        None,
        0,
        0,
        (),
    )


def test_option_spec_is_frozen() -> None:
    spec = OptionSpec("--flag", "flag")
    with pytest.raises(FrozenInstanceError):
        spec.name = "--changed"


def test_parsed_fields_and_mutability_match_the_contract() -> None:
    parsed = Parsed({"--name": "value"}, ["-"])
    parsed.values["--name"] = "changed"
    parsed.args.append("next")
    assert parsed.values == {"--name": "changed"}
    assert parsed.args == ["-", "next"]


def test_absent_options_keep_arbitrary_defaults_and_exact_keys() -> None:
    sentinel = object()
    specs = [
        OptionSpec("--flag", "flag", default=False),
        OptionSpec("--count", "int", default=sentinel, lo=0, hi=9),
        OptionSpec("--mode", "choice", default=None, choices=("a", "b")),
        OptionSpec("--text", "text", default="initial"),
    ]
    parsed = parse_options(["input.fa"], specs)
    assert parsed.values == {
        "--flag": False,
        "--count": sentinel,
        "--mode": None,
        "--text": "initial",
    }
    assert list(parsed.values) == ["--flag", "--count", "--mode", "--text"]
    assert parsed.args == ["input.fa"]


def test_supplied_flag_is_true_without_consuming_following_arguments() -> None:
    parsed = parse_options(["--verbose", "input.fa"], [OptionSpec("--verbose", "flag")])
    assert parsed.values == {"--verbose": True}
    assert parsed.args == ["input.fa"]


def test_positional_arguments_keep_order_and_hyphen_is_positional() -> None:
    parsed = parse_options(["first", "-", "last"], [])
    assert parsed.values == {}
    assert parsed.args == ["first", "-", "last"]


def test_equals_value_splits_only_at_the_first_equals_sign() -> None:
    parsed = parse_options(
        ["--text=left=middle=right"],
        [OptionSpec("--text", "text")],
    )
    assert parsed.values == {"--text": "left=middle=right"}


def test_equals_form_accepts_an_empty_text_value() -> None:
    parsed = parse_options(["--text="], [OptionSpec("--text", "text", default="fallback")])
    assert parsed.values == {"--text": ""}


def test_separate_form_consumes_help_token_as_text_value() -> None:
    parsed = parse_options(["--text", "--help"], [OptionSpec("--text", "text")])
    assert parsed.values == {"--text": "--help"}
    assert parsed.args == []


def test_separate_form_consumes_a_leading_dash_value() -> None:
    parsed = parse_options(["--text", "-leading"], [OptionSpec("--text", "text")])
    assert parsed.values == {"--text": "-leading"}


def test_text_values_preserve_whitespace_unicode_and_punctuation() -> None:
    value = "  café	猫!  "
    parsed = parse_options(
        [f"--text={value}"],
        [OptionSpec("--text", "text")],
    )
    assert parsed.values["--text"] == value


def test_last_supplied_option_value_wins() -> None:
    parsed = parse_options(
        ["--count=2", "--count", "7"],
        [OptionSpec("--count", "int", lo=0, hi=9)],
    )
    assert parsed.values == {"--count": 7}


def test_invalid_earlier_value_is_ignored_when_overwritten_by_valid_value() -> None:
    parsed = parse_options(
        ["--count=not-an-int", "--count=3"],
        [OptionSpec("--count", "int", lo=0, hi=9)],
    )
    assert parsed.values == {"--count": 3}


def test_invalid_final_repeated_value_is_rejected() -> None:
    with pytest.raises(SeqError) as exc_info:
        parse_options(
            ["--count=3", "--count=not-an-int"],
            [OptionSpec("--count", "int", lo=0, hi=9)],
        )
    assert exc_info.value.message == (
        'invalid --count "not-an-int" (expected an integer from 0 to 9)'
    )


@pytest.mark.parametrize(
    ("argv", "expected"),
    [
        (["-x"], 'unknown option "-x"'),
        (["--"], 'unknown option "--"'),
        (["--help"], 'unknown option "--help"'),
        (["--bad\noption"], 'unknown option "--bad?option"'),
        (["--badé"], 'unknown option "--bad?"'),
    ],
)
def test_unknown_option_tokens_have_exact_sanitized_diagnostics(
    argv: list[str], expected: str
) -> None:
    with pytest.raises(SeqError) as exc_info:
        parse_options(argv, [])
    assert exc_info.value.message == expected


def test_flag_written_with_equals_is_an_unknown_option() -> None:
    with pytest.raises(SeqError) as exc_info:
        parse_options(["--flag=yes"], [OptionSpec("--flag", "flag")])
    assert exc_info.value.message == 'unknown option "--flag=yes"'


def test_last_non_flag_without_a_value_has_exact_diagnostic() -> None:
    with pytest.raises(SeqError) as exc_info:
        parse_options(["--width"], [OptionSpec("--width", "int", lo=0, hi=1000)])
    assert exc_info.value.message == "--width needs a value"


def test_first_syntax_error_wins_over_a_later_syntax_error() -> None:
    with pytest.raises(SeqError) as exc_info:
        parse_options(
            ["--unknown", "--flag=yes"],
            [OptionSpec("--flag", "flag")],
        )
    assert exc_info.value.message == 'unknown option "--unknown"'


def test_leftmost_flag_equals_error_wins_over_a_later_unknown_option() -> None:
    with pytest.raises(SeqError) as exc_info:
        parse_options(
            ["--flag=yes", "--unknown"],
            [OptionSpec("--flag", "flag")],
        )
    assert exc_info.value.message == 'unknown option "--flag=yes"'


def test_syntax_errors_are_found_before_deferred_value_errors() -> None:
    with pytest.raises(SeqError) as exc_info:
        parse_options(
            ["--count=bad", "--unknown"],
            [OptionSpec("--count", "int", lo=0, hi=9)],
        )
    assert exc_info.value.message == 'unknown option "--unknown"'


@pytest.mark.parametrize(
    ("raw_value", "expected"),
    [
        ("0", 0),
        ("3", 3),
        ("-3", -3),
        ("999999999", 999999999),
        ("-999999999", -999999999),
    ],
)
def test_integer_parser_accepts_contract_lexemes_and_endpoints(
    raw_value: str, expected: int
) -> None:
    parsed = parse_options(
        [f"--count={raw_value}"],
        [OptionSpec("--count", "int", lo=-999999999, hi=999999999)],
    )
    assert parsed.values["--count"] == expected


@pytest.mark.parametrize(
    ("raw_value", "expected"),
    [
        ("2", 2),
        ("4", 4),
    ],
)
def test_integer_bounds_are_inclusive(raw_value: str, expected: int) -> None:
    parsed = parse_options(
        [f"--count={raw_value}"],
        [OptionSpec("--count", "int", lo=2, hi=4)],
    )
    assert parsed.values["--count"] == expected


@pytest.mark.parametrize(
    ("raw_value", "shown_value"),
    [
        ("-0", "-0"),
        ("+3", "+3"),
        ("007", "007"),
        ("", ""),
        (" ", "?"),
        ("٣", "?"),
        ("3.0", "3.0"),
        ("3e0", "3e0"),
        ("1234567890", "1234567890"),
        ("-4", "-4"),
        ("5", "5"),
    ],
)
def test_integer_parser_rejects_invalid_lexemes_and_out_of_range_values(
    raw_value: str, shown_value: str
) -> None:
    with pytest.raises(SeqError) as exc_info:
        parse_options(
            [f"--count={raw_value}"],
            [OptionSpec("--count", "int", lo=-3, hi=4)],
        )
    assert exc_info.value.message == (
        f'invalid --count "{shown_value}" (expected an integer from -3 to 4)'
    )


def test_rejected_integer_value_is_sanitized_in_its_diagnostic() -> None:
    with pytest.raises(SeqError) as exc_info:
        parse_options(
            ["--count=1\n2"],
            [OptionSpec("--count", "int", lo=0, hi=9)],
        )
    assert exc_info.value.message == 'invalid --count "1?2" (expected an integer from 0 to 9)'


@pytest.mark.parametrize(
    ("choices", "raw_value", "expected"),
    [
        (("only",), "other", "only"),
        (("red", "blue"), "green", "red or blue"),
        (("a", "b", "c"), "x", "a, b or c"),
        (("a", "b", "c", "d", "e", "f"), "x", "a, b, c, d, e or f"),
    ],
)
def test_choice_errors_use_exact_list_punctuation(
    choices: tuple[str, ...], raw_value: str, expected: str
) -> None:
    with pytest.raises(SeqError) as exc_info:
        parse_options(
            [f"--mode={raw_value}"],
            [OptionSpec("--mode", "choice", choices=choices)],
        )
    assert exc_info.value.message == f"invalid --mode \"{raw_value}\" (expected {expected})"


def test_choice_matching_is_case_sensitive() -> None:
    with pytest.raises(SeqError) as exc_info:
        parse_options(
            ["--mode=FAST"],
            [OptionSpec("--mode", "choice", choices=("fast", "slow"))],
        )
    assert exc_info.value.message == 'invalid --mode "FAST" (expected fast or slow)'


def test_frame_choice_diagnostic_matches_the_contract_exactly() -> None:
    with pytest.raises(SeqError) as exc_info:
        parse_options(
            ["--frame=4"],
            [OptionSpec("--frame", "choice", choices=("1", "2", "3", "-1", "-2", "-3"))],
        )
    assert exc_info.value.message == (
        'invalid --frame "4" (expected 1, 2, 3, -1, -2 or -3)'
    )


def test_validation_errors_follow_spec_order_not_argv_order() -> None:
    specs = [
        OptionSpec("--first", "int", lo=0, hi=9),
        OptionSpec("--second", "int", lo=0, hi=9),
    ]
    with pytest.raises(SeqError) as exc_info:
        parse_options(["--second=bad", "--first=also-bad"], specs)
    assert exc_info.value.message == (
        'invalid --first "also-bad" (expected an integer from 0 to 9)'
    )

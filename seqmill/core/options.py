from dataclasses import dataclass

from seqmill.core.seq import SeqError, show


@dataclass(frozen=True)
class OptionSpec:
    name: str
    kind: str
    default: object = None
    lo: int = 0
    hi: int = 0
    choices: tuple[str, ...] = ()


@dataclass
class Parsed:
    values: dict
    args: list[str]


def _is_integer_text(value: str) -> bool:
    digits = value.removeprefix("-")
    return (
        value == "0"
        or (
            1 <= len(digits) <= 9
            and digits[0] in "123456789"
            and all(char in "0123456789" for char in digits)
        )
    )


def _choice_list(choices: tuple[str, ...]) -> str:
    if len(choices) < 2:
        return choices[0] if choices else ""
    if len(choices) == 2:
        return f"{choices[0]} or {choices[1]}"
    return f"{', '.join(choices[:-1])} or {choices[-1]}"


def parse_options(argv: list[str], specs: list[OptionSpec]) -> Parsed:
    by_name = {spec.name: spec for spec in specs}
    values: dict[str, object] = {spec.name: spec.default for spec in specs}
    supplied: dict[str, str] = {}
    args: list[str] = []

    index = 0
    while index < len(argv):
        token = argv[index]
        if not token.startswith("-") or token == "-":
            args.append(token)
            index += 1
            continue

        has_value = "=" in token
        if has_value:
            name, raw_value = token.split("=", 1)
        else:
            name, raw_value = token, ""

        spec = by_name.get(name)
        if spec is None or (spec.kind == "flag" and has_value):
            raise SeqError(f"unknown option {show(token)}")

        if spec.kind == "flag":
            values[name] = True
        else:
            if not has_value:
                if index + 1 >= len(argv):
                    raise SeqError(f"{name} needs a value")
                index += 1
                raw_value = argv[index]
            supplied[name] = raw_value
            values[name] = raw_value

        index += 1

    for spec in specs:
        if spec.name not in supplied:
            continue

        raw_value = supplied[spec.name]
        if spec.kind == "int":
            if not _is_integer_text(raw_value):
                raise SeqError(
                    f'invalid {spec.name} {show(raw_value)} '
                    f"(expected an integer from {spec.lo} to {spec.hi})"
                )
            parsed_value = int(raw_value)
            if not spec.lo <= parsed_value <= spec.hi:
                raise SeqError(
                    f'invalid {spec.name} {show(raw_value)} '
                    f"(expected an integer from {spec.lo} to {spec.hi})"
                )
            values[spec.name] = parsed_value
        elif spec.kind == "choice" and raw_value not in spec.choices:
            expected = _choice_list(spec.choices)
            raise SeqError(f"invalid {spec.name} {show(raw_value)} (expected {expected})")

    return Parsed(values=values, args=args)

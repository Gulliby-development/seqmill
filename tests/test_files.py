import builtins
from dataclasses import fields, is_dataclass
from inspect import signature
from io import BytesIO, StringIO
from pathlib import Path
from typing import BinaryIO, Self, TextIO

import pytest

from seqmill.core.seq import FormatError, Record, SeqError
from seqmill.formats import files
from seqmill.formats.files import Io, Loaded


def make_io(data: bytes = b"") -> Io:
    return Io(BytesIO(data), StringIO(), StringIO())


def test_files_api_reuses_core_types_and_has_exact_container_fields() -> None:
    assert files.Record is Record
    assert files.FormatError is FormatError
    assert files.SeqError is SeqError

    assert is_dataclass(Io)
    assert [field.name for field in fields(Io)] == ["stdin", "stdout", "stderr"]
    assert [field.type for field in fields(Io)] == [BinaryIO, TextIO, TextIO]
    assert is_dataclass(Loaded)
    assert [field.name for field in fields(Loaded)] == ["display", "kind", "records"]
    assert [field.type for field in fields(Loaded)] == [str, str, list[Record]]

    assert list(signature(files.display_path).parameters) == ["path"]
    assert signature(files.display_path).return_annotation is str
    assert list(signature(files.split_lines).parameters) == ["text"]
    assert signature(files.split_lines).return_annotation == list[str]
    read_text_signature = signature(files.read_text)
    assert list(read_text_signature.parameters) == ["io", "path"]
    assert read_text_signature.parameters["io"].annotation is Io
    assert read_text_signature.parameters["path"].annotation is str
    assert read_text_signature.return_annotation == tuple[str, str]
    check_stdin_signature = signature(files.check_stdin_once)
    assert list(check_stdin_signature.parameters) == ["paths"]
    assert check_stdin_signature.parameters["paths"].annotation == list[str]
    assert check_stdin_signature.return_annotation is None
    load_records_signature = signature(files.load_records)
    assert list(load_records_signature.parameters) == ["io", "path"]
    assert load_records_signature.parameters["io"].annotation is Io
    assert load_records_signature.parameters["path"].annotation is str
    assert load_records_signature.return_annotation is Loaded


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("", []),
        ("\ufeff", []),
        ("\n", [""]),
        ("a\n", ["a"]),
        ("a\n\n", ["a", ""]),
        ("a", ["a"]),
        ("a\r\n", ["a"]),
        ("a\r", ["a\r"]),
        ("a\r\r\n", ["a\r"]),
        ("\ufeff>a\r\nAC", [">a", "AC"]),
        ("\ufeff\ufeffa", ["\ufeffa"]),
        ("a\rb\u2028c\u0085d", ["a\rb\u2028c\u0085d"]),
        ("a\n\ufeffb", ["a", "\ufeffb"]),
    ],
)
def test_split_lines_preserves_only_the_specified_lf_and_bom_rules(
    text: str, expected: list[str]
) -> None:
    assert files.split_lines(text) == expected


@pytest.mark.parametrize(
    ("path", "expected"),
    [
        ("-", "<stdin>"),
        ("C:\\data\\reads.fa", "C:/data/reads.fa"),
        ("./a/../b", "./a/../b"),
        ("folder/file.fa", "folder/file.fa"),
    ],
)
def test_display_path_only_maps_stdin_and_slashes(path: str, expected: str) -> None:
    assert files.display_path(path) == expected


@pytest.mark.parametrize(
    ("paths", "message"),
    [
        ([], None),
        (["a"], None),
        (["-"], None),
        (["a", "a"], None),
        (["-", "a"], None),
        (["-", "-"], "standard input given twice"),
        (["a", "-", "b", "-", "-"], "standard input given twice"),
    ],
)
def test_check_stdin_once_counts_only_duplicate_stdin_paths(
    paths: list[str], message: str | None
) -> None:
    if message is None:
        assert files.check_stdin_once(paths) is None
    else:
        with pytest.raises(SeqError, match="^standard input given twice$"):
            files.check_stdin_once(paths)


@pytest.mark.parametrize(
    ("data", "expected"),
    [
        (b"", ("<stdin>", "")),
        (b"\xef\xbb\xbf>a\r\nAC", ("<stdin>", "\ufeff>a\r\nAC")),
        (b">\xc3\xa9\r\nAC", ("<stdin>", ">\u00e9\r\nAC")),
    ],
)
def test_read_text_decodes_utf8_without_removing_bom_or_newlines(
    data: bytes, expected: tuple[str, str]
) -> None:
    io = make_io(data)

    assert files.read_text(io, "-") == expected


@pytest.mark.parametrize(
    ("source", "data"),
    [
        ("stdin", b">s\nAC\xff\n"),
        ("stdin", b"\x80"),
        ("stdin", b"\xc3("),
        ("file", b">s\nAC\xff\n"),
        ("file", b"\xff"),
        ("file", b"\xc3("),
    ],
)
def test_read_text_rejects_malformed_utf8_from_binary_channels(
    tmp_path: Path, source: str, data: bytes
) -> None:
    io = make_io(data if source == "stdin" else b"")
    path = "-"

    if source == "file":
        path = str(tmp_path / "invalid.fa")
        Path(path).write_bytes(data)

    with pytest.raises(SeqError) as caught:
        files.read_text(io, path)

    expected_display = "<stdin>" if path == "-" else path.replace("\\", "/")
    assert caught.value.message == f"{expected_display}: input is not valid UTF-8"


def test_read_text_maps_missing_file_to_display_diagnostic(tmp_path: Path) -> None:
    path = str(tmp_path / "missing.fa")

    with pytest.raises(SeqError) as caught:
        files.read_text(make_io(), path)

    expected_display = path.replace("\\", "/")
    assert caught.value.message == f"{expected_display}: cannot read the file"


def test_read_text_maps_directory_to_display_diagnostic(tmp_path: Path) -> None:
    path = str(tmp_path)

    with pytest.raises(SeqError) as caught:
        files.read_text(make_io(), path)

    expected_display = path.replace("\\", "/")
    assert caught.value.message == f"{expected_display}: cannot read the file"


def test_read_text_maps_embedded_nul_path_to_display_diagnostic(tmp_path: Path) -> None:
    path = str(tmp_path / "bad\x00.fa")

    with pytest.raises(SeqError) as caught:
        files.read_text(make_io(), path)

    expected_display = path.replace("\\", "/")
    assert caught.value.message == f"{expected_display}: cannot read the file"


def test_read_text_maps_an_oserror_during_binary_read(monkeypatch: pytest.MonkeyPatch) -> None:
    class FailingReader:
        def __enter__(self) -> Self:
            return self

        def __exit__(self, exc_type: object, exc_value: object, traceback: object) -> bool:
            return False

        def read(self) -> bytes:
            raise OSError("read failed")

    def open_failing_reader(path: str, mode: str) -> FailingReader:
        assert path == "unreadable.fa"
        assert mode == "rb"
        return FailingReader()

    monkeypatch.setattr(builtins, "open", open_failing_reader)

    with pytest.raises(SeqError) as caught:
        files.read_text(make_io(), "unreadable.fa")

    assert caught.value.message == "unreadable.fa: cannot read the file"


def test_read_text_opens_the_original_path_in_binary_mode(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    path = tmp_path / "reads.fa"
    path.write_bytes(b">a\r\nAC")
    original_open = builtins.open
    calls: list[tuple[object, str]] = []

    def recording_open(file: object, mode: str = "r", *args: object, **kwargs: object):
        calls.append((file, mode))
        return original_open(file, mode, *args, **kwargs)

    monkeypatch.setattr(builtins, "open", recording_open)

    display, text = files.read_text(make_io(), str(path))

    assert display == str(path).replace("\\", "/")
    assert text == ">a\r\nAC"
    assert len(calls) == 1
    opened_path, opened_mode = calls[0]
    assert Path(str(opened_path)).resolve() == path.resolve()
    assert opened_mode == "rb"


@pytest.mark.parametrize(
    ("data", "expected_kind", "expected_records"),
    [
        (b"", "fasta", []),
        (
            b"\xef\xbb\xbf>a\r\nac",
            "fasta",
            [Record("a", "", "AC")],
        ),
        (
            b"@a desc\nac\n+ignored\n!~\n",
            "fastq",
            [Record("a", "desc", "AC", "!~")],
        ),
    ],
)
def test_load_records_detects_format_and_preserves_loaded_metadata(
    data: bytes, expected_kind: str, expected_records: list[Record]
) -> None:
    loaded = files.load_records(make_io(data), "-")

    assert loaded == Loaded("<stdin>", expected_kind, expected_records)


def test_load_records_does_not_skip_a_leading_blank_line_for_detection() -> None:
    with pytest.raises(SeqError) as caught:
        files.load_records(make_io(b"\n@a\nA\n+\n!"), "-")

    assert caught.value.message == "<stdin>:2: sequence data before the first header"


def test_load_records_wraps_fasta_error_with_display_and_absolute_line() -> None:
    with pytest.raises(SeqError) as caught:
        files.load_records(make_io(b">a\nACGR\n"), "-")

    assert caught.value.message == '<stdin>:2: invalid character "R" in sequence'


def test_load_records_wraps_fastq_error_with_display_and_absolute_line() -> None:
    with pytest.raises(SeqError) as caught:
        files.load_records(make_io(b"@a\nA\nx\n!"), "-")

    assert caught.value.message == "<stdin>:3: expected + line"


def test_load_records_decodes_before_reporting_apparent_format_errors() -> None:
    with pytest.raises(SeqError) as caught:
        files.load_records(make_io(b"not a header\n\xff"), "-")

    assert caught.value.message == "<stdin>: input is not valid UTF-8"


def test_helpers_leave_caller_streams_open_and_outputs_untouched() -> None:
    io = make_io(b">a\nAC\n")

    assert files.read_text(io, "-") == ("<stdin>", ">a\nAC\n")
    io.stdin.seek(0)
    assert files.load_records(io, "-") == Loaded("<stdin>", "fasta", [Record("a", "", "AC")])
    assert io.stdin.closed is False
    assert io.stdout.getvalue() == ""
    assert io.stderr.getvalue() == ""

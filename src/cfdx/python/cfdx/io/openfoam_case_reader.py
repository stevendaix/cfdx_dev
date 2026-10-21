"""Native, dependency-free reader for common OpenFOAM ASCII case files.

The reader intentionally supports a documented subset of OpenFOAM syntax. It
never guesses binary payloads as text: unsupported encodings raise a clear
error and remain available in the original case directory.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class OpenFOAMFieldSummary:
    """Metadata for a field available at one OpenFOAM time."""

    name: str
    dimensions: tuple[int, ...] | None
    boundary_types: dict[str, str]


@dataclass(frozen=True)
class OpenFOAMTimeSummary:
    """Fields available in a time directory."""

    name: str
    time: float
    fields: tuple[OpenFOAMFieldSummary, ...]


@dataclass(frozen=True)
class OpenFOAMCaseSummary:
    """Discoverable configuration and result metadata for a case."""

    path: Path
    application: str | None
    times: tuple[OpenFOAMTimeSummary, ...]
    dictionaries: tuple[str, ...]


_TOKEN_RE = re.compile(
    r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|[{}();\[\]]|[^\s{}();\[\]]+'
)
_COMMENT_RE = re.compile(
    r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|//[^\n]*|/\*.*?\*/', re.DOTALL
)
_NUMERIC_RE = re.compile(r"^[+-]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[eE][+-]?\d+)?$")


def _strip_comments(text: str) -> str:
    """Remove comments without altering comment-like text inside quotes."""
    def replace(match: re.Match[str]) -> str:
        token = match.group(0)
        if token.startswith(('"', "'")):
            return token
        return " " * len(token)
    return _COMMENT_RE.sub(replace, text)


def _tokens(text: str) -> list[str]:
    """Tokenize the ASCII subset used by OpenFOAM dictionaries and fields."""
    cleaned = _strip_comments(text)
    if re.search(r"\bformat\s+binary\s*;", cleaned):
        raise ValueError(
            "Binary OpenFOAM files are not supported by the native Python reader"
        )
    if re.search(r"(?m)^\s*#", cleaned):
        raise ValueError(
            "OpenFOAM preprocessor directives are not supported by the native Python reader"
        )
    return _TOKEN_RE.findall(cleaned)


def _atom(token: str) -> Any:
    if len(token) >= 2 and token[0] == token[-1] and token[0] in {'"', "'"}:
        return token[1:-1]
    if _NUMERIC_RE.match(token):
        value = float(token)
        if (
            math.isfinite(value)
            and value.is_integer()
            and not any(c in token.lower() for c in (".", "e"))
        ):
            return int(value)
        return value
    if token in {"true", "false"}:
        return token == "true"
    return token


def _parse_group(tokens: list[str], index: int) -> tuple[list[Any], int]:
    """Parse a parenthesized list and return the next token index."""
    if index >= len(tokens) or tokens[index] not in {"(", "["}:
        raise ValueError("Expected a list delimiter while parsing OpenFOAM list")
    closing = ")" if tokens[index] == "(" else "]"
    values: list[Any] = []
    index += 1
    while index < len(tokens) and tokens[index] != closing:
        if tokens[index] in {"(", "["}:
            nested, index = _parse_group(tokens, index)
            values.append(nested)
        elif tokens[index] in {"{", "}", ";"}:
            raise ValueError(f"Unexpected token {tokens[index]!r} in OpenFOAM list")
        else:
            values.append(_atom(tokens[index]))
            index += 1
    if index >= len(tokens):
        raise ValueError("Unterminated OpenFOAM list")
    return values, index + 1


def _value(tokens: list[str]) -> Any:
    if not tokens:
        return None
    if tokens[0] in {"(", "["}:
        parsed, next_index = _parse_group(tokens, 0)
        if next_index == len(tokens):
            return parsed
    if len(tokens) == 1:
        return _atom(tokens[0])
    if tokens[0] == "uniform":
        if len(tokens) > 1 and tokens[1] == "(":
            vector, next_index = _parse_group(tokens, 1)
            if next_index == len(tokens):
                return {"uniform": vector}
        if len(tokens) == 2:
            return {"uniform": _atom(tokens[1])}
    return " ".join(str(_atom(token)) for token in tokens)


def _parse_dictionary_tokens(
    tokens: list[str], index: int = 0, *, nested: bool = False
) -> tuple[dict[str, Any], int]:
    result: dict[str, Any] = {}
    while index < len(tokens):
        if tokens[index] == "}":
            if not nested:
                raise ValueError("Unexpected '}' in OpenFOAM dictionary")
            return result, index + 1
        if tokens[index] in {"{", ")", ";"}:
            raise ValueError(
                f"Unexpected token {tokens[index]!r} in OpenFOAM dictionary"
            )
        key = str(_atom(tokens[index]))
        index += 1
        if index < len(tokens) and tokens[index] == "{":
            child, index = _parse_dictionary_tokens(tokens, index + 1, nested=True)
            result[key] = child
            if index < len(tokens) and tokens[index] == ";":
                index += 1
            continue

        value_tokens: list[str] = []
        depth = 0
        while index < len(tokens):
            token = tokens[index]
            if token in {"(", "["}:
                depth += 1
            elif token in {")", "]"}:
                depth -= 1
                if depth < 0:
                    raise ValueError(
                        f"Unmatched closing list delimiter in value for {key!r}"
                    )
            if token == ";" and depth == 0:
                index += 1
                break
            if token == "}" and depth == 0:
                # Some OpenFOAM entries omit the trailing semicolon.
                break
            value_tokens.append(token)
            index += 1
        result[key] = _value(value_tokens)
    if nested:
        raise ValueError("Unterminated OpenFOAM dictionary block")
    return result, index


def _read_dictionary_file(path: Path) -> dict[str, Any]:
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError(f"Not an ASCII OpenFOAM dictionary: {path}") from exc
    tokens = _tokens(text)
    # The FoamFile header is a regular dictionary; remove it before parsing the
    # user entries so callers receive application settings, not file metadata.
    if tokens[:1] == ["FoamFile"]:
        # Parse the body following the balanced FoamFile block.
        brace = tokens.index("{", 1)
        depth = 1
        index = brace + 1
        while index < len(tokens) and depth:
            if tokens[index] == "{":
                depth += 1
            elif tokens[index] == "}":
                depth -= 1
            index += 1
        tokens = tokens[index:]
    parsed, _ = _parse_dictionary_tokens(tokens)
    return parsed


def _time_directories(path: Path) -> list[tuple[float, Path]]:
    times: list[tuple[float, Path]] = []
    for child in path.iterdir():
        if not child.is_dir():
            continue
        try:
            value = float(child.name)
        except ValueError:
            continue
        if math.isfinite(value):
            times.append((value, child))
    return sorted(times, key=lambda item: item[0])


def _read_field_file(path: Path, *, cell_count: int | None = None) -> dict[str, Any]:
    with path.open("rb") as stream:
        prefix = stream.read(4096).decode("ascii", errors="ignore")
    if re.search(r"\bformat\s+binary\s*;", _strip_comments(prefix)):
        raise ValueError(
            "Binary OpenFOAM files are not supported by the native Python reader"
        )
    text = path.read_text(encoding="utf-8")
    tokens = _tokens(text)
    if not tokens or tokens[0] != "FoamFile":
        raise ValueError(f"Missing FoamFile header in {path}")
    header, _ = _parse_dictionary_tokens(tokens, 0)
    foam_header = header.get("FoamFile", {})
    source_class = (
        str(foam_header.get("class", "")) if isinstance(foam_header, dict) else ""
    )
    if not source_class:
        raise ValueError(f"Missing FoamFile class in {path}")

    # The parser above consumes the whole file as a dictionary. For fields,
    # parse the known top-level entries directly, preserving nested boundary
    # dictionaries while extracting dimensions/internalField.
    body_start = 0
    brace = tokens.index("{", 1)
    depth = 1
    body_start = brace + 1
    while body_start < len(tokens) and depth:
        if tokens[body_start] == "{":
            depth += 1
        elif tokens[body_start] == "}":
            depth -= 1
        body_start += 1
    body = tokens[body_start:]
    entries, _ = _parse_dictionary_tokens(body)
    dimensions_raw = entries.get("dimensions")
    dimensions = (
        tuple(int(v) for v in dimensions_raw)
        if isinstance(dimensions_raw, list)
        else ()
    )

    internal_tokens: list[str] | None = None
    for i, token in enumerate(body):
        if token == "internalField":
            end = i + 1
            depth = 0
            while end < len(body):
                if body[end] == "(":
                    depth += 1
                elif body[end] == ")":
                    depth -= 1
                elif body[end] == ";" and depth == 0:
                    break
                end += 1
            internal_tokens = body[i + 1 : end]
            break
    if internal_tokens is None:
        raise ValueError(f"Missing internalField in {path}")

    uniform = bool(internal_tokens and internal_tokens[0] == "uniform")
    if uniform:
        values_tokens = internal_tokens[1:]
        if values_tokens and values_tokens[0] == "(":
            parsed_values, _ = _parse_group(values_tokens, 0)
            raw = np.asarray(parsed_values, dtype=np.float64)
        elif values_tokens:
            raw = np.asarray(float(_atom(values_tokens[0])), dtype=np.float64)
        else:
            raise ValueError(f"Malformed uniform internalField in {path}")
        if cell_count is not None:
            raw = np.broadcast_to(raw, (cell_count,) + raw.shape).copy()
    else:
        list_start = next(
            (i for i, token in enumerate(internal_tokens) if token == "("), None
        )
        if list_start is None:
            raise ValueError(f"Unsupported internalField encoding in {path}")
        parsed_values, _ = _parse_group(internal_tokens, list_start)
        raw = np.asarray(parsed_values, dtype=np.float64)

    boundary: dict[str, Any] = {}
    boundary_fields = entries.get("boundaryField", {})
    if isinstance(boundary_fields, dict):
        for patch, condition in boundary_fields.items():
            if isinstance(condition, dict):
                boundary[str(patch)] = condition

    return {
        "name": path.name,
        "class": source_class,
        "time": 0.0,
        "dimensions": dimensions,
        "internal": raw,
        "boundary": boundary,
        "uniform": uniform,
    }


def inspect_openfoam_case(case_path: str | Path) -> OpenFOAMCaseSummary:
    """Inspect dictionaries, numeric time directories, and field metadata."""
    path = Path(case_path).expanduser().resolve()
    if not path.is_dir():
        raise NotADirectoryError(f"OpenFOAM case directory does not exist: {path}")

    control_path = path / "system" / "controlDict"
    application = None
    if control_path.is_file():
        application = (
            str(_read_dictionary_file(control_path).get("application", "")) or None
        )

    dictionary_paths: list[str] = []
    for folder in ("system", "constant"):
        root = path / folder
        if root.is_dir():
            dictionary_paths.extend(
                str(item.relative_to(path))
                for item in sorted(root.rglob("*"))
                if item.is_file() and not item.name.startswith(".")
            )

    time_summaries = []
    for value, directory in _time_directories(path):
        summaries = []
        for file_path in sorted(directory.iterdir()):
            if not file_path.is_file() or file_path.name.startswith("."):
                continue
            try:
                data = _read_field_file(file_path)
            except (OSError, ValueError, TypeError, IndexError, UnicodeError):
                continue
            summaries.append(
                OpenFOAMFieldSummary(
                    name=file_path.name,
                    dimensions=data["dimensions"] or None,
                    boundary_types={
                        patch: str(condition.get("type"))
                        for patch, condition in data["boundary"].items()
                        if condition.get("type") is not None
                    },
                )
            )
        time_summaries.append(
            OpenFOAMTimeSummary(
                name=directory.name,
                time=value,
                fields=tuple(summaries),
            )
        )
    return OpenFOAMCaseSummary(
        path=path,
        application=application,
        times=tuple(time_summaries),
        dictionaries=tuple(dictionary_paths),
    )


def read_openfoam_field(
    case_path: str | Path, time: str | float, field_name: str
) -> dict[str, Any]:
    """Read a single field's internal and boundary values from an ASCII case."""
    path = Path(case_path).expanduser().resolve()
    if not path.is_dir():
        raise NotADirectoryError(f"OpenFOAM case directory does not exist: {path}")
    if time == "latest":
        directories = _time_directories(path)
        if not directories:
            raise ValueError(f"No numeric OpenFOAM time directories found in {path}")
        time_value, directory = directories[-1]
    else:
        time_value = float(time)
        directory = path / str(time)
        if not directory.is_dir():
            # Normalize values such as 1.0 to a directory named "1".
            matches = [
                item for value, item in _time_directories(path) if value == time_value
            ]
            if not matches:
                raise FileNotFoundError(
                    f"OpenFOAM time directory {time!r} not found in {path}"
                )
            directory = matches[0]
    field_path = directory / field_name
    data = _read_field_file(field_path)
    data["time"] = time_value
    return data


def read_openfoam_time_fields(
    case_path: str | Path,
    time: str | float = "latest",
    *,
    cell_count: int | None = None,
) -> dict[str, dict[str, Any]]:
    """Read all ASCII field files from a selected time directory.

    Uniform scalar/vector fields are expanded to cell_count when supplied.
    Tensor and non-volume fields are retained as source data for explicit gap
    reporting by the adapter; they are never coerced to scalar/vector fields.
    """
    path = Path(case_path).expanduser().resolve()
    if not path.is_dir():
        raise NotADirectoryError(f"OpenFOAM case directory does not exist: {path}")
    if time == "latest":
        directories = _time_directories(path)
        if not directories:
            raise ValueError(f"No numeric OpenFOAM time directories found in {path}")
        time_value, directory = directories[-1]
    else:
        time_value = float(time)
        matches = [item for value, item in _time_directories(path) if value == time_value]
        if not matches:
            raise FileNotFoundError(f"OpenFOAM time directory {time!r} not found in {path}")
        time_value, directory = matches[0]

    fields: dict[str, dict[str, Any]] = {}
    for field_path in sorted(directory.iterdir()):
        if not field_path.is_file() or field_path.name.startswith("."):
            continue
        try:
            data = _read_field_file(field_path, cell_count=cell_count)
        except (OSError, ValueError, TypeError, IndexError, UnicodeError) as exc:
            raise ValueError(
                f"Cannot decode field {field_path.name} at time {directory.name}: {exc}"
            ) from exc
        data["time"] = time_value
        fields[field_path.name] = data
    return fields


def read_openfoam_dictionaries(case_path: str | Path) -> dict[str, Any]:
    """Parse known OpenFOAM setup dictionaries into JSON-compatible values."""
    path = Path(case_path).expanduser().resolve()
    if not path.is_dir():
        raise NotADirectoryError(f"OpenFOAM case directory does not exist: {path}")
    candidates = (
        "system/controlDict",
        "system/fvSchemes",
        "system/fvSolution",
        "system/fvOptions",
        "system/fvConstraints",
        "system/decomposeParDict",
        "constant/transportProperties",
        "constant/physicalProperties",
        "constant/thermophysicalProperties",
        "constant/turbulenceProperties",
        "constant/momentumTransport",
        "constant/radiationProperties",
        "constant/g",
    )
    dictionaries: dict[str, Any] = {}
    for relative in candidates:
        file_path = path / relative
        if file_path.is_file():
            try:
                dictionaries[relative] = _read_dictionary_file(file_path)
            except (OSError, ValueError, TypeError, UnicodeError) as exc:
                raise ValueError(
                    f"Cannot parse OpenFOAM dictionary {relative}: {exc}"
                ) from exc
    return dictionaries

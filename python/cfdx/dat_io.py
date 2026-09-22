"""Reader for the native CFDX DAT numerical checkpoint format.

The reader is intentionally solver-independent: it discovers the fields present
in a checkpoint and exposes them to the GUI/viewer without assuming that a
particular physics model is enabled.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math


@dataclass(frozen=True)
class DatField:
    name: str
    dimension: int
    values: list[float]

    @property
    def is_vector(self) -> bool:
        return self.dimension > 1

    def component(self, index: int) -> list[float]:
        if index < 0 or index >= self.dimension:
            raise IndexError(index)
        return self.values[index::self.dimension]


@dataclass(frozen=True)
class DatRestart:
    version: int
    cells: int
    iteration: int
    time: float
    fields: dict[str, DatField]


def _read_tokens(path: Path) -> list[str]:
    try:
        text = path.read_text(encoding="utf-8")
    except OSError:
        raise
    tokens = text.split()
    if not tokens:
        raise ValueError("empty DAT restart")
    return tokens


def read_dat_restart(path: str | Path) -> DatRestart:
    """Parse a CFDX DAT checkpoint and discover all cell fields."""
    path = Path(path)
    tokens = _read_tokens(path)
    pos = 0

    def take(label: str) -> str:
        nonlocal pos
        if pos >= len(tokens):
            raise ValueError(f"unexpected end of DAT while reading {label}")
        value = tokens[pos]
        pos += 1
        return value

    if take("DAT magic") != "CFDX-DAT":
        raise ValueError("unsupported DAT format")
    try:
        version = int(take("DAT version"))
    except ValueError as exc:
        raise ValueError("invalid DAT version") from exc
    if version not in (1, 2):
        raise ValueError(f"unsupported DAT version {version}")

    if take("cells key") != "cells":
        raise ValueError("expected 'cells'")
    try:
        cells = int(take("cell count"))
    except ValueError as exc:
        raise ValueError("invalid cell count") from exc
    if cells < 0:
        raise ValueError("cell count must be non-negative")

    if take("iteration key") != "iteration":
        raise ValueError("expected 'iteration'")
    try:
        iteration = int(take("iteration"))
    except ValueError as exc:
        raise ValueError("invalid iteration") from exc
    if iteration < 0:
        raise ValueError("iteration must be non-negative")

    if take("time key") != "time":
        raise ValueError("expected 'time'")
    try:
        time = float(take("time"))
    except ValueError as exc:
        raise ValueError("invalid time") from exc
    if not math.isfinite(time):
        raise ValueError("time must be finite")

    fields: dict[str, DatField] = {}
    while pos < len(tokens):
        if take("field key") != "field":
            raise ValueError("expected 'field'")
        name = take("field name")
        if not name or name in fields:
            raise ValueError(f"duplicate or empty field name: {name!r}")
        try:
            dimension = int(take(f"{name} dimension"))
        except ValueError as exc:
            raise ValueError(f"invalid dimension for field {name!r}") from exc
        if dimension <= 0:
            raise ValueError(f"field {name!r} dimension must be positive")

        count = cells * dimension
        values: list[float] = []
        for _ in range(count):
            try:
                value = float(take(f"{name} values"))
            except ValueError as exc:
                raise ValueError(f"invalid value in field {name!r}") from exc
            if not math.isfinite(value):
                raise ValueError(f"non-finite value in field {name!r}")
            values.append(value)
        fields[name] = DatField(name, dimension, values)

    if not fields:
        raise ValueError("DAT restart contains no fields")
    return DatRestart(version, cells, iteration, time, fields)

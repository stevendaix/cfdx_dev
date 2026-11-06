"""Versioned line protocol for committed live point-probe samples."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass

_MARKER = "CFDX_PROBE_SAMPLE_V1"
_HEX = re.compile(r"[0-9a-fA-F]+\Z")


@dataclass(frozen=True)
class ProbeStreamSample:
    """One committed native probe sample received from solver stdout."""

    name: str
    field: str
    unit: str
    iteration: int
    time: float
    value: float


def parse_probe_stream_line(line: str) -> ProbeStreamSample | None:
    """Parse a ``CFDX_PROBE_SAMPLE_V1`` line; ignore unrelated solver output.

    Probe names are UTF-8 hex-encoded by the native solver so spaces, tabs and
    other text characters cannot alter the tab-delimited transport framing.
    """
    if line != _MARKER and not line.startswith(_MARKER + "\t"):
        return None
    columns = line.split("\t")
    if len(columns) != 7 or columns[0] != _MARKER:
        raise ValueError("malformed CFDX_PROBE_SAMPLE_V1 record: expected 7 columns")
    name_hex, field, unit, iteration_text, time_text, value_text = columns[1:]
    if not _HEX.fullmatch(name_hex) or len(name_hex) % 2:
        raise ValueError("malformed CFDX_PROBE_SAMPLE_V1 probe name encoding")
    try:
        name = bytes.fromhex(name_hex).decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError("probe stream name is not valid UTF-8") from exc
    if not name:
        raise ValueError("probe stream name must not be empty")
    if not field or not unit:
        raise ValueError("probe stream field and unit must not be empty")
    try:
        iteration = int(iteration_text)
        time_value = float(time_text)
        value = float(value_text)
    except ValueError as exc:
        raise ValueError("probe stream iteration/time/value must be numeric") from exc
    if iteration < 0:
        raise ValueError("probe stream iteration must not be negative")
    if not math.isfinite(time_value) or not math.isfinite(value):
        raise ValueError("probe stream time and value must be finite")
    return ProbeStreamSample(name, field, unit, iteration, time_value, value)

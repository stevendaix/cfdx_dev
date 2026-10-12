"""Quantitative point-probe validation contracts for GUI post-processing."""
from __future__ import annotations
from dataclasses import dataclass
from collections.abc import Sequence
from pathlib import Path
import math

# Schema written by the native probe exporter (`cfdx::io::write_probe_csv`) and
# read back here, so a file produced by a solver run is consumable by the
# application layer without a translation step.
PROBE_CSV_MAGIC = "cfdx-probe-csv"
PROBE_CSV_VERSION = 1
_PROBE_CSV_COLUMNS = ("probe", "iteration", "time", "value")


@dataclass(frozen=True)
class ProbeSample:
    iteration: int
    time: float
    value: float

@dataclass(frozen=True)
class ProbeSeries:
    name: str
    point: tuple[float, ...]
    samples: tuple[ProbeSample, ...]
    field: str | None = None
    unit: str | None = None

    def validate_monotonic(self) -> None:
        for a, b in zip(self.samples, self.samples[1:]):
            if b.iteration < a.iteration or b.time < a.time:
                raise ValueError("probe samples must be monotonic")

    def max_relative_error(self, reference: Sequence[float]) -> float:
        if len(reference) != len(self.samples):
            raise ValueError("reference and probe sample lengths differ")
        self.validate_monotonic()
        maximum = 0.0
        for sample, expected in zip(self.samples, reference):
            if not math.isfinite(sample.value) or not math.isfinite(float(expected)):
                raise ValueError("probe validation requires finite values")
            scale = max(abs(float(expected)), 1.0)
            maximum = max(maximum, abs(sample.value - float(expected)) / scale)
        return maximum

    def validate_against(self, reference: Sequence[float], *, atol: float = 1e-12, rtol: float = 1e-8) -> float:
        error = self.max_relative_error(reference)
        for sample, expected in zip(self.samples, reference):
            if abs(sample.value - float(expected)) > atol + rtol * abs(float(expected)):
                raise ValueError(
                    f"probe {self.name!r} exceeds tolerance at iteration {sample.iteration}"
                )
        return error


def _format_probe_csv_value(value: float) -> str:
    """Canonical probe-CSV float token shared with the C++ writer.

    ``repr(float)`` yields the shortest representation that round-trips exactly,
    using Python notation (decimal for 10^-4 <= |v| < 10^16, scientific
    otherwise, with a trailing ``.0`` for integer magnitudes). The C++ writer
    reproduces this verbatim, so the two implementations emit byte-identical
    values. See ``CANONICAL_FLOATS`` / ``test_canonical_float_format_matches_cpp``.
    """
    if not math.isfinite(value):
        raise ValueError("probe CSV values must be finite")
    return repr(float(value))


def _format_probe_metadata_line(
    name: str, point: Sequence[float], field: str, unit: str
) -> str:
    """Serialize a ``# probe`` metadata line byte-for-byte with the C++ writer."""
    coords = ",".join(_format_probe_csv_value(coordinate) for coordinate in point)
    return f"# probe {name} at {coords} field={field} unit={unit}"


def _parse_probe_metadata_line(body: str) -> tuple[str, tuple[float, float, float], str, str] | None:
    """Parse the text after the leading ``#`` of a ``# probe`` metadata line.

    Returns ``(name, (x, y, z), field, unit)`` or ``None`` for anything that does
    not match the ``# probe <name> at <x>,<y>,<z> field=<f> unit=<u>`` schema.
    """
    prefix = "probe "
    if not body.startswith(prefix):
        return None
    rest = body[len(prefix):]
    field_key = " field="
    unit_key = " unit="
    field_index = rest.rfind(field_key)
    unit_index = rest.rfind(unit_key)
    if field_index == -1 or unit_index == -1 or unit_index < field_index:
        return None
    unit = rest[unit_index + len(unit_key):]
    field = rest[field_index + len(field_key):unit_index]
    head = rest[:field_index]
    at_key = " at "
    at_index = head.rfind(at_key)
    if at_index == -1:
        return None
    name = head[:at_index]
    coords = head[at_index + len(at_key):].split(",")
    if len(coords) != 3:
        return None
    try:
        point = tuple(float(coordinate) for coordinate in coords)  # type: ignore[return-value]
    except ValueError:
        return None
    return name, point, field, unit


def write_probe_csv(path: str | Path, series: Sequence["ProbeSeries"]) -> Path:
    """Write probe series as the native probe CSV schema.

    Rows are ordered by probe name then iteration, matching
    ``cfdx::io::write_probe_csv``, so the same run produces the same bytes from
    either language.
    """
    if not series:
        raise ValueError("probe CSV export requires at least one probe series")
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)

    rows: list[tuple[str, int, float, float]] = []
    for probe in series:
        if not probe.name:
            raise ValueError("probe name must not be empty")
        if any(character in probe.name for character in ",\r\n"):
            raise ValueError(
                f"probe name must not contain a comma or newline: {probe.name!r}"
            )
        for sample in probe.samples:
            if sample.iteration < 0:
                raise ValueError("probe iteration must be non-negative")
            if not math.isfinite(sample.time) or not math.isfinite(sample.value):
                raise ValueError(
                    f"probe {probe.name!r} has a non-finite sample"
                )
            rows.append((probe.name, sample.iteration, sample.time, sample.value))
    rows.sort(key=lambda row: (row[0], row[1]))

    lines = [
        f"# {PROBE_CSV_MAGIC} v{PROBE_CSV_VERSION}",
        "# columns: " + ",".join(_PROBE_CSV_COLUMNS),
    ]
    by_name: dict[str, "ProbeSeries"] = {probe.name: probe for probe in series}
    sampled_names = sorted({row[0] for row in rows})
    for name in sampled_names:
        probe = by_name.get(name)
        if probe is not None and probe.point and probe.field is not None and probe.unit is not None:
            lines.append(_format_probe_metadata_line(name, probe.point, probe.field, probe.unit))
    for name, iteration, time_value, value in rows:
        lines.append(
            f"{name},{iteration},{_format_probe_csv_value(time_value)},"
            f"{_format_probe_csv_value(value)}"
        )
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return destination


def read_probe_csv(path: str | Path) -> tuple[ProbeSeries, ...]:
    """Read the native probe CSV schema back into probe series.

    ``# probe`` metadata lines record each probe's location, sampled field and
    unit, so a series comes back with ``point``/``field``/``unit`` populated when
    the file carries them (a legacy file without metadata yields empty values).
    Series are returned in probe-name order.
    """
    source = Path(path)
    text = source.read_text(encoding="utf-8")

    magic: str | None = None
    columns_ok = False
    metadata: dict[str, tuple[tuple[float, ...], str | None, str | None]] = {}
    data_lines: list[str] = []
    for line in text.splitlines():
        if line.startswith("#"):
            body = line[1:].strip()
            if body.startswith(PROBE_CSV_MAGIC) and magic is None:
                parts = body.split()
                if len(parts) != 2 or parts[0] != PROBE_CSV_MAGIC:
                    raise ValueError("malformed probe CSV header")
                try:
                    version = int(parts[1].removeprefix("v"))
                except ValueError as exc:
                    raise ValueError(f"malformed probe CSV header: {body}") from exc
                if version != PROBE_CSV_VERSION:
                    raise ValueError(
                        f"unsupported probe CSV version: {parts[1]}"
                    )
                magic = parts[0]
            elif body.startswith("columns:"):
                expected_columns = "columns: " + ",".join(_PROBE_CSV_COLUMNS)
                if body != expected_columns:
                    raise ValueError(
                        "probe CSV column header mismatch: expected "
                        + ",".join(_PROBE_CSV_COLUMNS)
                    )
                columns_ok = True
            elif body.startswith("probe "):
                parsed = _parse_probe_metadata_line(body)
                if parsed is not None:
                    name, point, field, unit = parsed
                    metadata[name] = (point, field, unit)
            continue
        if line.strip():
            data_lines.append(line)

    if magic is None:
        raise ValueError(f"not a CFDX probe CSV: {source}")
    if not columns_ok:
        raise ValueError("probe CSV is missing the # columns: header")

    grouped: dict[str, list[ProbeSample]] = {}
    for line in data_lines:
        fields = line.split(",")
        if len(fields) != len(_PROBE_CSV_COLUMNS):
            raise ValueError(f"probe CSV row must have four columns: {line!r}")
        name, iteration_text, time_text, value_text = fields
        if not name:
            raise ValueError("probe CSV row has an empty probe name")
        try:
            iteration = int(iteration_text)
        except ValueError as exc:
            raise ValueError(f"probe CSV row has a non-integer iteration: {line!r}") from exc
        if iteration < 0:
            raise ValueError("probe CSV iteration must be non-negative")
        try:
            time_value = float(time_text)
            value = float(value_text)
        except ValueError as exc:
            raise ValueError(f"probe CSV row has a non-numeric value: {line!r}") from exc
        if not math.isfinite(time_value) or not math.isfinite(value):
            raise ValueError(f"probe CSV row has a non-finite value: {line!r}")
        grouped.setdefault(name, []).append(ProbeSample(iteration, time_value, value))

    if not grouped:
        raise ValueError("probe CSV contains no samples")

    series = []
    for name in sorted(grouped):
        point, field, unit = metadata.get(name, ((), None, None))
        probe = ProbeSeries(name, point, tuple(grouped[name]), field=field, unit=unit)
        probe.validate_monotonic()
        series.append(probe)
    return tuple(series)

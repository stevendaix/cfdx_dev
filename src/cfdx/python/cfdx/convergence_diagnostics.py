"""Evidence-first convergence trends for the CFDX monitoring UI.

This module deliberately never evaluates a solver's scientific pass/fail gate.
A trend is a visual diagnostic, computed only from values reported by the real
solver. The independent convergence certificate remains authoritative.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from statistics import median

from .monitors import MonitorSeries

_SPARK = "▁▂▃▄▅▆▇█"


@dataclass(frozen=True, slots=True)
class MetricTrend:
    name: str
    group: str
    last: float
    count: int
    trend: str
    decades_per_iteration: float | None
    sparkline: str
    # This is source telemetry, not a computed physical convergence verdict.
    source: str = "solver telemetry"


def metric_group(name: str) -> str:
    """Do not confuse linear algebra, conservation, or update norms."""
    key = name.casefold()
    if key in {"continuity_normalized", "global_mass_imbalance", "local_mass_imbalance"} or key.startswith("mass_imbalance"):
        return "Conservation"
    if "flux_velocity_mismatch" in key or "conservation" in key:
        return "Conservation"
    if "linear" in key or key == "pressure_relative_residual" or key.startswith("ksp_"):
        return "Linear solve"
    if "change" in key or "update" in key:
        return "State updates"
    if "residual" in key or key.startswith("res_"):
        return "Nonlinear / equation"
    return "Other reported"


def _trend(values: list[tuple[int, float]]) -> tuple[str, float | None]:
    if len(values) < 4:
        return "INSUFFICIENT", None
    # Theil-Sen robust median slope is less sensitive to transient spikes
    # than a least-squares fit. Floors protect log10(0) without claiming a
    # scientifically meaningful residual tolerance.
    points = [(iteration, math.log10(max(abs(value), 1e-300))) for iteration, value in values]
    slopes = [
        (points[j][1] - points[i][1]) / (points[j][0] - points[i][0])
        for i in range(len(points))
        for j in range(i + 1, len(points))
        if points[j][0] > points[i][0]
    ]
    if not slopes:
        return "INSUFFICIENT", None
    slope = float(median(slopes))
    log_values = [y for _, y in points]
    reversals = sum(
        (b - a) * (c - b) < 0
        for a, b, c in zip(log_values, log_values[1:], log_values[2:])
    )
    if reversals >= 3 and max(log_values) - min(log_values) >= 0.5:
        return "OSCILLATING", slope
    if slope > 0.02:
        return "GROWING", slope
    if slope < -0.02:
        return "DECREASING", slope
    return "STAGNATING", slope


def _sparkline(values: list[float]) -> str:
    if not values:
        return ""
    logs = [math.log10(max(abs(value), 1e-300)) for value in values]
    low, high = min(logs), max(logs)
    if high == low:
        return _SPARK[len(_SPARK) // 2] * len(logs)
    return "".join(_SPARK[min(len(_SPARK) - 1, int((v - low) / (high - low) * (len(_SPARK) - 1)))] for v in logs)


def analyze_monitor_series(
    series: MonitorSeries | None,
    *,
    window: int = 12,
    max_metrics: int = 48,
) -> tuple[MetricTrend, ...]:
    """Analyze a bounded recent window, preserving missing-data semantics.

    Interleaved stdout/stderr samples can contain incomplete rows.  A metric
    is only analyzed where it was explicitly reported; no interpolation and
    no reuse of a preceding measurement are permitted.
    """
    if not isinstance(window, int) or isinstance(window, bool) or not 4 <= window <= 256:
        raise ValueError("window must be an integer between 4 and 256")
    if not isinstance(max_metrics, int) or isinstance(max_metrics, bool) or not 1 <= max_metrics <= 1000:
        raise ValueError("max_metrics must be between 1 and 1000")
    if series is None:
        return ()
    collected: dict[str, list[tuple[int, float]]] = {}
    # Bound the work to the most recent 4096 snapshots; the run center has
    # its own full monitoring history for archival export.
    for sample in series.samples[-4096:]:
        if not isinstance(sample.iteration, int) or isinstance(sample.iteration, bool) or sample.iteration < 0:
            raise ValueError("invalid monitor iteration")
        for name, raw in sample.values.items():
            if not isinstance(name, str) or not name:
                raise ValueError("invalid metric name")
            if isinstance(raw, bool) or not isinstance(raw, (float, int)) or not math.isfinite(raw):
                raise ValueError(f"invalid/non-finite solver telemetry: {name}")
            if name == "CFL":
                continue  # CFL is not a residual, and no residual trend applies.
            if name not in collected and len(collected) >= max_metrics:
                continue
            values = collected.setdefault(name, [])
            if values and sample.iteration == values[-1][0]:
                values[-1] = (sample.iteration, float(raw))
            else:
                values.append((sample.iteration, float(raw)))
                if len(values) > window:
                    del values[:-window]
    trends: list[MetricTrend] = []
    for name, values in collected.items():
        if not values:
            continue
        label, slope = _trend(values)
        trends.append(
            MetricTrend(
                name=name,
                group=metric_group(name),
                last=values[-1][1],
                count=len(values),
                trend=label,
                decades_per_iteration=slope,
                sparkline=_sparkline([value for _, value in values]),
            )
        )
    group_order = {v: i for i, v in enumerate((
        "Nonlinear / equation", "Conservation", "Linear solve",
        "State updates", "Other reported",
    ))}
    return tuple(sorted(trends, key=lambda item: (group_order[item.group], item.name)))

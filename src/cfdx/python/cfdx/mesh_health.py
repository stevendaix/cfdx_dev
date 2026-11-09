"""Bounded, verifiable volume-mesh health summaries, independent of Qt.

The native CFDX finite-volume skewness convention is used. Thresholds are
*user-adjustable visual alerts*, never numerical solver accept/reject gates.
Cells without an internal face are unassessed rather than silently scored 0.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .volume_mesh_quality import VolumeSkewnessReport, analyze_cfdx_volume_skewness

_SKEW_BINS = (0.0, 0.1, 0.25, 0.5, 1.0, 2.0, math.inf)
_ANGLE_BINS = (0.0, 10.0, 30.0, 45.0, 60.0, 70.0, 80.0, 90.0, 180.0, math.inf)


@dataclass(frozen=True, slots=True)
class MetricDistribution:
    label: str
    units: str
    assessed: int
    unavailable: int
    minimum: float | None
    p50: float | None
    p95: float | None
    maximum: float | None
    alert_limit: float
    above_alert: int
    intervals: tuple[tuple[str, int], ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "label": self.label,
            "units": self.units,
            "assessed": self.assessed,
            "unavailable": self.unavailable,
            "minimum": self.minimum,
            "p50": self.p50,
            "p95": self.p95,
            "maximum": self.maximum,
            "alert_limit": self.alert_limit,
            "above_alert": self.above_alert,
            "intervals": [
                {"range": interval, "count": count}
                for interval, count in self.intervals
            ],
        }


@dataclass(frozen=True, slots=True)
class MeshHealthSummary:
    cell_count: int
    internal_face_count: int
    skewness: MetricDistribution
    non_orthogonality: MetricDistribution
    worst_cells: tuple[tuple[int, float | None, float | None], ...]
    source_file: str
    file_signature: tuple[int, int, int, int]

    def to_dict(self) -> dict[str, object]:
        return {
            "format": "CFDX-MESH-HEALTH",
            "schema_version": 1,
            "source_file": self.source_file,
            "file_signature": list(self.file_signature),
            "cell_count": self.cell_count,
            "internal_face_count": self.internal_face_count,
            "skewness": self.skewness.to_dict(),
            "non_orthogonality_deg": self.non_orthogonality.to_dict(),
            "worst_cells": [
                {"cell_id": cell, "skewness": skew, "non_orthogonality_deg": angle}
                for cell, skew, angle in self.worst_cells
            ],
            "solver_mesh_acceptance": "not_evaluated",
        }


def _distribution(
    values: np.ndarray,
    *,
    name: str,
    unit: str,
    limit: float,
    bins: tuple[float, ...],
) -> MetricDistribution:
    if (
        not isinstance(limit, (int, float))
        or isinstance(limit, bool)
        or not math.isfinite(limit)
        or limit <= 0
    ):
        raise ValueError(f"{name} alert threshold must be a positive finite number")
    data = np.asarray(values, dtype=float)
    if data.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    # Missing values use NaN, and infinities are never legitimate measurements.
    if np.isinf(data).any():
        raise ValueError(f"{name} contains infinite values")
    valid = data[np.isfinite(data)]
    if np.any(valid < 0):
        raise ValueError(f"{name} contains negative values")
    distribution = []
    for lower, upper in zip(bins[:-1], bins[1:], strict=True):
        count = int(np.count_nonzero((valid >= lower) & (valid < upper)))
        label = f"{lower:g}–{upper:g}" if math.isfinite(upper) else f">={lower:g}"
        distribution.append((label, count))
    return MetricDistribution(
        label=name,
        units=unit,
        assessed=int(valid.size),
        unavailable=int(data.size - valid.size),
        minimum=float(np.min(valid)) if valid.size else None,
        p50=float(np.percentile(valid, 50)) if valid.size else None,
        p95=float(np.percentile(valid, 95)) if valid.size else None,
        maximum=float(np.max(valid)) if valid.size else None,
        alert_limit=float(limit),
        above_alert=int(np.count_nonzero(valid > limit)),
        intervals=tuple(distribution),
    )


def summarize_mesh_health(
    report: VolumeSkewnessReport,
    *,
    source_file: str | Path,
    file_signature: tuple[int, int, int, int],
    skewness_alert: float = 0.5,
    non_orthogonality_alert_deg: float = 70.0,
    top_n: int = 20,
) -> MeshHealthSummary:
    """Summarize only evaluated cells, preserving source identity and bounds."""
    if isinstance(top_n, bool) or not isinstance(top_n, int) or not 1 <= top_n <= 100:
        raise ValueError("top_n must be between 1 and 100")
    if len(file_signature) != 4 or any(not isinstance(v, int) for v in file_signature):
        raise ValueError("invalid source signature")
    scores = np.asarray(report.cell_skewness, dtype=float)
    angles = np.asarray(report.cell_non_orthogonality_deg, dtype=float)
    if scores.shape != angles.shape:
        raise ValueError("mesh metrics must share the same cell index range")
    skew = _distribution(
        scores,
        name="cell maximum internal-face skewness",
        unit="1",
        limit=skewness_alert,
        bins=_SKEW_BINS,
    )
    angle = _distribution(
        angles,
        name="cell maximum face non-orthogonality",
        unit="deg",
        limit=non_orthogonality_alert_deg,
        bins=_ANGLE_BINS,
    )
    evaluated = np.flatnonzero(np.isfinite(scores))
    worst_order = evaluated[np.argsort(scores[evaluated], kind="stable")[::-1]][:top_n]
    worst = tuple(
        (
            int(i),
            float(scores[i]),
            float(angles[i]) if math.isfinite(angles[i]) else None,
        )
        for i in worst_order
    )
    return MeshHealthSummary(
        cell_count=int(scores.size),
        internal_face_count=int(report.internal_face_count),
        skewness=skew,
        non_orthogonality=angle,
        worst_cells=worst,
        source_file=str(source_file),
        file_signature=file_signature,
    )


def analyze_mesh_health_file(
    path: str | Path,
    *,
    expected_signature: tuple[int, int, int, int] | None = None,
    skewness_alert: float = 0.5,
    non_orthogonality_alert_deg: float = 70.0,
) -> MeshHealthSummary:
    """Run the existing guarded HDF5 loader; reject files modified mid-read."""
    source = Path(path)
    stat = source.stat()
    signature = (
        int(stat.st_dev),
        int(stat.st_ino),
        int(stat.st_size),
        int(stat.st_mtime_ns),
    )
    if expected_signature is not None and expected_signature != signature:
        raise RuntimeError("mesh has changed since it was inspected")
    report = analyze_cfdx_volume_skewness(
        source,
        expected_signature=signature,
        threshold=skewness_alert,
    )
    return summarize_mesh_health(
        report,
        source_file=source,
        file_signature=signature,
        skewness_alert=skewness_alert,
        non_orthogonality_alert_deg=non_orthogonality_alert_deg,
    )


def mesh_health_text(summary: MeshHealthSummary) -> str:
    """Readable, source-grounded report for the GUI and CLI."""
    lines = [
        "Volume mesh diagnostic — CFDX native finite-volume conventions",
        f"Cells: {summary.cell_count:,} | internal faces: {summary.internal_face_count:,}",
        "Alert limits are visualization choices, NOT solver quality gates.",
    ]
    for distribution in (summary.skewness, summary.non_orthogonality):
        lines.append(
            f"{distribution.label} [{distribution.units}] | assessed "
            f"{distribution.assessed:,}/{summary.cell_count:,}, "
            f"unavailable {distribution.unavailable:,}; "
            f"p95 {distribution.p95 if distribution.p95 is not None else 'n/a'}; "
            f"max {distribution.maximum if distribution.maximum is not None else 'n/a'}; "
            f"> {distribution.alert_limit:g}: {distribution.above_alert:,}"
        )
        highest_count = max((count for _, count in distribution.intervals), default=0)
        for label, count in distribution.intervals:
            bar = "█" * (round(22 * count / highest_count) if highest_count else 0)
            lines.append(f"  {label:>10} | {count:>7,} {bar}")
    lines.append("Worst assessed cells (cell ID : skewness, non-orthogonality deg):")
    lines.extend(
        f"  {cell}: {skew:.6g}, {angle if angle is not None else 'n/a'}"
        for cell, skew, angle in summary.worst_cells
    )
    lines.append("External physical/mesh qualification: NOT EVALUATED")
    return "\n".join(lines)

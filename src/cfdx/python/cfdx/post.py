"""Small, dependency-free post-processing readers and reductions."""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class MonitorPoint:
    time: float
    values: dict[str, float]


def read_monitor_csv(path: Path) -> tuple[MonitorPoint, ...]:
    """Read a CSV with a time column and numeric monitor columns."""
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or "time" not in reader.fieldnames:
            raise ValueError("monitor CSV must contain a time column")
        points = []
        for row in reader:
            points.append(
                MonitorPoint(
                    time=float(row["time"]),
                    values={
                        key: float(value)
                        for key, value in row.items()
                        if key != "time" and value not in (None, "")
                    },
                )
            )
    return tuple(points)


def reduce_values(values: Iterable[float]) -> dict[str, float]:
    """Return count, min, max and arithmetic mean for a scalar sequence."""
    data = tuple(float(value) for value in values)
    if not data:
        raise ValueError("cannot reduce an empty sequence")
    return {
        "count": float(len(data)),
        "min": min(data),
        "max": max(data),
        "mean": sum(data) / len(data),
    }

"""Simple point probes over structured scalar samples."""
from __future__ import annotations

from collections.abc import Sequence


def nearest_probe(
    coordinates: Sequence[Sequence[float]],
    values: Sequence[float],
    point: Sequence[float],
) -> float:
    """Return the value at the nearest sample point."""
    if len(coordinates) != len(values) or not coordinates:
        raise ValueError("coordinates and values must be non-empty and aligned")
    if not point:
        raise ValueError("probe point must not be empty")
    dimensions = len(point)
    if any(len(coord) != dimensions for coord in coordinates):
        raise ValueError("coordinate dimensions must match the probe")
    index = min(
        range(len(coordinates)),
        key=lambda i: sum((coordinates[i][j] - point[j]) ** 2 for j in range(dimensions)),
    )
    return float(values[index])

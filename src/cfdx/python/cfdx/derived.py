"""Dependency-free scalar post-processing primitives."""
from __future__ import annotations

from collections.abc import Mapping, Sequence


def derived_magnitude(values: Mapping[str, Sequence[float]]) -> tuple[float, ...]:
    """Compute Euclidean magnitude from vector component arrays."""
    if not values:
        raise ValueError("at least one component is required")
    names = tuple(values)
    size = len(values[names[0]])
    if any(len(values[name]) != size for name in names):
        raise ValueError("component lengths must match")
    return tuple(
        sum(float(values[name][i]) ** 2 for name in names) ** 0.5
        for i in range(size)
    )


def derived_sum(*fields: Sequence[float]) -> tuple[float, ...]:
    """Compute pointwise sum of equally-sized scalar fields."""
    if not fields:
        raise ValueError("at least one field is required")
    size = len(fields[0])
    if any(len(field) != size for field in fields):
        raise ValueError("field lengths must match")
    return tuple(sum(float(field[i]) for field in fields) for i in range(size))

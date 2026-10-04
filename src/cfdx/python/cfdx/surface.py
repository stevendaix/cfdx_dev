"""Generic surface reduction primitives for forces and fluxes."""
from __future__ import annotations

from collections.abc import Iterable, Sequence


def surface_integral(values: Sequence[float], areas: Sequence[float]) -> float:
    """Compute a scalar surface integral from face values and face areas."""
    if len(values) != len(areas):
        raise ValueError("values and areas must have equal lengths")
    if not values:
        raise ValueError("surface integral requires at least one face")
    return sum(float(value) * float(area) for value, area in zip(values, areas))


def vector_surface_integral(
    values: Iterable[Sequence[float]], areas: Sequence[float]
) -> tuple[float, ...]:
    """Integrate vector-valued face quantities component-wise."""
    rows = tuple(tuple(float(x) for x in row) for row in values)
    if not rows or len(rows) != len(areas):
        raise ValueError("values and areas must describe the same non-empty faces")
    width = len(rows[0])
    if width == 0 or any(len(row) != width for row in rows):
        raise ValueError("vector components must have consistent width")
    return tuple(
        sum(row[j] * float(area) for row, area in zip(rows, areas))
        for j in range(width)
    )

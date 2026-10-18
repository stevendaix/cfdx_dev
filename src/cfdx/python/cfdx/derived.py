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
        sum(float(values[name][i]) ** 2 for name in names) ** 0.5 for i in range(size)
    )


def derived_vector_magnitude(vector: Sequence[Sequence[float]]) -> tuple[float, ...]:
    """Compute the Euclidean magnitude of each cell's vector value.

    The production solver stores ``U`` as one cell field with three components
    per cell (``Field<double, Location::CELL>(n, "U", "m/s", 3)`` in
    ``apps/production_solver/main.cpp``), so ``vector`` is the per-cell array
    ``[[u_x_0, u_y_0, u_z_0], [u_x_1, u_y_1, u_z_1], ...]`` and the result has
    one magnitude per cell. Computing from the vector array directly (rather
    than from separately named component scalars) keeps the derived field
    aligned with how the dataset actually stores the velocity.
    """
    if not vector:
        raise ValueError("at least one cell is required")
    components = tuple(vector)
    if any(len(cell) != len(components[0]) for cell in components):
        raise ValueError("vector component lengths must match")
    return tuple(sum(float(value) ** 2 for value in cell) ** 0.5 for cell in components)


def derived_sum(*fields: Sequence[float]) -> tuple[float, ...]:
    """Compute pointwise sum of equally-sized scalar fields."""
    if not fields:
        raise ValueError("at least one field is required")
    size = len(fields[0])
    if any(len(field) != size for field in fields):
        raise ValueError("field lengths must match")
    return tuple(sum(float(field[i]) for field in fields) for i in range(size))

"""Validation helpers for bounded in-memory mesh preview topology."""

from __future__ import annotations

from typing import Any

import numpy as np


def _validated_offsets(
    offsets: Any, *, name: str, expected_count: int, terminal: int
) -> np.ndarray:
    values = np.asarray(offsets)
    if values.ndim != 1 or values.size != expected_count + 1:
        raise ValueError(f"{name} must contain {expected_count + 1} offsets")
    if values.size == 0 or int(values[0]) != 0:
        raise ValueError(f"{name} must start at zero")
    if np.any(values[1:] < values[:-1]):
        raise ValueError(f"{name} must be monotonically non-decreasing")
    if int(values[-1]) != terminal:
        raise ValueError(f"{name} terminal offset must equal {terminal}")
    return values


def validate_surface_topology(
    points: Any,
    face_vertices: Any,
    face_offsets: Any,
    patch_face_ids: Any,
    patch_offsets: Any,
    *,
    n_faces: int,
    n_patches: int,
) -> None:
    """Validate all CSR ranges before passing arrays to VTK/PyVista."""
    coordinates = np.asarray(points)
    vertices = np.asarray(face_vertices)
    patch_ids = np.asarray(patch_face_ids)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("points must have shape (n_points, 3)")
    if not np.isfinite(coordinates).all():
        raise ValueError("point coordinates must all be finite")
    if vertices.ndim != 1:
        raise ValueError("face connectivity must be one-dimensional")
    if vertices.size and (
        int(vertices.min()) < 0 or int(vertices.max()) >= coordinates.shape[0]
    ):
        raise ValueError("face connectivity contains an out-of-range point index")

    offsets = _validated_offsets(
        face_offsets,
        name="face_offsets",
        expected_count=n_faces,
        terminal=int(vertices.size),
    )
    if np.any(np.diff(offsets) < 3):
        raise ValueError("every surface face must contain at least three vertices")

    if patch_ids.ndim != 1:
        raise ValueError("patch_face_ids must be one-dimensional")
    if patch_ids.size and (
        int(patch_ids.min()) < 0 or int(patch_ids.max()) >= n_faces
    ):
        raise ValueError("patch_face_ids contains an out-of-range face index")
    _validated_offsets(
        patch_offsets,
        name="patch_face_offsets",
        expected_count=n_patches,
        terminal=int(patch_ids.size),
    )


def validate_cell_topology(
    cell_faces: Any,
    cell_offsets: Any,
    face_offsets: Any,
    *,
    n_cells: int,
    n_faces: int,
    max_vtk_entries: int = 2_000_000,
) -> None:
    """Validate cell-to-face CSR and bound the expanded VTK cell array."""
    face_ids = np.asarray(cell_faces)
    if face_ids.ndim != 1:
        raise ValueError("cell_faces must be one-dimensional")
    if face_ids.size and (
        int(face_ids.min()) < 0 or int(face_ids.max()) >= n_faces
    ):
        raise ValueError("cell_faces contains an out-of-range face index")
    offsets = _validated_offsets(
        cell_offsets,
        name="cell_offsets",
        expected_count=n_cells,
        terminal=int(face_ids.size),
    )
    if np.any(np.diff(offsets) < 1):
        raise ValueError("every cell must reference at least one face")

    face_sizes = np.diff(np.asarray(face_offsets))
    expanded_entries = n_cells
    for start in range(0, face_ids.size, 65_536):
        referenced = face_ids[start : start + 65_536]
        expanded_entries += int(referenced.size + face_sizes[referenced].sum())
        if expanded_entries > max_vtk_entries:
            raise ValueError(
                "expanded VTK cell connectivity exceeds the 2M-entry safety budget"
            )

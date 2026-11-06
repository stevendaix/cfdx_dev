"""Malformed HDF5 connectivity must be rejected before VTK construction."""

from __future__ import annotations

import numpy as np
import pytest
from cfdx.mesh_topology import validate_cell_topology, validate_surface_topology


def _valid_surface():
    return (
        np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        np.asarray([0, 1, 2], dtype=np.int64),
        np.asarray([0, 3], dtype=np.int64),
        np.asarray([0], dtype=np.int64),
        np.asarray([0, 1], dtype=np.int64),
    )


def test_valid_surface_topology_is_accepted():
    validate_surface_topology(
        *_valid_surface(),
        n_faces=1,
        n_patches=1,
    )


@pytest.mark.parametrize(
    ("array_index", "replacement", "message"),
    [
        (1, np.asarray([0, 1, 2, 99]), "out-of-range point index"),
        (2, np.asarray([0, 2]), "terminal offset"),
        (2, np.asarray([0, 4]), "terminal offset"),
        (2, np.asarray([1, 3]), "start at zero"),
        (3, np.asarray([1]), "out-of-range face index"),
        (4, np.asarray([0, 2]), "terminal offset"),
    ],
)
def test_invalid_surface_connectivity_is_rejected(array_index, replacement, message):
    arrays = list(_valid_surface())
    arrays[array_index] = replacement
    with pytest.raises(ValueError, match=message):
        validate_surface_topology(*arrays, n_faces=1, n_patches=1)


def test_surface_rejects_non_finite_coordinates_and_degenerate_faces():
    points, vertices, face_offsets, patch_ids, patch_offsets = _valid_surface()
    points[0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        validate_surface_topology(
            points,
            vertices,
            face_offsets,
            patch_ids,
            patch_offsets,
            n_faces=1,
            n_patches=1,
        )

    with pytest.raises(ValueError, match="at least three"):
        validate_surface_topology(
            _valid_surface()[0],
            np.asarray([0, 1]),
            np.asarray([0, 2]),
            patch_ids,
            patch_offsets,
            n_faces=1,
            n_patches=1,
        )


def test_valid_cell_topology_is_accepted():
    validate_cell_topology(
        np.asarray([0, 1, 2, 3], dtype=np.int64),
        np.asarray([0, 4], dtype=np.int64),
        np.asarray([0, 3, 6, 9, 12], dtype=np.int64),
        n_cells=1,
        n_faces=4,
    )


@pytest.mark.parametrize(
    ("cell_faces", "cell_offsets", "message"),
    [
        ([0, 4], [0, 2], "out-of-range face index"),
        ([0, 1], [0, 1], "terminal offset"),
        ([0, 1], [1, 2], "start at zero"),
        ([0, 1], [0, 2, 1], "monotonically"),
        ([], [0, 0], "at least one face"),
    ],
)
def test_invalid_cell_connectivity_is_rejected(cell_faces, cell_offsets, message):
    with pytest.raises(ValueError, match=message):
        validate_cell_topology(
            np.asarray(cell_faces, dtype=np.int64),
            np.asarray(cell_offsets, dtype=np.int64),
            np.asarray([0, 3, 6, 9, 12], dtype=np.int64),
            n_cells=len(cell_offsets) - 1,
            n_faces=4,
        )


def test_cell_topology_bounds_expanded_vtk_connectivity():
    with pytest.raises(ValueError, match="expanded VTK cell connectivity"):
        validate_cell_topology(
            np.zeros(1000, dtype=np.int64),
            np.asarray([0, 1000], dtype=np.int64),
            np.asarray([0, 2000], dtype=np.int64),
            n_cells=1,
            n_faces=1,
            max_vtk_entries=1000,
        )

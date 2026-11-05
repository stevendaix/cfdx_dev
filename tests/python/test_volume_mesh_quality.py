"""Regression tests for volumetric skewness and cell adjacency."""

import h5py
import numpy as np
import pytest
from cfdx.volume_mesh_quality import (
    analyze_cfdx_volume_skewness,
    analyze_volume_skewness,
)


def two_tetrahedra():
    points = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 2.0],
        ]
    )
    # Shared face (1,2,3), plus three exterior faces per tetrahedron.
    faces = [
        [1, 2, 3],
        [0, 2, 1],
        [0, 1, 3],
        [0, 3, 2],
        [4, 2, 1],
        [4, 1, 3],
        [4, 3, 2],
    ]
    face_vertices = np.asarray([v for face in faces for v in face], dtype=np.int64)
    face_offsets = np.arange(0, 3 * (len(faces) + 1), 3, dtype=np.int64)
    owner = np.asarray([0, 0, 0, 0, 1, 1, 1], dtype=np.int64)
    neighbour = np.asarray([1, -1, -1, -1, -1, -1, -1], dtype=np.int64)
    cell_faces = np.asarray([0, 1, 2, 3, 0, 4, 5, 6], dtype=np.int64)
    cell_offsets = np.asarray([0, 4, 8], dtype=np.int64)
    return (
        points,
        face_vertices,
        face_offsets,
        owner,
        neighbour,
        cell_faces,
        cell_offsets,
    )


def test_volume_skewness_reports_internal_face_and_symmetric_neighbours():
    report = analyze_volume_skewness(*two_tetrahedra(), threshold=0.01)
    assert report.internal_face_count == 1
    assert report.cell_neighbours == ((1,), (0,))
    assert report.face_skewness.shape == (7,)
    assert report.face_non_orthogonality_deg.shape == (7,)
    assert np.isfinite(report.cell_centres).all()
    assert report.cell_non_orthogonality_deg.shape == (2,)
    assert np.isfinite(report.cell_non_orthogonality_deg).all()
    assert np.isnan(report.face_skewness[1:]).all()
    assert report.cell_skewness[0] == pytest.approx(report.cell_skewness[1])
    assert report.bad_cell_ids == (0, 1)


def test_volume_skewness_rejects_invalid_owner_neighbour_indices():
    arrays = list(two_tetrahedra())
    arrays[3] = np.asarray([4, 0, 0, 0, 1, 1, 1])
    with pytest.raises(ValueError, match="owner cell index"):
        analyze_volume_skewness(*arrays)


def test_volume_skewness_rejects_inconsistent_cell_face_incidence():
    arrays = list(two_tetrahedra())
    arrays[5] = np.asarray([0, 1, 2, 3, 1, 4, 5, 6], dtype=np.int64)
    with pytest.raises(ValueError, match="face list disagrees"):
        analyze_volume_skewness(*arrays)


def test_volume_skewness_rejects_boundary_only_mesh():
    arrays = list(two_tetrahedra())
    arrays[1] = arrays[1][:12]
    arrays[2] = np.asarray([0, 3, 6, 9, 12], dtype=np.int64)
    arrays[3] = np.zeros(4, dtype=np.int64)
    arrays[4] = np.full(4, -1, dtype=np.int64)
    arrays[5] = np.asarray([0, 1, 2, 3], dtype=np.int64)
    arrays[6] = np.asarray([0, 4], dtype=np.int64)
    with pytest.raises(ValueError, match="no internal faces"):
        analyze_volume_skewness(*arrays)


def test_hdf5_volume_skewness_preflights_and_rechecks_file(tmp_path):
    path = tmp_path / "two-tets.cfdx.h5"
    names = (
        "points",
        "face_vertices",
        "face_offsets",
        "owner",
        "neighbour",
        "cell_faces",
        "cell_offsets",
    )
    with h5py.File(path, "w") as h5:
        for name, values in zip(names, two_tetrahedra(), strict=True):
            h5.create_dataset(name, data=values)
    report = analyze_cfdx_volume_skewness(path)
    assert report.internal_face_count == 1
    stat = path.stat()
    signature = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns)
    with pytest.raises(RuntimeError, match="changed since metadata"):
        analyze_cfdx_volume_skewness(path, expected_signature=signature[:-1] + (0,))


def test_volume_quality_uses_cfdx_core_skewness_normalization():
    report = analyze_volume_skewness(*two_tetrahedra(), threshold=0.01)
    owner, neighbour = report.cell_centres
    line = neighbour - owner
    face_centre = np.asarray([1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0])
    projection = np.dot(face_centre - owner, line) / np.dot(line, line)
    closest = owner + projection * line
    expected = np.linalg.norm(face_centre - closest) / np.linalg.norm(line)
    assert report.face_skewness[0] == pytest.approx(expected)


def test_hdf5_volume_limits_are_checked_before_payload_reads(tmp_path):
    path = tmp_path / "too-many-faces.h5"
    with h5py.File(path, "w") as h5:
        h5.create_dataset("points", data=np.zeros((4, 3), dtype=np.float64))
        h5.create_dataset("face_vertices", data=np.zeros(3, dtype=np.int64))
        h5.create_dataset("face_offsets", data=np.asarray([0, 3], dtype=np.int64))
        h5.create_dataset(
            "owner", shape=(500_001,), dtype=np.int8, chunks=(1024,)
        )
        h5.create_dataset(
            "neighbour", shape=(500_001,), dtype=np.int8, chunks=(1024,)
        )
        h5.create_dataset("cell_faces", data=np.asarray([0, 0, 0, 0]))
        h5.create_dataset("cell_offsets", data=np.asarray([0, 4]))
    with pytest.raises(ValueError, match="500,000-face eager-load limit"):
        analyze_cfdx_volume_skewness(path)


def test_hdf5_volume_rejects_non_integer_connectivity_before_read(tmp_path):
    path = tmp_path / "float-connectivity.h5"
    arrays = two_tetrahedra()
    names = (
        "points",
        "face_vertices",
        "face_offsets",
        "owner",
        "neighbour",
        "cell_faces",
        "cell_offsets",
    )
    with h5py.File(path, "w") as h5:
        for name, values in zip(names, arrays, strict=True):
            if name == "cell_faces":
                values = values.astype(np.float64)
            h5.create_dataset(name, data=values)
    with pytest.raises(ValueError, match="cell_faces must be a one-dimensional integer"):
        analyze_cfdx_volume_skewness(path)


@pytest.mark.skipif(
    __import__("importlib").util.find_spec("PySide6") is None
    or __import__("importlib").util.find_spec("pyvista") is None
    or __import__("importlib").util.find_spec("pyvistaqt") is None,
    reason="optional Qt/PyVista runtime",
)
def test_volume_skewness_dialog_moves_cut_and_slices_local_neighbours(tmp_path):
    from cfdx.gui import create_application
    from cfdx.volume_skewness_view import VolumeSkewnessDialog
    from PySide6.QtWidgets import QSlider

    path = tmp_path / "two-tets.cfdx.h5"
    names = (
        "points",
        "face_vertices",
        "face_offsets",
        "owner",
        "neighbour",
        "cell_faces",
        "cell_offsets",
    )
    with h5py.File(path, "w") as h5:
        for name, values in zip(names, two_tetrahedra(), strict=True):
            h5.create_dataset(name, data=values)

    app = create_application(["cfdx-volume-skewness-dialog-test"])
    dialog = VolumeSkewnessDialog(path)
    slider = dialog.findChild(QSlider, "volume_skewness.slice_position")
    assert slider is not None
    slider.setValue(250)
    bounds = dialog._grid.bounds
    expected_z = bounds[4] + 0.25 * (bounds[5] - bounds[4])
    assert dialog._cut_axis == "z"
    assert dialog._cut_origin[2] == pytest.approx(expected_z)
    assert dialog._cell_selector.count() == 2
    dialog._cell_selector.setCurrentIndex(0)
    assert dialog._details.text().startswith("Cell ")
    dialog.close()
    app.processEvents()
    app.quit()

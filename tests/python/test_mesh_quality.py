"""Surface quality indicators are bounded and geometrically explicit."""

from __future__ import annotations

import h5py
import numpy as np
import pytest

from cfdx.mesh_quality import (
    analyze_cfdx_surface_file,
    analyze_surface_quality,
)


def test_planar_square_has_unit_area_zero_warpage_and_equal_edges():
    points = np.asarray(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
    )
    report = analyze_surface_quality(
        points, np.asarray([0, 1, 2, 3]), np.asarray([0, 4])
    )

    assert report.face_count == 1
    assert report.valid_face_count == 1
    assert report.degenerate_face_ids == ()
    assert report.face_area_min == pytest.approx(1.0)
    assert report.face_area_max == pytest.approx(1.0)
    assert report.edge_length_ratio_min == pytest.approx(1.0)
    assert report.edge_length_ratio_max == pytest.approx(1.0)
    assert report.warpage_ratio_max == pytest.approx(0.0, abs=1e-14)
    assert report.worst_warpage_face_ids == (0,)


def test_warped_quad_reports_nonzero_plane_deviation():
    points = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.2],
            [0.0, 1.0, 0.0],
        ]
    )
    report = analyze_surface_quality(
        points, np.asarray([0, 1, 2, 3]), np.asarray([0, 4])
    )

    assert report.valid_face_count == 1
    assert report.warpage_ratio_max > 0.0
    assert report.warpage_ratio_p95 == pytest.approx(report.warpage_ratio_max)


def test_degenerate_faces_and_zero_length_edges_are_reported():
    points = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )
    vertices = np.asarray([0, 1, 2, 0, 0, 3])
    offsets = np.asarray([0, 3, 6])
    report = analyze_surface_quality(points, vertices, offsets)

    assert report.face_count == 2
    assert report.valid_face_count == 0
    assert report.degenerate_face_ids == (0, 1)
    assert report.zero_length_edge_face_ids == (1,)
    assert report.face_area_min is None
    assert report.warpage_ratio_max is None


def test_quality_metrics_are_scale_invariant_for_small_meshes():
    points = 1e-9 * np.asarray(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
    )
    report = analyze_surface_quality(
        points, np.asarray([0, 1, 2, 3]), np.asarray([0, 4])
    )

    assert report.valid_face_count == 1
    assert report.face_area_min == pytest.approx(1e-18)
    assert report.edge_length_ratio_max == pytest.approx(1.0)


@pytest.mark.parametrize(
    ("points", "vertices", "offsets", "message"),
    [
        (np.zeros((3, 2)), [0, 1, 2], [0, 3], "shape"),
        (np.zeros((3, 3)), [0.0, 1.0, 2.0], [0, 3], "integer"),
        (np.zeros((3, 3)), [0, 1, 2], [1, 3], "start at zero"),
        (np.zeros((3, 3)), [0, 1, 2], [0, 2], "at least three"),
        (np.zeros((3, 3)), [0, 1, 4], [0, 3], "out-of-range"),
    ],
)
def test_invalid_surface_arrays_fail_closed(points, vertices, offsets, message):
    with pytest.raises(ValueError, match=message):
        analyze_surface_quality(
            np.asarray(points), np.asarray(vertices), np.asarray(offsets)
        )


def test_surface_quality_limits_are_enforced():
    points = np.asarray(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    )
    with pytest.raises(ValueError, match="limits must be positive"):
        analyze_surface_quality(
            points, np.asarray([0, 1, 2]), np.asarray([0, 3]), max_faces=0
        )
    with pytest.raises(ValueError, match="connectivity"):
        analyze_surface_quality(
            points,
            np.asarray([0, 1, 2]),
            np.asarray([0, 3]),
            max_connectivity_entries=2,
        )



def test_hdf5_surface_analysis_respects_inspected_file_identity(tmp_path):
    path = tmp_path / "surface.cfdx.h5"
    with h5py.File(path, "w") as h5:
        h5.create_dataset(
            "points",
            data=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        )
        h5.create_dataset("face_vertices", data=[0, 1, 2, 3])
        h5.create_dataset("face_offsets", data=[0, 4])

    stat = path.stat()
    signature = (
        int(stat.st_dev),
        int(stat.st_ino),
        int(stat.st_size),
        int(stat.st_mtime_ns),
    )
    report = analyze_cfdx_surface_file(path, expected_signature=signature)
    assert report.valid_face_count == 1
    assert report.face_area_min == pytest.approx(1.0)

    with pytest.raises(RuntimeError, match="changed since metadata inspection"):
        analyze_cfdx_surface_file(path, expected_signature=(0, 0, 0, 0))


def test_hdf5_surface_analysis_refuses_oversized_chunks_and_payload(tmp_path):
    path = tmp_path / "surface-limits.h5"
    with h5py.File(path, "w") as h5:
        h5.create_dataset("points", shape=(10, 3), dtype="f8", chunks=(10, 3))
        h5.create_dataset("face_vertices", data=[0, 1, 2])
        h5.create_dataset("face_offsets", data=[0, 3])

    with pytest.raises(ValueError, match="chunk exceeds"):
        analyze_cfdx_surface_file(path, max_block_bytes=64)

    with pytest.raises(ValueError, match="budget"):
        analyze_cfdx_surface_file(path, max_source_bytes=32)

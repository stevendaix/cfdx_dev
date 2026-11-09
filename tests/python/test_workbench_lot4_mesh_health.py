import numpy as np
import pytest
from cfdx.mesh_health import (
    analyze_mesh_health_file,
    mesh_health_text,
    summarize_mesh_health,
)
from cfdx.volume_mesh_quality import VolumeSkewnessReport


def _report():
    return VolumeSkewnessReport(
        cell_skewness=np.array([0.1, 0.6, np.nan, 1.2]),
        face_skewness=np.array([0.1, 0.6, 1.2]),
        face_non_orthogonality_deg=np.array([5.0, 79.0, 88.0]),
        cell_non_orthogonality_deg=np.array([5.0, 79.0, np.nan, 88.0]),
        cell_centres=np.zeros((4, 3)),
        cell_neighbours=((), (), (), ()),
        worst_cell_ids=(3, 1, 0),
        bad_cell_ids=(1, 3),
        threshold=0.5,
        internal_face_count=3,
        max_skewness=1.2,
        p95_skewness=1.14,
    )


def test_mesh_health_preserves_unassessed_cells_and_reports_histograms():
    snapshot = summarize_mesh_health(
        _report(), source_file="mesh.cfdx.h5", file_signature=(1, 2, 3, 4)
    )
    assert snapshot.cell_count == 4
    assert snapshot.skewness.assessed == 3
    assert snapshot.skewness.unavailable == 1
    assert snapshot.skewness.above_alert == 2
    assert snapshot.non_orthogonality.above_alert == 2
    assert snapshot.non_orthogonality.maximum == 88
    assert sum(count for _, count in snapshot.skewness.intervals) == 3
    assert snapshot.worst_cells[0][0] == 3
    assert snapshot.to_dict()["solver_mesh_acceptance"] == "not_evaluated"
    assert "not volume" not in mesh_health_text(snapshot)
    assert "NOT EVALUATED" in mesh_health_text(snapshot)


@pytest.mark.parametrize("limit", [0.0, float("nan"), float("inf"), -3])
def test_bad_user_quality_threshold_never_accepted(limit):
    with pytest.raises(ValueError, match="threshold"):
        summarize_mesh_health(
            _report(),
            source_file="mesh",
            file_signature=(0, 1, 2, 3),
            skewness_alert=limit,
        )


def test_negative_and_infinite_mesh_metrics_are_rejected():
    original = _report()
    original.cell_skewness[0] = float("inf")
    with pytest.raises(ValueError, match="infinite"):
        summarize_mesh_health(original, source_file="mesh", file_signature=(0, 1, 2, 3))
    original.cell_skewness[0] = -0.1
    with pytest.raises(ValueError, match="negative"):
        summarize_mesh_health(original, source_file="mesh", file_signature=(0, 1, 2, 3))


def test_hdf5_health_snapshot_is_from_real_mesh(tmp_path):
    h5py = pytest.importorskip("h5py")
    from test_volume_mesh_quality import two_tetrahedra

    path = tmp_path / "two.cfdx.h5"
    names = (
        "points",
        "face_vertices",
        "face_offsets",
        "owner",
        "neighbour",
        "cell_faces",
        "cell_offsets",
    )
    with h5py.File(path, "w") as handle:
        for name, array in zip(names, two_tetrahedra(), strict=True):
            handle.create_dataset(name, data=array)
    result = analyze_mesh_health_file(path)
    assert result.cell_count == 2
    assert result.internal_face_count == 1
    assert result.skewness.assessed == 2
    assert result.skewness.maximum == pytest.approx(result.worst_cells[0][1])
    st = path.stat()
    with pytest.raises(RuntimeError, match="changed"):
        analyze_mesh_health_file(
            path, expected_signature=(st.st_dev, st.st_ino, st.st_size, 0)
        )

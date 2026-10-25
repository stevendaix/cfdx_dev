"""Executable topology regressions for the NASA/TMR NACA C-grid importer."""
from pathlib import Path
import gzip
import sys

import h5py
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from import_tmr_naca0012 import convert, read_plot3d_2d, main
from generate_2d_validation_meshes import write_mesh


def miniature_c_grid():
    # Two paired wake edges on each side, a four-edge airfoil, one radial
    # layer. This is independent of NASA's ni/nj ratios and is deliberately
    # small enough for topology counts to be an exact independent oracle.
    x = np.array([[3, 2, 1, .5, 0, .5, 1, 2, 3],
                  [3, 2, 1, 0, -1, 0, 1, 2, 3]], dtype=float)
    y = np.array([[0, 0, 0, -.1, 0, .1, 0, 0, 0],
                  [-1, -1, -1, -1, 0, 1, 1, 1, 1]], dtype=float)
    return x, y


def plot3d_text(x, y):
    nj, ni = x.shape
    return f"1\n{ni} {nj}\n" + "\n".join(str(v) for v in np.r_[x.ravel(), y.ravel()])


@pytest.mark.parametrize("compressed", [False, True])
def test_reads_official_single_block_header_and_gzip(tmp_path, compressed):
    x, y = miniature_c_grid()
    path = tmp_path / ("n0012_9-2.p2dfmt.gz" if compressed else "n0012_9-2.p2dfmt")
    text = plot3d_text(x, y)
    if compressed:
        with gzip.open(path, "wt", encoding="ascii") as f:
            f.write(text)
    else:
        path.write_text(text, encoding="ascii")
    actual_x, actual_y = read_plot3d_2d(path)
    np.testing.assert_array_equal(actual_x, x)
    np.testing.assert_array_equal(actual_y, y)


@pytest.mark.parametrize("mutation,diagnostic", [
    (lambda s: s.replace("1\n", "2\n", 1), "single-block"),
    (lambda s: s + "\n123", "numeric tokens"),
    (lambda s: s.replace("3.0", "nan", 1), "finite"),
])
def test_rejects_unsupported_or_malformed_plot3d(tmp_path, mutation, diagnostic):
    x, y = miniature_c_grid()
    path = tmp_path / "n0012_9-2.p2dfmt"
    path.write_text(mutation(plot3d_text(x, y)), encoding="ascii")
    with pytest.raises(ValueError, match=diagnostic):
        read_plot3d_2d(path)


def test_wake_is_internal_and_extrusion_is_empty(tmp_path):
    points, cells, boundaries = convert(*miniature_c_grid())
    assert len(points) == 30  # 18 plane points minus 3 duplicate pairs, twice.
    assert len(cells) == 8
    assert {name: len(faces) for name, faces in boundaries.items()} == {
        "airfoil": 4, "farfield": 10, "front": 8, "back": 8,
    }
    path = tmp_path / "miniature.h5"
    write_mesh(path, points, cells, boundaries)
    with h5py.File(path, "r") as f:
        neighbours = f["neighbour"][:]
        assert np.count_nonzero(neighbours >= 0) == 9  # 7 adjacent + 2 joined wake.
        assert np.count_nonzero(neighbours == -1) == 30
        patch_ids = f["patch_face_ids"][:]
        assert len(np.unique(patch_ids)) == len(patch_ids)
        assert set(patch_ids) == set(np.flatnonzero(neighbours == -1))
        patch_types = {v.split(":")[0]: int(v.split(":")[-1])
                       for v in f.attrs["boundary_patches"].split(";")}
        assert patch_types == {"airfoil": 0, "farfield": 7, "front": 6, "back": 6}
        points = f["points"][:]
        face_vertices, face_offsets = f["face_vertices"][:], f["face_offsets"][:]
        centres = np.array([points[face_vertices[face_offsets[i]:face_offsets[i+1]]].mean(axis=0)
                            for i in range(len(neighbours))])
        on_wake = (centres[:, 0] > 1) & (np.abs(centres[:, 1]) < 1e-14) & (centres[:, 2] == .5)
        assert np.count_nonzero(on_wake) == 2
        assert (neighbours[on_wake] >= 0).all()
        # Signed polygon area need not share a winding convention; its
        # nonzero magnitude is independently checked for every input cell.
        for cell in cells:
            q = points[np.asarray(cell[:4]), :2]
            area = .5 * np.sum(q[:, 0] * np.roll(q[:, 1], -1) - q[:, 1] * np.roll(q[:, 0], -1))
            assert abs(area) > 0


def test_rejects_nonmatching_wake():
    x, y = miniature_c_grid()
    y[0, -1] = .01
    with pytest.raises(ValueError, match="coincident"):
        convert(x, y)


def test_command_records_provenance_without_physical_promotion(tmp_path, monkeypatch):
    source = tmp_path / "n0012_9-2.p2dfmt"
    output = tmp_path / "grid.h5"
    source.write_text(plot3d_text(*miniature_c_grid()), encoding="ascii")
    monkeypatch.setattr(sys, "argv", ["import_tmr_naca0012", "--input", str(source), "--output", str(output)])
    assert main() == 0
    with h5py.File(output, "r") as f:
        assert f.attrs["validation_wake_pairs"] == 2
        assert f.attrs["validation_wall_faces"] == 4
        assert f.attrs["validation_status"] == "imported-provenance-only"
        assert len(f.attrs["validation_source_sha256"]) == 64


def test_partial_wake_corruption_cannot_become_solid_wall():
    x, y = miniature_c_grid()
    y[0, -3] = .01  # break only the trailing-edge pair, after two good pairs.
    with pytest.raises(ValueError, match="trailing edge"):
        convert(x, y)

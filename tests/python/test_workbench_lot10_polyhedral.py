"""LOT10 evidence: true bounded volume tiles and convex polyhedral conservative transfer."""

from __future__ import annotations

import json

import h5py
import numpy as np
import pytest
from cfdx.convex_polyhedral_remap import conservative_convex_remap
from cfdx.simulation_compare import IncompatibleSimulations
from cfdx.volume_partitions import read_volume_partition
from test_workbench_lot9_industrial import make_cartesian_case


def _sheared(path):
    with h5py.File(path, "r+") as h5:
        points = h5["points"][:]
        # A genuine affine skew: non-Cartesian convex polyhedra, exact volumes.
        points[:, 1] += 0.27 * points[:, 0]
        points[:, 2] += 0.13 * points[:, 0] + 0.07 * points[:, 1]
        h5["points"][:] = points
    return path


def _tetra_case(path, vertices, cells, values):
    """Write independent closed tetrahedra, with actual triangular face topology."""
    global_points = []
    faces = []
    fo = [0]
    cell_faces = []
    co = [0]
    for cell in cells:
        offset = len(global_points)
        global_points.extend([vertices[k] for k in cell])
        for face in ((0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)):
            cell_faces.append(len(fo) - 1)
            faces.extend([offset + f for f in face])
            fo.append(len(faces))
        co.append(len(cell_faces))
    with h5py.File(path, "w") as h5:
        for name, data in {
            "points": global_points,
            "face_offsets": fo,
            "face_vertices": faces,
            "cell_faces": cell_faces,
            "cell_offsets": co,
        }.items():
            h5.create_dataset(
                name,
                data=np.asarray(data, dtype=float if name == "points" else np.int64),
            )
        h5.attrs["physical_time"] = 2.0
        d = h5.create_dataset("fields/scalar/T", data=np.asarray(values, dtype=float))
        d.attrs["units"] = "K"
    return path


def test_convex_remap_affine_sheared_hexahedra(tmp_path):
    a = _sheared(
        make_cartesian_case(
            tmp_path / "a.h5",
            [(0, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 1)],
            [1, 3, 5, 7],
        )
    )
    b = _sheared(make_cartesian_case(tmp_path / "b.h5", [(0, 0.5), (0.5, 1)], [2, 6]))
    mapped, report = conservative_convex_remap(a, b, "/fields/scalar/T")
    np.testing.assert_allclose(mapped, [2, 6], atol=1e-8)
    assert report.relative_conservation_error < 1e-8
    assert report.method == "convex_polyhedron_overlap_first_order_qhull"
    assert not report.certificate
    with pytest.raises(IncompatibleSimulations, match="non-Cartesian"):
        from cfdx.conservative_remap import conservative_box_remap

        conservative_box_remap(a, b, "/fields/scalar/T")


def test_skewed_tetrahedron_to_four_tetrahedra_conserves_mass(tmp_path):
    verts = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.2], [0.0, 1.0, 0.1], [0.1, 0.1, 1.5]]
    inner = np.mean(verts, axis=0).tolist()
    a = _tetra_case(tmp_path / "src.h5", verts, [(0, 1, 2, 3)], [42.0])
    b = _tetra_case(
        tmp_path / "dst.h5",
        verts + [inner],
        [(4, 1, 2, 3), (0, 4, 2, 3), (0, 1, 4, 3), (0, 1, 2, 4)],
        [42.0, 42.0, 42.0, 42.0],
    )
    mapped, report = conservative_convex_remap(a, b, "/fields/scalar/T")
    np.testing.assert_allclose(mapped, [42.0] * 4, rtol=1e-8)
    assert report.relative_conservation_error < 1e-8
    back, backwards = conservative_convex_remap(b, a, "/fields/scalar/T")
    assert back.tolist() == pytest.approx([42.0], rel=1e-8)
    assert backwards.relative_conservation_error < 1e-8


def test_rejects_incomplete_convex_domain(tmp_path):
    a = make_cartesian_case(tmp_path / "a.h5", [(0, 0.5), (0.5, 1)], [1.0, 2.0])
    b = make_cartesian_case(tmp_path / "b.h5", [(0, 0.5), (0.6, 1)], [1.0, 2.0])
    with pytest.raises(
        IncompatibleSimulations, match="unequal domain volume|not fully covered"
    ):
        conservative_convex_remap(a, b, "/fields/scalar/T")


def test_rejects_invalid_open_polyhedral_cell(tmp_path):
    path = make_cartesian_case(tmp_path / "bad.h5", [(0, 1)], [1.0])
    other = make_cartesian_case(tmp_path / "good.h5", [(0, 1)], [1.0])
    with h5py.File(path, "r+") as h5:
        # Keep the same offsets, but duplicate a face (breaks edge manifold).
        h5["cell_faces"][5] = h5["cell_faces"][0]
    with pytest.raises(IncompatibleSimulations, match="repeated|open|nonmanifold"):
        conservative_convex_remap(path, other, "/fields/scalar/T")


def test_convex_rejects_nonplanar_and_symlink(tmp_path):
    a = make_cartesian_case(tmp_path / "a.h5", [(0, 1)], [1.0])
    b = make_cartesian_case(tmp_path / "b.h5", [(0, 1)], [1.0])
    alias = tmp_path / "alias.h5"
    alias.symlink_to(a)
    with pytest.raises(IncompatibleSimulations, match="symbolic"):
        conservative_convex_remap(alias, b, "/fields/scalar/T")
    with h5py.File(a, "r+") as h5:
        h5["points"][0, 2] = 0.2
    with pytest.raises(IncompatibleSimulations, match="nonplanar|concave|volume|hull"):
        conservative_convex_remap(a, b, "/fields/scalar/T")


def test_volume_partition_returns_actual_orig_ids_and_values(tmp_path):
    path = make_cartesian_case(
        tmp_path / "mesh.h5",
        [(0, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 1)],
        [10, 20, 30, 40],
    )
    part = read_volume_partition(path, 1, count=2, field="/fields/scalar/T")
    assert part.original_cell_ids.tolist() == [1, 2]
    assert part.original_point_ids.tolist() == list(range(8, 24))
    np.testing.assert_array_equal(part.field_values, [20, 30])
    assert part.fraction == 0.5
    assert len(part.points) == 16
    assert part.polyhedron_stream[0] == 6
    assert len(part.polyhedron_stream) == 62
    assert int(part.polyhedron_stream.max()) <= 16  # stream includes face counts
    assert part.total_cells == 4


def test_volume_partition_detects_invalids_and_budgets(tmp_path):
    path = make_cartesian_case(tmp_path / "mesh.h5", [(0, 1)], [25])
    with pytest.raises(ValueError, match="range"):
        read_volume_partition(path, 1, count=1)
    with pytest.raises(ValueError, match="count|budget"):
        read_volume_partition(path, 0, count=513)
    with pytest.raises(ValueError, match="scalar"):
        with h5py.File(path, "r+") as h5:
            h5["fields/scalar/T"][0] = np.nan
        read_volume_partition(path, 0, count=1, field="/fields/scalar/T")
    with h5py.File(path, "r+") as h5:
        h5["cell_faces"][0] = 9999
    with pytest.raises(ValueError, match="faces"):
        read_volume_partition(path, 0, count=1)


def test_volume_partition_refuses_symlink_and_export_cli(tmp_path, capsys):
    from cfdx.volume_partitions import main

    path = make_cartesian_case(tmp_path / "mesh.h5", [(0, 1)], [25])
    alias = tmp_path / "alias.h5"
    alias.symlink_to(path)
    with pytest.raises(ValueError, match="symbolic"):
        read_volume_partition(alias, 0, count=1)
    assert main([str(path), "--start", "0", "--count", "1"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["cells_in_tile"] == 1
    assert payload["total_cells"] == 1
    assert payload["scientific_convergence_certificate"] is False
    assert "not_device_certified" in payload["renderer"]


@pytest.mark.skipif(
    __import__("importlib").util.find_spec("pyvista") is None, reason="PyVista optional"
)
def test_volume_partition_reconstructs_native_vtk_polyhedron(tmp_path):
    import pyvista as pv

    source = make_cartesian_case(tmp_path / "cube.h5", [(0, 1)], [5.0])
    piece = read_volume_partition(source, 0, count=1, field="/fields/scalar/T")
    grid = piece.as_pyvista(pv)
    assert grid.n_cells == 1
    assert grid.n_points == 8
    assert grid.cell_data["CFDX_original_cell_id"].tolist() == [0]
    assert grid.cell_data["/fields/scalar/T"].tolist() == pytest.approx([5.0])
    assert grid.volume == pytest.approx(1.0)


@pytest.mark.skipif(
    __import__("importlib").util.find_spec("PySide6") is None, reason="Qt optional"
)
def test_volume_tile_editor_emits_real_request_without_starting_solver():
    from cfdx.industrial_tools_panel import IndustrialToolsPanel
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication(["cfdx-volume-tile-test"])
    panel = IndustrialToolsPanel()
    sent = []
    panel.volume_partition_requested.connect(lambda *args: sent.append(args))
    panel.volume_file.setText("/tmp/test.cfdx.h5")
    panel.volume_start.setValue(12)
    panel.volume_count.setValue(48)
    panel.volume_scalar.setText("/fields/scalar/T")
    panel.request_volume_partition()
    assert sent == [("/tmp/test.cfdx.h5", 12, 48, "/fields/scalar/T")]
    panel.close()
    app.processEvents()


def test_volume_partition_blocks_invalid_nonmanifold_edges(tmp_path):
    path = make_cartesian_case(tmp_path / "badpoly.h5", [(0, 1)], [300])
    with h5py.File(path, "r+") as h5:
        # Keep cell face IDs unique but replace a face with a different face's
        # vertex data, breaking edge incidence and producing an open cell.
        h5["face_vertices"][20:24] = h5["face_vertices"][0:4]
    with pytest.raises(ValueError, match="nonmanifold|open"):
        read_volume_partition(path, 0, count=1)


@pytest.mark.skipif(
    __import__("importlib").util.find_spec("vtk") is None, reason="VTK Python optional"
)
def test_actual_polyhedron_tile_is_accepted_by_native_vtk(tmp_path):
    """No renderer or GPU required: verify actual VTK polyhedron cell topology."""
    import vtk

    path = make_cartesian_case(tmp_path / "vtk_polyhedron.h5", [(0, 1)], [300.0])
    tile = read_volume_partition(path, 0, count=1)
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(len(tile.points))
    for index, point in enumerate(tile.points):
        points.SetPoint(index, *map(float, point))
    face_array = vtk.vtkCellArray()
    stream = tile.polyhedron_stream.tolist()
    cursor = 1
    for _ in range(stream[0]):
        length = stream[cursor]
        vertices = stream[cursor + 1 : cursor + 1 + length]
        cursor += length + 1
        face_array.InsertNextCell(length)
        for vertex in vertices:
            face_array.InsertCellPoint(vertex)
    assert cursor == len(stream)
    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)
    grid.Allocate(1)
    grid.InsertNextCell(
        vtk.VTK_POLYHEDRON, len(tile.points), list(range(len(tile.points))), face_array
    )
    assert grid.GetNumberOfCells() == 1
    cell = grid.GetCell(0)
    assert cell.GetCellType() == vtk.VTK_POLYHEDRON
    assert cell.GetNumberOfEdges() == 12
    assert cell.GetBounds() == pytest.approx((0, 1, 0, 1, 0, 1))

"""LOT9 proof contracts: conservative overlap, bounded LOD and native timings."""

from __future__ import annotations

import json

import h5py
import numpy as np
import pytest
from cfdx.conservative_remap import conservative_box_remap
from cfdx.simulation_compare import IncompatibleSimulations, compare_run_reports
from cfdx.surface_streaming import read_surface_preview


def make_cartesian_case(path, intervals, values, *, offset=0.0, time=1.0, units="K"):
    """Write independently closed Cartesian hexes with six real quad faces."""
    points = []
    vertices = []
    fo = [0]
    cell_faces = []
    co = [0]
    owners = []
    neighbours = []
    for cell, (lo, hi) in enumerate(intervals):
        orig = len(points)
        points.extend(
            [
                [lo, 0, 0],
                [hi, 0, 0],
                [hi, 1, 0],
                [lo, 1, 0],
                [lo, 0, 1],
                [hi, 0, 1],
                [hi, 1, 1],
                [lo, 1, 1],
            ]
        )
        for face in (
            (0, 3, 2, 1),
            (4, 5, 6, 7),
            (0, 1, 5, 4),
            (3, 7, 6, 2),
            (0, 4, 7, 3),
            (1, 2, 6, 5),
        ):
            cell_faces.append(len(fo) - 1)
            vertices.extend([orig + i for i in face])
            fo.append(len(vertices))
            owners.append(cell)
            neighbours.append(-1)
        co.append(len(cell_faces))
    with h5py.File(path, "w") as h5:
        for key, val in {
            "points": points,
            "face_vertices": vertices,
            "face_offsets": fo,
            "cell_faces": cell_faces,
            "cell_offsets": co,
            "owner": owners,
            "neighbour": neighbours,
        }.items():
            h5.create_dataset(
                key, data=np.asarray(val, dtype="f8" if key == "points" else "i8")
            )
        h5.create_dataset("patch_face_ids", data=list(range(len(owners))))
        h5.create_dataset("patch_face_offsets", data=[0, len(owners)])
        h5.attrs["boundary_patches"] = "walls:WALL:0:6"
        h5.attrs["physical_time"] = time
        field = h5.create_dataset(
            "fields/scalar/T", data=np.asarray(values, dtype=float) + offset
        )
        field.attrs["units"] = units
    return path


def test_conservative_remap_coarsen_and_refine(tmp_path):
    a = make_cartesian_case(
        tmp_path / "fine.h5",
        [(0, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 1)],
        [1, 3, 5, 7],
    )
    b = make_cartesian_case(tmp_path / "coarse.h5", [(0, 0.5), (0.5, 1)], [2, 6])
    mapped, report = conservative_box_remap(a, b, "/fields/scalar/T")
    np.testing.assert_allclose(mapped, [2, 6], atol=1e-14)
    assert report.relative_conservation_error < 1e-13
    assert report.linf_difference < 1e-13
    assert report.source_integral == pytest.approx(4)
    assert not report.certificate
    mapped_reverse, back = conservative_box_remap(b, a, "/fields/scalar/T")
    np.testing.assert_allclose(mapped_reverse, [2, 2, 6, 6])
    assert back.mapped_integral == pytest.approx(back.source_integral)


def test_conservative_remap_rejects_partial_domain_and_skew(tmp_path):
    a = make_cartesian_case(tmp_path / "a.h5", [(0, 0.5), (0.5, 1)], [0, 1])
    b = make_cartesian_case(tmp_path / "b.h5", [(0, 0.5), (0.6, 1)], [0, 1])
    with pytest.raises(IncompatibleSimulations, match="volumes|fully covered"):
        conservative_box_remap(a, b, "/fields/scalar/T")
    b = make_cartesian_case(tmp_path / "c.h5", [(0, 0.5), (0.5, 1)], [0, 1])
    with h5py.File(b, "r+") as h:
        h["points"][0, 1] = 0.3
    with pytest.raises(IncompatibleSimulations, match="non-Cartesian|corner|face"):
        conservative_box_remap(a, b, "/fields/scalar/T")
    b = make_cartesian_case(tmp_path / "d.h5", [(0, 0.5), (0.5, 1)], [0, 1], time=2.0)
    with pytest.raises(IncompatibleSimulations, match="physical times"):
        conservative_box_remap(a, b, "/fields/scalar/T")


def test_surface_sample_bounded_source_indices_and_coverage(tmp_path):
    intervals = [(i / 24, (i + 1) / 24) for i in range(24)]
    p = make_cartesian_case(tmp_path / "mesh.h5", intervals, [1] * 24)
    preview = read_surface_preview(p, max_faces=12, max_vertices=200)
    assert preview.n_boundary_faces == 144
    assert 0 < preview.coverage < 0.1
    assert len(preview.original_face_ids) == 12
    assert preview.points.shape[1] == 3
    assert len(preview.faces) > 12
    assert int(preview.faces.max()) < len(preview.points)
    with pytest.raises(ValueError, match="budget"):
        read_surface_preview(p, max_faces=0)
    alias = tmp_path / "alias.h5"
    alias.symlink_to(p)
    with pytest.raises(ValueError, match="symbolic"):
        read_surface_preview(alias)


def test_surface_preview_corrupted_topology_is_not_rendered(tmp_path):
    p = make_cartesian_case(tmp_path / "mesh.h5", [(0, 1)], [2])
    with h5py.File(p, "r+") as h:
        h["patch_face_ids"][0] = 1000
    with pytest.raises(ValueError, match="invalid face"):
        read_surface_preview(p)


def _write_run(folder, seconds):
    folder.mkdir()
    (folder / "execution.json").write_text(
        json.dumps(
            {
                "format": "CFDX-EXECUTION",
                "schema_version": 1,
                "process_exit_code": 1,
                "converged": False,
                "iterations": 1,
                "elapsed_wall_seconds": seconds,
            }
        )
    )
    (folder / "convergence.json").write_text(
        json.dumps(
            {
                "format": "CFDX-CONVERGENCE",
                "schema_version": 1,
                "converged": False,
                "iterations": 1,
                "history": [{"iteration": 1, "mass_balance": 0.1}],
            }
        )
    )
    return folder


def test_wall_time_from_native_artifacts_and_bad_values(tmp_path):
    a = _write_run(tmp_path / "a", 3.0)
    b = _write_run(tmp_path / "b", 4.75)
    assert compare_run_reports(a, b)["elapsed_wall_seconds"] == {
        "reference": 3.0,
        "candidate": 4.75,
        "delta": 1.75,
    }
    data = json.loads((b / "execution.json").read_text())
    data["elapsed_wall_seconds"] = -2
    (b / "execution.json").write_text(json.dumps(data))
    with pytest.raises(IncompatibleSimulations, match="elapsed_wall_seconds"):
        compare_run_reports(a, b)
    del data["elapsed_wall_seconds"]
    (b / "execution.json").write_text(json.dumps(data))
    assert (
        compare_run_reports(a, b)["elapsed_wall_seconds"]
        == "not_recorded_in_native_reports"
    )


def test_surface_preview_scales_with_sample_budget_not_boundary_face_count(tmp_path):
    p = make_cartesian_case(tmp_path / "huge_boundary.h5", [(0, 1)], [2])
    with h5py.File(p, "r+") as h:
        del h["patch_face_ids"]
        del h["patch_face_offsets"]
        n = 1_000_000
        # Synthetic I/O stress only; repeated face references are NOT an
        # actual million-face mesh and are not physical CFD qualification.
        h.create_dataset(
            "patch_face_ids",
            shape=(n,),
            dtype="i8",
            chunks=(4096,),
            compression="gzip",
            fillvalue=0,
        )
        h.create_dataset("patch_face_offsets", data=[0, n])
    sampled = read_surface_preview(p, max_faces=128, max_vertices=512)
    assert len(sampled.original_face_ids) == 128
    assert sampled.n_boundary_faces == 1_000_000
    assert sampled.coverage == pytest.approx(128 / 1_000_000)
    assert sampled.points.shape == (4, 3)


def test_headless_diagnostic_cli_reports_no_physical_certification(tmp_path, capsys):
    from cfdx.conservative_remap import main as remap_cli
    from cfdx.surface_streaming import main as surface_cli

    a = make_cartesian_case(tmp_path / "a.h5", [(0, 1)], [315])
    b = make_cartesian_case(tmp_path / "b.h5", [(0, 0.5), (0.5, 1)], [315, 315])
    assert remap_cli([str(a), str(b), "/fields/scalar/T"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data["source_integral"] == pytest.approx(315)
    assert data["relative_conservation_error"] < 1e-12
    assert not data["certificate"]
    assert surface_cli([str(b), "--max-faces", "3"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data["total_boundary_faces"] == 12
    assert data["sampled_boundary_faces"] == 3
    assert not data["full_volume_available"]


def test_conservative_remap_refuses_symlink_input(tmp_path):
    a = make_cartesian_case(tmp_path / "a.h5", [(0, 1)], [1])
    b = make_cartesian_case(tmp_path / "b.h5", [(0, 1)], [1])
    link = tmp_path / "link.h5"
    link.symlink_to(a)
    with pytest.raises(IncompatibleSimulations, match="symbolic"):
        conservative_box_remap(link, b, "/fields/scalar/T")

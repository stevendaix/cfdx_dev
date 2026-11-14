"""Geometric/IO regression for the published OpenFOAM 2D BFS topology.

These tests verify mesh integrity and source BC mapping, not physical validity.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from cfdx.io.adapters.openfoam import OpenFOAMAdapter
from cfdx.io.schema import BCType

from scripts import generate_th03_openfoam_blockmesh as b


def _fake_ascii_polymesh(root: Path, mesh: dict) -> None:
    """Write a fixture with actual OpenFOAM ASCII list syntax and face ordering."""
    poly = root / "constant" / "polyMesh"
    poly.mkdir(parents=True)
    n_faces = mesh["n_faces"]
    internal = np.flatnonzero(mesh["neighbour"] >= 0).tolist()
    by_patch = [
        (patch["name"], patch["type"], patch["face_ids"])
        for patch in mesh["boundary_patch_assignment"]
    ]
    ordered = internal + [fid for _, _, faces in by_patch for fid in faces]
    assert len(ordered) == n_faces
    assert len(set(ordered)) == n_faces
    header = "FoamFile\n{ version 2.0; format ascii; class vectorField; }\n"

    def write_list(name: str, values: list[str]) -> None:
        (poly / name).write_text(
            header + f"{len(values)}\n(\n" + "\n".join(values) + "\n)\n"
        )

    write_list(
        "points",
        ["(" + " ".join(f"{x:.16g}" for x in pt) + ")" for pt in mesh["points"]],
    )

    def face(fid: int) -> str:
        first, last = mesh["face_offsets"][fid : fid + 2]
        ids = mesh["face_vertices"][first:last]
        return str(len(ids)) + "(" + " ".join(map(str, ids)) + ")"

    write_list("faces", [face(fid) for fid in ordered])
    write_list("owner", [str(mesh["owner"][fid]) for fid in ordered])
    write_list("neighbour", [str(mesh["neighbour"][fid]) for fid in internal])
    patches = []
    start = len(internal)
    for name, oftype, fids in by_patch:
        if name in ("inlet", "outlet"):
            oftype = "patch"
        patches.append(
            f"{name} {{ type {oftype}; nFaces {len(fids)}; startFace {start}; }}"
        )
        start += len(fids)
    (poly / "boundary").write_text(
        "FoamFile\n{ format ascii; }\n"
        + f"{len(patches)}\n(\n"
        + "\n".join(patches)
        + "\n)\n"
    )


def test_six_block_geometry_has_matching_interfaces():
    mesh = b.make_structured_mesh("coarse")
    assert mesh["n_cells"] == sum(
        block["nx"] * block["ny"] for block in b.block_specs("coarse")
    )
    assert mesh["n_cells"] == 1280
    assert mesh["n_faces"] > mesh["n_cells"]
    assert mesh["minimum_cell_area_m2"] > 0
    np.testing.assert_allclose(
        np.min(mesh["points"], axis=0), [-130 * b.H, 0, 0], atol=1.0e-13
    )
    np.testing.assert_allclose(
        np.max(mesh["points"], axis=0), [50 * b.H, 9 * b.H, b.H], atol=1.0e-13
    )
    assert int(np.count_nonzero(mesh["neighbour"] < 0)) == sum(
        len(p["face_ids"]) for p in mesh["boundary_patch_assignment"]
    )
    seen = set()
    for patch in mesh["boundary_patch_assignment"]:
        for fid in patch["face_ids"]:
            assert mesh["neighbour"][fid] == -1
            assert fid not in seen
            seen.add(fid)
    assert any(
        p["name"] == "lowerWallStartup" and len(p["face_ids"]) > 0
        for p in mesh["boundary_patch_assignment"]
    )
    x_h = np.unique(np.round(mesh["points"][:, 0] / b.H, 8))
    y_h = np.unique(np.round(mesh["points"][:, 1] / b.H, 8))
    assert all(np.any(np.isclose(x_h, v)) for v in (-130, -110, 0, 8, 50))
    assert all(np.any(np.isclose(y_h, v)) for v in (0, 1, 3.125, 9))


def test_full_six_block_mesh_counts_and_stretch():
    mesh = b.make_structured_mesh("full")
    assert mesh["n_cells"] == 20540
    assert len(mesh["boundary_patch_assignment"]) == 8
    assert all(len(p["face_ids"]) > 0 for p in mesh["boundary_patch_assignment"])
    assert mesh["minimum_cell_area_m2"] > 0


def test_real_ascii_poly_mesh_roundtrip_and_patch_semantics(tmp_path: Path):
    original = b.make_structured_mesh("smoke")
    _fake_ascii_polymesh(tmp_path, original)
    read = b.load_openfoam_polymesh(tmp_path)
    assert read["n_cells"] == original["n_cells"]
    assert read["n_faces"] == original["n_faces"]
    np.testing.assert_allclose(
        read["points"], original["points"], rtol=1.0e-14, atol=1.0e-14
    )
    assert {p["name"]: p["type"] for p in read["boundary_patch_assignment"]} == {
        p["name"]: (p["name"] if p["name"] in ("inlet", "outlet") else p["type"])
        for p in original["boundary_patch_assignment"]
    }
    output = tmp_path / "imported.h5"
    summary = b.write_case(
        output, read, mesh_source="fixture ASCII polyMesh", max_iterations=2
    )
    assert summary["cells"] == original["n_cells"]
    with h5py.File(output) as h5:
        setup = json.loads(h5.attrs["case_setup_json"])
        assert setup["turbulence_model"] == "k_omega_sst"
        # Verify *physical* outlet semantics, not just OpenFOAM patch names.
        # OUTLET_PRESSURE requires PRESSURE_OUTLET in the production adapter.
        outlet = next(
            p for p in setup["boundary_conditions"] if p["patch_name"] == "outlet"
        )
        assert outlet["type"] == BCType.PRESSURE_OUTLET.value
        assert outlet["value_type"] == "outlet_pressure"
        assert outlet["pressure"] == 0.0
        assert setup["materials"][0]["dynamic_viscosity"] == pytest.approx(b.RHO * b.NU)
        assert {p["patch_name"] for p in setup["boundary_conditions"]} == set(
            b.PATCH_TYPES
        )
        assert "lowerWallStartup:" in h5.attrs["boundary_patches"]
        assert "lowerWall:" in h5.attrs["boundary_patches"]
        assert len(h5["owner"]) == original["n_faces"]
        assert int(np.count_nonzero(h5["neighbour"][:] < 0)) == sum(
            summary["patch_faces"].values()
        )


def test_generic_openfoam_patch_is_not_mapped_to_a_wall():
    adapter = OpenFOAMAdapter()
    assert adapter._map_patch_type("patch") == BCType.UNKNOWN
    assert adapter._map_patch_type("symmetryPlane") == BCType.SYMMETRY
    assert adapter._map_patch_type("pressureOutlet") == BCType.PRESSURE_OUTLET


def test_import_rejects_invalid_patch_and_binary(tmp_path: Path):
    mesh = b.make_structured_mesh("smoke")
    _fake_ascii_polymesh(tmp_path, mesh)
    boundary = tmp_path / "constant/polyMesh/boundary"
    boundary.write_text(
        boundary.read_text().replace("lowerWallStartup", "incorrectPatch")
    )
    with pytest.raises(ValueError, match="unexpected OpenFOAM patch"):
        b.load_openfoam_polymesh(tmp_path)
    boundary.write_text(
        boundary.read_text().replace("incorrectPatch", "lowerWallStartup")
    )
    points = tmp_path / "constant/polyMesh/points"
    points.write_text(points.read_text().replace("format ascii;", "format binary;"))
    with pytest.raises(ValueError, match="mesh parse failed"):
        b.load_openfoam_polymesh(tmp_path)


def test_independent_mesh_quality_check_passes_real_geometry():
    from scripts import audit_th03_openfoam_mesh as inspect

    result = inspect.audit(b.make_structured_mesh("coarse"))
    assert result["status"] == "PASS_GEOMETRY_ONLY"
    assert result["qualification"] == "NOT_QUALIFIED"
    assert result["flipped_internal_faces"] == 0
    assert result["flipped_boundary_faces"] == 0
    assert result["max_nonorthogonality_deg"] < 90


def test_quality_check_fails_closed_on_inverted_internal_face():
    from scripts import audit_th03_openfoam_mesh as inspect

    mesh = b.make_structured_mesh("smoke")
    f = int(np.flatnonzero(mesh["neighbour"] >= 0)[10])
    a = int(mesh["face_offsets"][f])
    z = int(mesh["face_offsets"][f + 1])
    mesh["face_vertices"][a:z] = mesh["face_vertices"][a:z][::-1]
    result = inspect.audit(mesh)
    assert result["status"] == "FAIL_GEOMETRY"
    assert result["flipped_internal_faces"] >= 1


def test_empty_patches_are_only_extrusion_caps():
    from scripts import audit_th03_openfoam_mesh as inspect

    for resolution in ("smoke", "coarse"):
        result = inspect.audit(b.make_structured_mesh(resolution))
        assert result["status"] == "PASS_GEOMETRY_ONLY"
        assert result["invalid_extrusion_patch_faces"] == 0


def test_mesh_audit_rejects_swapped_empty_and_wall_labels():
    from scripts import audit_th03_openfoam_mesh as inspect

    mesh = b.make_structured_mesh("smoke")
    patches = {p["name"]: p["face_ids"] for p in mesh["boundary_patch_assignment"]}
    # Keep patch coverage perfect: an index/duplicate-only audit would PASS.
    patches["front"][0], patches["lowerWall"][0] = (
        patches["lowerWall"][0],
        patches["front"][0],
    )
    result = inspect.audit(mesh)
    assert result["status"] == "FAIL_GEOMETRY"
    assert result["invalid_extrusion_patch_faces"] == 2
    with pytest.raises(ValueError, match="EMPTY/wall patch geometry"):
        b.write_case(
            Path("never_created.h5"), mesh, mesh_source="intentionally invalid fixture"
        )
    assert any("incompatible extrusion patch" in issue for issue in result["issues"])


def test_krylov_config_is_exported_and_rejects_invalid_entries(tmp_path):
    with pytest.raises(ValueError, match="invalid pressure Krylov controls"):
        b.make_case_setup(5, 0.25, 0.1, pressure_restart=0)
    with pytest.raises(ValueError, match="invalid pressure Krylov controls"):
        b.make_case_setup(5, 0.25, 0.1, pressure_scaling="bad")
    mesh = b.make_structured_mesh("smoke")
    dest = tmp_path / "case.h5"
    b.write_case(
        dest,
        mesh,
        mesh_source="test",
        pressure_restart=192,
        pressure_scaling="symmetric_diagonal",
    )
    with h5py.File(dest) as h:
        settings = json.loads(h.attrs["case_setup_json"])["numerics"]["raw_settings"]
    assert settings["pressure_gmres_restart"] == "192"
    assert settings["pressure_linear_scaling"] == "symmetric_diagonal"


def test_cli_requires_native_case_unless_reconstruction_is_explicit(tmp_path, monkeypatch):
    output = tmp_path / "diagnostic.cfdx.h5"
    monkeypatch.setattr(
        sys, "argv", ["generate_th03_openfoam_blockmesh.py", "--output", str(output)]
    )
    with pytest.raises(SystemExit) as error:
        b.main()
    assert error.value.code == 2
    assert not output.exists()

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "generate_th03_openfoam_blockmesh.py",
            "--output",
            str(output),
            "--resolution",
            "smoke",
            "--allow-reconstructed-mesh",
            "--max-iterations",
            "1",
        ],
    )
    assert b.main() == 0
    with h5py.File(output) as case:
        assert case.attrs["source_format"] == "openfoam-blockmesh-reconstruction"
        setup = json.loads(case.attrs["case_setup_json"])
        assert setup["source_metadata"]["physical_qualification"] == "NOT_QUALIFIED"

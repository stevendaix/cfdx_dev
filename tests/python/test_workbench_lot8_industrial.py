"""LOT8: contracts for full cell map, native BC, streaming, comparison and CI."""

from __future__ import annotations

import hashlib
import json
import struct

import h5py
import numpy as np
import pytest
from cfdx.cell_field_map import CellFieldMap, CellMapError
from cfdx.native_boundary_bridge import (
    NativeBoundaryError,
    build_native_boundary_payload,
    write_native_boundary_case,
)
from cfdx.progressive_fields import field_catalog, stream_field_statistics
from cfdx.qualification_center import STATUSES
from cfdx.qualification_evidence import qualification_evidence_manifest
from cfdx.simulation_compare import (
    IncompatibleSimulations,
    compare_run_reports,
    compare_scalar_fields,
)


def _binary_fixture(path, complete=True, malformed=False):
    data = (
        b"CFDXMAP1" + struct.pack("<QQId", 3, 6, len(b"/mesh.h5"), 1e-6) + b"/mesh.h5"
    )
    data += b"FRME" + struct.pack("<QB", 4, 1)
    data += (
        struct.pack("<ddHB", 1e-8, 0.0, 0, 1)
        + struct.pack("<ddHB", 1e-3, 0.25, 3, 1)
        + struct.pack("<ddHB", np.nan, np.nan, 0, 0)
    )
    data += b"FAIL" if malformed else b"DONE"
    if complete:
        data += b"CFND" + struct.pack("<Q", 4)
    path.write_bytes(data)
    return path


def test_full_cell_map_stable_scale_and_unknown_cells(tmp_path):
    with CellFieldMap(_binary_fixture(tmp_path / "map.bin")) as m:
        assert m.complete and m.n_cells == 3 and m.frame().iteration == 4
        scalar, clim, _title = m.scalar("momentum")
        assert clim == (-2, 4)
        assert scalar[0] == pytest.approx(-2)
        assert scalar[1] == pytest.approx(3)
        assert np.isnan(scalar[2])
        assert m.summary()["above_reference"] == 1
        assert not m.summary()["convergence_certificate"]


def test_full_cell_map_bad_marker_and_symlink(tmp_path):
    source = _binary_fixture(tmp_path / "bad.bin", malformed=True)
    with pytest.raises(CellMapError, match="frame"):
        CellFieldMap(source)
    source = _binary_fixture(tmp_path / "valid.bin")
    alias = tmp_path / "alias.bin"
    alias.symlink_to(source)
    with pytest.raises(CellMapError, match="safely"):
        CellFieldMap(alias)


def _h5(path, values=(1.0, 2.0, 3.0), *, offset=0.0, units="Pa"):
    with h5py.File(path, "w") as h:
        for name, val in {
            "points": [[0, 0, 0], [1, 0, 0]],
            "face_vertices": [0, 1],
            "face_offsets": [0, 2],
            "owner": [0],
            "neighbour": [-1],
            "cell_faces": [0, 0, 0],
            "cell_offsets": [0, 1, 2, 3],
        }.items():
            h.create_dataset(name, data=val)
        d = h.create_dataset("fields/scalar/p", data=np.asarray(values) + offset)
        d.attrs["units"] = units
        h.attrs["physical_time"] = 2.0
    return path


def test_streamed_stats_and_comparison(tmp_path):
    a, b = _h5(tmp_path / "a.h5"), _h5(tmp_path / "b.h5", offset=2)
    assert field_catalog(a) == ("/fields/scalar/p",)
    stats = stream_field_statistics(a, "/fields/scalar/p", max_block_bytes=1024)
    assert stats.count == 3 and stats.mean == 2 and stats.units == "Pa"
    result = compare_scalar_fields(a, b, "/fields/scalar/p", block_bytes=1024)
    assert result.linf_difference == 2
    assert result.rms_difference == 2
    assert not result.certificate
    with h5py.File(b, "r+") as h:
        h["points"][0, 0] = 3
    with pytest.raises(IncompatibleSimulations, match="topology|geometry"):
        compare_scalar_fields(a, b, "/fields/scalar/p")


def test_streaming_rejects_untrusted_fields_and_units(tmp_path):
    a, b = _h5(tmp_path / "a.h5"), _h5(tmp_path / "b.h5", units="bar")
    with pytest.raises(IncompatibleSimulations, match="units"):
        compare_scalar_fields(a, b, "/fields/scalar/p")
    with h5py.File(a, "r+") as h:
        h["fields/scalar/alias"] = h5py.SoftLink("/fields/scalar/p")
    assert field_catalog(a) == ("/fields/scalar/p",)


def test_native_bridge_preserves_case_setup_and_multi_field(tmp_path):
    patch = {
        "patch_name": "inlet",
        "type": "inlet",
        "value_type": "fixed",
        "velocity_vector": [5.0, 0, 0],
        "pressure": 0.0,
        "sst_scalar_conditions": [],
        "raw_params": {"unchanged": "yes"},
    }
    baseline = {
        "boundary_conditions": [patch],
        "physics_model": "incompressible",
        "numerics": {"coupled_solver": "SIMPLE"},
        "materials": [{"density": 1.2}],
    }
    editor = {
        "inlet": {
            "type": "inlet",
            "fields": ["velocity", "pressure", "temperature", "k", "omega"],
            "field_conditions": {
                "velocity": {
                    "velocity_type": "FIXED_VALUE",
                    "velocity_value": [1, 2, 3],
                },
                "pressure": {"scalar_type": "ZERO_GRADIENT", "value": 0, "gradient": 0},
                "temperature": {
                    "scalar_type": "FIXED_VALUE",
                    "value": 300,
                    "gradient": 0,
                },
                "k": {"scalar_type": "FIXED_VALUE", "value": 0.001, "gradient": 0},
                "omega": {"scalar_type": "FIXED_VALUE", "value": 10, "gradient": 0},
            },
        }
    }
    data = build_native_boundary_payload(baseline, editor)
    bc = data["boundary_conditions"][0]
    assert data["numerics"] == baseline["numerics"]
    assert bc["field_conditions"]["velocity"]["value"] == [1, 2, 3]
    assert bc["field_conditions"]["temperature"]["value"] == [300]
    assert {c["field"] for c in bc["sst_scalar_conditions"]} == {"k", "omega"}
    src = tmp_path / "input.cfdx.h5"
    with h5py.File(src, "w") as h:
        h.attrs["case_setup_json"] = json.dumps(baseline)
    dest = tmp_path / "output.cfdx.h5"
    write_native_boundary_case(src, dest, editor)
    assert dest.exists() and src.exists()
    with h5py.File(src, "r") as h:
        assert json.loads(h.attrs["case_setup_json"]) == baseline
    with h5py.File(dest, "r") as h:
        assert (
            json.loads(h.attrs["case_setup_json"])["boundary_conditions"][0][
                "field_conditions"
            ]["pressure"]["type"]
            == "ZERO_GRADIENT"
        )
    with pytest.raises(NativeBoundaryError, match="destination exists"):
        write_native_boundary_case(src, dest, editor)


def test_native_bridge_rejects_unsafe_turbulence(tmp_path):
    original = {
        "boundary_conditions": [
            {"patch_name": "wall", "type": "wall", "sst_scalar_conditions": []}
        ]
    }
    editor = {
        "wall": {
            "type": "wall",
            "fields": ["omega"],
            "scalar_type": "ZERO_GRADIENT",
            "value": 0,
        }
    }
    with pytest.raises(NativeBoundaryError, match="wall omega"):
        build_native_boundary_payload(original, editor)


def test_ci_manifest_does_not_promote_registry_pass_without_proof(tmp_path):
    reg = tmp_path / "reg.json"
    reg.write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "campaign": "T",
                "status_policy": {s: s for s in STATUSES},
                "cases": [
                    {
                        "id": "SEP-BFS",
                        "domain": "Separation",
                        "case": "BFS",
                        "status": "PASS",
                    }
                ],
            }
        )
    )
    root = tmp_path / "proofs"
    root.mkdir()
    sha = "a" * 40
    before = qualification_evidence_manifest(reg, root, sha)
    assert before["cases"][0]["proof_status"] == "UNVERIFIED"
    evidence = root / "bfs.json"
    evidence.write_text("physical measurements")
    proof = {
        "case_id": "SEP-BFS",
        "commit_sha": sha,
        "status": "PASS",
        "qualifier": "physical",
        "artifacts": ["bfs.json"],
        "artifact_sha256": {
            "bfs.json": hashlib.sha256(evidence.read_bytes()).hexdigest()
        },
        "checks": [
            {
                "name": "reattachment relative error",
                "error": 0.02,
                "limit": 0.05,
                "passed": True,
            }
        ],
    }
    (root / "SEP-BFS.json").write_text(json.dumps(proof))
    after = qualification_evidence_manifest(reg, root, sha)
    assert after["physically_evidenced"] == 1
    assert (
        qualification_evidence_manifest(reg, root, "b" * 40)["physically_evidenced"]
        == 0
    )
    proof["checks"][0]["error"] = 0.06
    (root / "SEP-BFS.json").write_text(json.dumps(proof))
    assert qualification_evidence_manifest(reg, root, sha)["physically_evidenced"] == 0


def test_comparison_rejects_mismatched_number_of_cells(tmp_path):
    a, b = _h5(tmp_path / "a.h5"), _h5(tmp_path / "b.h5")
    with h5py.File(b, "r+") as h:
        del h["fields/scalar/p"]
        h.create_dataset("fields/scalar/p", data=[1.0, 2.0])
        h["fields/scalar/p"].attrs["units"] = "Pa"
    with pytest.raises(IncompatibleSimulations, match="cell field size"):
        compare_scalar_fields(a, b, "/fields/scalar/p")


def test_comparison_rejects_nonfinite_physical_time(tmp_path):
    a, b = _h5(tmp_path / "a.h5"), _h5(tmp_path / "b.h5")
    with h5py.File(b, "r+") as h:
        h.attrs["physical_time"] = float("nan")
    with pytest.raises(IncompatibleSimulations, match="physical time"):
        compare_scalar_fields(a, b, "/fields/scalar/p")


def _run_fixture(path, residual=0.2, force=1.0):
    path.mkdir()
    (path / "execution.json").write_text(
        json.dumps(
            {
                "format": "CFDX-EXECUTION",
                "schema_version": 1,
                "converged": False,
                "iterations": 4,
                "process_exit_code": 1,
            }
        )
    )
    (path / "convergence.json").write_text(
        json.dumps(
            {
                "format": "CFDX-CONVERGENCE",
                "schema_version": 1,
                "iterations": 4,
                "converged": False,
                "history": [
                    {
                        "iteration": 4,
                        "momentum_equation_residual_rhs_relative": residual,
                        "continuity_normalized": 0.001,
                    }
                ],
            }
        )
    )
    (path / "wall_forces.json").write_text(
        json.dumps(
            {
                "format": "CFDX-WALL-FORCES",
                "schema_version": 1,
                "method": "quadratic_wls_3d",
                "normal_convention": "body_to_fluid",
                "pressure_reference_pa": 0,
                "moment_center_m": [0, 0, 0],
                "force_definition": "fluid_on_body = -p*n_body + tau*n_body",
                "wall_patches": ["wall"],
                "pressure_force_n": [force, 0, 0],
                "viscous_force_n": [0, 0, 0],
                "total_force_n": [force, 0, 0],
                "total_moment_nm": [0, 0, 0],
            }
        )
    )
    return path


def test_run_telemetry_and_dimensional_forces_comparison_without_false_vv(tmp_path):
    a = _run_fixture(tmp_path / "a", residual=0.2, force=1.0)
    b = _run_fixture(tmp_path / "b", residual=0.1, force=2.0)
    report = compare_run_reports(a, b)
    assert report["final_solver_metrics"]["momentum_equation_residual_rhs_relative"][
        "delta"
    ] == pytest.approx(-0.1)
    assert report["wall_forces"]["vectors"]["total_force_n"]["delta"] == [1.0, 0.0, 0.0]
    assert report["comparison_certificate"] is False
    assert report["elapsed_wall_seconds"] == "not_recorded_in_native_reports"
    with open(b / "wall_forces.json", "r+") as stream:
        data = json.loads(stream.read())
        stream.seek(0)
        stream.truncate()
        data["normal_convention"] = "wrong"
        stream.write(json.dumps(data))
    with pytest.raises(IncompatibleSimulations, match="incompatible"):
        compare_run_reports(a, b)


def test_progressive_catalog_rejects_symlinked_group(tmp_path):
    case = tmp_path / "case.h5"
    with h5py.File(case, "w") as h:
        h.create_dataset("real/p", data=[1.0])
        h.create_group("fields")
        h["fields/scalar"] = h5py.SoftLink("/real")
    with pytest.raises(ValueError, match="hard-linked"):
        field_catalog(case)

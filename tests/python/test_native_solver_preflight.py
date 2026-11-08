"""GUI and native production solver must agree about runnable case schema."""
from __future__ import annotations

import json
from unittest.mock import patch

import h5py
import pytest

from cfdx import CFDXSession
from cfdx.case import Case
from cfdx.case_io import read_case, save_case
from cfdx.native_solver_preflight import (
    check_native_production_case,
    is_native_production_solver,
    native_production_command,
)


def test_detects_only_native_executable():
    assert is_native_production_solver("/opt/bin/cfdx_production_solver")
    assert is_native_production_solver("C:/CFDX/cfdx_production_solver.exe")
    assert not is_native_production_solver("custom.py")


def test_gui_only_case_fails_native_preflight(tmp_path):
    case = CFDXSession()
    case.case.set_boundary("inlet", type="inlet", fields=["velocity"])
    path = save_case(case, tmp_path / "case.cfdx.h5")
    with pytest.raises(ValueError, match="case_setup_json absent"):
        check_native_production_case(path)


def test_native_schema_can_be_loaded_without_gui_data(tmp_path):
    path = tmp_path / "case.cfdx.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["case_setup_json"] = json.dumps(
            {
                "boundary_conditions": [
                    {"patch_name": "in", "type": "INLET", "value_type": "FIXED"}
                ]
            }
        )
    check_native_production_case(path)
    command = native_production_command(
        "cfdx_production_solver", path, output_dir="results", max_iterations=300
    )
    assert command == [
        "cfdx_production_solver",
        "--mesh",
        str(path),
        "--output-dir",
        "results",
        "--iterations",
        "300",
    ]


def test_dual_schema_with_gui_boundaries_rejected(tmp_path):
    case = CFDXSession()
    case.case.set_boundary(
        "in", type="inlet", fields=["velocity"], velocity_type="FIXED_VALUE"
    )
    path = save_case(case, tmp_path / "case.cfdx.h5")
    with h5py.File(path, "a") as h5:
        h5.attrs["case_setup_json"] = json.dumps(
            {"boundary_conditions": [{"patch_name": "in"}]}
        )
    with pytest.raises(ValueError, match="sans synchronisation prouvee"):
        check_native_production_case(path)


def test_native_rejects_empty_bc_invalid_json_iterations(tmp_path):
    path = tmp_path / "case.cfdx.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["case_setup_json"] = "{"
    with pytest.raises(ValueError, match="invalide"):
        check_native_production_case(path)
    with h5py.File(path, "a") as h5:
        h5.attrs["case_setup_json"] = json.dumps({"boundary_conditions": []})
    with pytest.raises(ValueError, match="aucune condition limite"):
        check_native_production_case(path)
    with pytest.raises(ValueError, match="max_iterations"):
        native_production_command(
            "cfdx_production_solver", path, output_dir="x", max_iterations=0
        )


def test_case_run_builds_native_flags_and_mpi(tmp_path):
    path = tmp_path / "case.cfdx.h5"
    with h5py.File(path, "w") as h5:
        h5.attrs["case_setup_json"] = json.dumps(
            {"boundary_conditions": [{"patch_name": "in"}]}
        )
    case = Case()
    case.execution.solver = "cfdx_production_solver"
    case.execution.output_dir = "my-results"
    case.execution.mpi_ranks = 2
    case.numerics["max_iterations"] = 850
    with patch("subprocess.call", return_value=0) as call:
        assert case.run(path) == 0
        assert call.call_args.args[0] == [
            "mpiexec",
            "-n",
            "2",
            "cfdx_production_solver",
            "--mesh",
            str(path),
            "--output-dir",
            "my-results",
            "--iterations",
            "850",
        ]


def test_gui_execution_paths_are_saved_and_loaded(tmp_path):
    session = CFDXSession()
    session.case.execution.output_dir = "custom-results"
    session.case.execution.mesh = "mesh.cfdx.h5"
    path = save_case(session, tmp_path / "roundtrip.cfdx.h5")
    restored = read_case(path)
    assert restored.case.execution.output_dir == "custom-results"
    assert restored.case.execution.mesh == "mesh.cfdx.h5"

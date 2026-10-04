from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest


# This test intentionally imports only the canonical I/O package tree.
def test_production_solver_explicit_cfdx_case_e2e(tmp_path: Path) -> None:
    solver_path = os.environ.get("CFDX_PRODUCTION_SOLVER")
    if not solver_path:
        pytest.skip("CFDX_PRODUCTION_SOLVER is provided by ctest")
    solver = Path(solver_path)
    assert solver.is_file()

    root = Path(__file__).resolve().parents[2]
    source = root / "tests" / "data" / "su2" / "mesh_NACA0012_inv.su2"
    config = root / "tests" / "data" / "su2" / "inv_NACA0012_basic.cfg"
    assert source.is_file()
    assert config.is_file()

    case_dir = tmp_path / "explicit_case"
    case_dir.mkdir()
    source_case = case_dir / source.name
    config_case = case_dir / source.with_suffix(".cfg").name
    source_case.write_bytes(source.read_bytes())
    config_case.write_bytes(config.read_bytes())
    case_path = case_dir / "case.cfdx.h5"

    from cfdx.io.converter import convert

    converted = convert(source_case, output=case_path)
    assert converted.case is not None
    selected = {
        (entry.family, entry.configuration_key)
        for entry in converted.case.numerics.selection.entries
    }
    assert ("gradient", "numerics.gradient.gauss") in selected
    assert ("convection", "numerics.convection.upwind") in selected
    assert ("pressure_velocity", "pressure_velocity.simple") in selected
    assert ("linear_solver", "linear.fgmres") in selected

    output_dir = tmp_path / "explicit_run"
    result = subprocess.run(
        [
            str(solver),
            "--mesh",
            str(case_path),
            "--output-dir",
            str(output_dir),
            "--iterations",
            "20",
        ],
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    diagnostics = result.stdout + result.stderr
    assert result.returncode == 0, diagnostics[-12000:]
    assert "Resolved numerical selections:" in diagnostics
    assert "scheme[gradient]=numerics.gradient.gauss" in diagnostics
    assert "scheme[convection]=numerics.convection.upwind" in diagnostics
    assert "scheme[pressure_velocity]=pressure_velocity.simple" in diagnostics
    assert "scheme[linear_solver]=linear.fgmres" in diagnostics
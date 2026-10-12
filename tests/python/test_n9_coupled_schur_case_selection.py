from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import h5py
import pytest

SCHUR_MODELS = ("block_local", "pcd", "lsc", "bfbt", "simple", "simplec")


def test_n9_coupled_schur_uses_case_selection_and_production_resolution(
    tmp_path: Path,
) -> None:
    """Prove requested -> resolved -> production dispatch through case.cfdx.h5.

    This is deliberately an integration/selection gate, not a physical
    qualification gate. The heavy N9 physical matrix remains responsible for
    convergence, residual, conservation and QoI acceptance.
    """
    solver_path = os.environ.get("CFDX_PRODUCTION_SOLVER")
    if not solver_path:
        pytest.skip("CFDX_PRODUCTION_SOLVER is provided by ctest")
    solver = Path(solver_path)
    assert solver.is_file()

    root = Path(__file__).resolve().parents[2]
    source = root / "tests" / "data" / "su2" / "mesh_NACA0012_inv.su2"
    assert source.is_file()

    from cfdx.io.converter import convert, write_result
    from cfdx.io.numerical_selection import build_numerical_selection

    for model in SCHUR_MODELS:
        case_path = tmp_path / f"n9_{model}.cfdx.h5"
        result = convert(source, solver="su2", output=None, write_reports=False)
        assert result.case is not None

        numerics = result.case.numerics
        numerics.coupled_solver = "coupled"
        numerics.coupled_schur = model
        numerics.linear_solver = "fgmres"
        numerics.gradient_operator = "green_gauss_cell"
        numerics.momentum_scheme = "first_order"
        numerics.selection = build_numerical_selection(numerics)

        assert not any(
            finding.severity.value == "unsupported_blocking"
            for finding in result.gap_analysis.findings
        )
        write_result(result, case_path, write_reports=False)

        with h5py.File(case_path, "r") as h5:
            setup = json.loads(h5.attrs["case_setup_json"])
            entries = {
                (entry["family"], entry["configuration_key"])
                for entry in setup["numerics"]["selection"]["entries"]
            }
            assert ("pressure_velocity", "pressure_velocity.coupled") in entries
            assert ("schur", f"schur.{model}") in entries
            report = h5.attrs["numerical_selection_report"]
            if isinstance(report, bytes):
                report = report.decode()
            assert f"scheme[schur]=schur.{model}" in report

        output_dir = tmp_path / f"run_{model}"
        completed = subprocess.run(
            [
                str(solver),
                "--mesh",
                str(case_path),
                "--output-dir",
                str(output_dir),
                "--iterations",
                "1",
            ],
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
        )
        diagnostics = completed.stdout + completed.stderr

        # One iteration is intentionally used here: this test verifies the
        # production selection boundary, not convergence of N9 physical cases.
        assert completed.returncode in (0, 1), diagnostics[-12000:]
        assert "Resolved numerical selections:" in diagnostics
        assert "scheme[pressure_velocity]=pressure_velocity.coupled" in diagnostics
        assert f"scheme[schur]=schur.{model}" in diagnostics

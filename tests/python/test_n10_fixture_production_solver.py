from __future__ import annotations

import json
import os
import math
import re
import subprocess
from pathlib import Path


def test_n10_verified_fixture_reaches_production_solver(tmp_path: Path) -> None:
    """Exercise a verified public fixture through the production solver.

    The mesh-to-CFDX conversion is deliberately performed by the production
    C++ mesh importer, not by Python format auto-detection. This keeps the N10
    fixture path aligned with the independently qualified importer and avoids
    introducing generic VTK support into the Python converter as part of this
    increment.

    This is numerical-support smoke evidence, not physical validation.
    """
    fixture_root_value = os.environ.get("CFDX_N10_FIXTURE_ROOT")
    solver_value = os.environ.get("CFDX_PRODUCTION_SOLVER")
    case_builder_value = os.environ.get("CFDX_N10_PRODUCTION_CASE_BUILDER")
    evidence_value = os.environ.get("CFDX_N10_SOLVER_EVIDENCE")
    if not fixture_root_value or not solver_value or not case_builder_value:
        import pytest

        pytest.skip(
            "N10 production-solver fixture environment is provided by the dedicated workflow"
        )

    fixture_root = Path(fixture_root_value)
    solver = Path(solver_value)
    case_builder = Path(case_builder_value)
    source = fixture_root / "meshio-vtk-unstructured"

    assert source.is_file()
    assert solver.is_file()
    assert case_builder.is_file()

    vtk_source = tmp_path / "06_unstructured.vtk"
    vtk_source.write_bytes(source.read_bytes())

    case_path = tmp_path / "unstructured.cfdx.h5"
    build_case = subprocess.run(
        [str(case_builder), str(vtk_source), str(case_path)],
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
        env=os.environ.copy(),
    )
    build_diagnostics = build_case.stdout + build_case.stderr
    assert build_case.returncode == 0, build_diagnostics[-12000:]
    assert case_path.is_file()

    output_dir = tmp_path / "run"
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
        timeout=60,
        check=False,
        env=os.environ.copy(),
    )
    diagnostics = result.stdout + result.stderr
    restart_path = output_dir / "restart.dat"
    assert result.returncode == 0, diagnostics[-12000:]

    trace_pattern = re.compile(
        r"^INCOMPRESSIBLE_ITER iter=(?P<iteration>\\d+)"
        r" momentum_residual_relative=(?P<momentum>[-+0-9.eE]+)"
        r" continuity_normalized=(?P<continuity>[-+0-9.eE]+)"
        r" velocity_change_inf=(?P<velocity>[-+0-9.eE]+)"
        r" pressure_change_inf=(?P<pressure>[-+0-9.eE]+)"
        r" alpha_u=(?P<alpha_u>[-+0-9.eE]+)"
        r" alpha_p=(?P<alpha_p>[-+0-9.eE]+)"
        r" pressure_linear_iterations=(?P<pressure_iterations>\\d+)"
        r" pressure_relative_residual=(?P<pressure_residual>[-+0-9.eE]+)"
        r" momentum_linear_iterations=(?P<momentum_iterations>\\d+)"
        r" flux_velocity_mismatch_linf=(?P<flux_mismatch>[-+0-9.eE]+)"
    )
    history: list[dict[str, float | int]] = []
    for line in diagnostics.splitlines():
        match = trace_pattern.match(line.strip())
        if not match:
            continue
        row: dict[str, float | int] = {
            "iteration": int(match.group("iteration")),
            "momentum_residual_relative": float(match.group("momentum")),
            "continuity_normalized": float(match.group("continuity")),
            "velocity_change_inf": float(match.group("velocity")),
            "pressure_change_inf": float(match.group("pressure")),
            "alpha_u": float(match.group("alpha_u")),
            "alpha_p": float(match.group("alpha_p")),
            "pressure_linear_iterations": int(match.group("pressure_iterations")),
            "pressure_relative_residual": float(match.group("pressure_residual")),
            "momentum_linear_iterations": int(match.group("momentum_iterations")),
            "flux_velocity_mismatch_linf": float(match.group("flux_mismatch")),
        }
        assert all(
            math.isfinite(float(value))
            for key, value in row.items()
            if key not in {"iteration", "pressure_linear_iterations", "momentum_linear_iterations"}
        )
        history.append(row)

    assert history, diagnostics[-12000:]
    assert [int(row["iteration"]) for row in history] == sorted(
        int(row["iteration"]) for row in history
    )
    assert "Resolved numerical selections:" in diagnostics
    assert "Converged YES" in diagnostics
    assert restart_path.is_file()

    evidence = {
        "fixture": "meshio-vtk-unstructured",
        "conversion_path": "production_cpp_mesh_importer",
        "case": str(case_path),
        "production_solver": str(solver),
        "execution_status": "PASS",
        "numerical_support": "execution_and_convergence_smoke_only",
        "physical_validation": "NOT_CLAIMED",
        "restart_artifact": True,
        "iteration_budget": 20,
        "quantitative_solver_evidence": {
            "trace_count": len(history),
            "iterations": history,
            "final_iteration": int(history[-1]["iteration"]),
            "max_momentum_residual_relative": max(
                float(row["momentum_residual_relative"]) for row in history
            ),
            "max_continuity_normalized": max(
                float(row["continuity_normalized"]) for row in history
            ),
            "final_momentum_residual_relative": float(
                history[-1]["momentum_residual_relative"]
            ),
            "final_continuity_normalized": float(history[-1]["continuity_normalized"]),
            "final_velocity_change_inf": float(history[-1]["velocity_change_inf"]),
            "final_pressure_change_inf": float(history[-1]["pressure_change_inf"]),
            "final_flux_velocity_mismatch_linf": float(
                history[-1]["flux_velocity_mismatch_linf"]
            ),
        },
        "case_builder_output_tail": build_diagnostics[-4000:],
        "solver_output_tail": diagnostics[-4000:],
    }
    report = (
        Path(evidence_value)
        if evidence_value
        else output_dir / "n10_numerical_support_smoke.json"
    )
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as directory:
        test_n10_verified_fixture_reaches_production_solver(Path(directory))

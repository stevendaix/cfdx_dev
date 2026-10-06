from __future__ import annotations

import json
import os
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
    assert result.returncode == 0, diagnostics[-12000:]
    assert "Resolved numerical selections:" in diagnostics
    assert "Converged YES" in diagnostics
    assert (output_dir / "restart.dat").is_file()

    evidence = {
        "fixture": "meshio-vtk-unstructured",
        "conversion_path": "production_cpp_mesh_importer",
        "case": str(case_path),
        "production_solver": str(solver),
        "execution_status": "PASS",
        "numerical_support": "execution_and_convergence_smoke_only",
        "physical_validation": "NOT_CLAIMED",
        "case_builder_output_tail": build_diagnostics[-4000:],
        "solver_output_tail": diagnostics[-4000:],
    }
    report = output_dir / "n10_numerical_support_smoke.json"
    report.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as directory:
        test_n10_verified_fixture_reaches_production_solver(Path(directory))

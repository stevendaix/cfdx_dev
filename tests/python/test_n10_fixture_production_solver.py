from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path


def test_n10_verified_fixture_reaches_production_solver(tmp_path: Path) -> None:
    """Exercise one independently verified fixture through the production solver.

    This is a numerical-support smoke gate, not a physical accuracy claim. The
    solver must actually load the converted CFDX case and report convergence;
    no tolerance is changed by the test.
    """
    fixture_root_value = os.environ.get("CFDX_N10_FIXTURE_ROOT")
    solver_value = os.environ.get("CFDX_PRODUCTION_SOLVER")
    if not fixture_root_value or not solver_value:
        import pytest

        pytest.skip("N10 production-solver fixture environment is provided by the dedicated workflow")
    fixture_root = Path(fixture_root_value)
    solver = Path(solver_value)
    # Use the independently verified 3-D VTK fixture here. The SU2 square is a
    # 2-D mesh with zero cell volume and is suitable for importer verification,
    # but not for the 3-D incompressible production solver.
    source = fixture_root / "meshio-vtk-unstructured"

    assert source.is_file()
    assert solver.is_file()

    case_dir = tmp_path / "case"
    case_dir.mkdir()
    source_case = case_dir / "unstructured.vtk"
    shutil.copyfile(source, source_case)

    case_path = case_dir / "unstructured.cfdx.h5"
    from cfdx.io.converter import convert

    converted = convert(source_case, output=case_path, write_reports=False)
    assert converted.case is not None

    # Build a deterministic closed-cavity execution case from the verified
    # 3-D mesh: all boundaries are no-slip walls and the first imported patch
    # is the moving lid. This is a numerical execution case, not a reference
    # solution or physical qualification claim.
    from cfdx.io.schema import BCType, BCValueType

    assert converted.case.boundary_conditions
    for boundary in converted.case.boundary_conditions:
        boundary.type = BCType.WALL
        boundary.value_type = BCValueType.WALL_NO_SLIP
        boundary.velocity_vector = [0.0, 0.0, 0.0]
        boundary.velocity_magnitude = 0.0

    moving_lid = converted.case.boundary_conditions[0]
    moving_lid.velocity_vector = [1.0, 0.0, 0.0]
    converted.case.initial_condition.velocity_vector = [0.0, 0.0, 0.0]

    from cfdx.io.converter import write_result

    write_result(converted, case_path, write_reports=False)
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
        "case": str(case_path),
        "production_solver": str(solver),
        "execution_status": "PASS",
        "numerical_support": "execution_and_convergence_smoke_only",
        "physical_validation": "NOT_CLAIMED",
        "solver_output_tail": diagnostics[-4000:],
    }
    report = output_dir / "n10_numerical_support_smoke.json"
    report.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as directory:
        test_n10_verified_fixture_reaches_production_solver(Path(directory))

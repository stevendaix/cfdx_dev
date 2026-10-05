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
    fixture_root = Path(os.environ["CFDX_N10_FIXTURE_ROOT"])
    solver = Path(os.environ["CFDX_PRODUCTION_SOLVER"])
    source = fixture_root / "meshio-su2-square"

    assert source.is_file()
    assert solver.is_file()

    case_dir = tmp_path / "case"
    case_dir.mkdir()
    source_case = case_dir / "square.su2"
    shutil.copyfile(source, source_case)

    case_path = case_dir / "square.cfdx.h5"
    from cfdx.io.converter import convert

    converted = convert(source_case, output=case_path)
    assert converted.case is not None
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
        "fixture": "meshio-su2-square",
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

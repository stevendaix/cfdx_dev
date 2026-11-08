"""A real solver is required for PASS; stubs and missing binaries cannot qualify."""

from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "scripts/qualify_coupled_matrix.py"
spec = importlib.util.spec_from_file_location("cfdx_coupled_comparator", SOURCE)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_naca_classification_requires_independent_solver_success_marker() -> None:
    assert module.classify_naca(0, "did something") == "UNVERIFIED_SUCCESS_NO_PHYSICAL_MARKER"
    assert module.classify_naca(0, "NACA_SOLVE_WALL accepted_iterations=305 converged=1 algorithm=COUPLED") == "SOLVER_CONVERGED_ONE_GRID_ONLY"
    assert module.classify_naca(1, "NACA_SOLVE_WALL converged=1") == "FAILED_OR_SHORT_DIAGNOSTIC_NOT_CONVERGED"
    assert module.classify_naca(124, "NACA_SOLVE_WALL converged=1") == "TIMEOUT_NOT_CONVERGED"


def test_fp_summary_uses_true_unscaled_residuals() -> None:
    content = (
        "COUPLED_NESTED_PCD_FP original_true_relative_rhs=0.35 accepted=1\n"
        "COUPLED_NESTED_PCD_FP original_true_relative_rhs=2e-9 accepted=1\n"
        "COUPLED_NESTED_PCD_FP original_true_relative_rhs=2e-2 accepted=0\n"
    )
    data = module.summarize_fp(content)
    assert data["observed_applications"] == 3
    assert data["true_relative_rhs_max"] == 0.35
    assert data["true_relative_rhs_above_1e_3"] == 2
    assert data["quality_rejected"] == 1


def test_missing_binaries_are_nonqualifying(tmp_path: Path) -> None:
    row = module.run_case(
        label="N9_PCD", exe=tmp_path / "does_not_exist", args=[],
        env={}, output=tmp_path / "run.log", timeout=2, naca=False,
    )
    assert row["status"] == "NOT_RUN_MISSING_BINARY_OR_MESH"
    assert row["physical_qualification"] == "NOT_QUALIFIED_BY_THIS_RUNNER"

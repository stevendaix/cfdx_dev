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


def test_n9_requires_all_real_post_oracle_cases() -> None:
    model = "block_local"
    line = "STARTING COUPLED block_local\n"
    for case in module.N9_STAGES:
        line += (
            f"COMPLETED COUPLED {model} {case} iterations=150 "
            "momentum_rhs_relative=3e-9 continuity_normalized=4e-9 "
            "mass_local_linf=1e-10 face_momentum_x_relative=2e-9 "
            "face_momentum_y_relative=3e-9\n"
        )
    passing = module.summarize_n9_physical_shard(line, model, 0)
    assert passing["physical_shard_gate_complete"]
    assert passing["completed_stages"] == list(module.N9_STAGES)
    assert not module.summarize_n9_physical_shard(line, model, 124)["physical_shard_gate_complete"]
    assert not module.summarize_n9_physical_shard(line.replace("Ghia_Re100", "Ghia_BOGUS"), model, 0)["physical_shard_gate_complete"]
    assert not module.summarize_n9_physical_shard(line.replace("4e-9", "4e-4"), model, 0)["physical_shard_gate_complete"]
    face_bad = line.replace("face_momentum_y_relative=3e-9", "face_momentum_y_relative=3e-3")
    assert not module.summarize_n9_physical_shard(face_bad, model, 0)["physical_shard_gate_complete"]
    assert not module.summarize_n9_physical_shard(line.replace("COMPLETED COUPLED block_local Poiseuille", "COMPLETED COUPLED pcd Poiseuille"), model, 0)["physical_shard_gate_complete"]
    assert not module.summarize_n9_physical_shard(line + line.splitlines(1)[-1], model, 0)["physical_shard_gate_complete"]


def test_n9_zero_exit_unverified_driver_is_not_physical_success(tmp_path: Path) -> None:
    # A fake executable with a successful exit must not be accepted as a
    # CFDX physical campaign (an old N9 synthetic driver did just that).
    executable = tmp_path / "fake_n9"
    executable.write_text("#!/bin/sh\necho 'STARTING COUPLED pcd'\nexit 0\n")
    executable.chmod(0o755)
    row = module.run_case(label="N9_bad", exe=executable, args=["--schur-model", "pcd"],
                          env={}, output=tmp_path / "n9.log", timeout=5, naca=False)
    assert row["exit_code"] == 0
    assert row["status"] == "N9_UNVERIFIED_ZERO_EXIT"
    assert not row["n9_physical_cases"]["physical_shard_gate_complete"]


def test_naca_requires_full_muscl_real_physical_gates() -> None:
    base = {
        "status": "OBSERVED_NOT_INDEPENDENTLY_CERTIFIED",
        "algorithm": "COUPLED",
        "final_convection_scheme": "MUSCL_CELL_LIMITED",
        "metrics": {
            "high_order_fraction": 1.0, "momentum_rhs_relative": 3e-9,
            "momentum_component_max": 4e-9, "face_momentum_x_relative": 2e-9,
            "face_momentum_y_relative": 5e-9, "continuity_normalized": 2e-9,
            "mass_local_linf": 1e-10, "nonfinite_mass_faces": 0.0,
            "nonfinite_velocities": 0.0,
        },
    }
    assert module.final_naca_muscl_contract_valid(base)
    for key, value in (("high_order_fraction", 0.8),
                       ("momentum_rhs_relative", 2e-7),
                       ("face_momentum_y_relative", 2e-7),
                       ("mass_local_linf", float("nan")),
                       ("nonfinite_mass_faces", 1.0)):
        bad = {**base, "metrics": {**base["metrics"], key: value}}
        assert not module.final_naca_muscl_contract_valid(bad)
    assert not module.final_naca_muscl_contract_valid({**base, "final_convection_scheme": "TVD_VANLEER"})
    assert not module.final_naca_muscl_contract_valid({**base, "algorithm": "SIMPLE"})

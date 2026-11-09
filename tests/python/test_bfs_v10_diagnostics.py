"""V10 independently tests the first blocker and fail-closed NASA observables."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src/cfdx/python"))
from diagnose_bfs_failure_v10 import classify_log
from qualify_bfs_driver_seegmiller_v10 import cf_and_reattachment, make_report


def test_linear_failure_is_specific_and_uses_true_relative():
    message = (
        "BFS_PHASE_V10 iteration=1 phase=COUPLED_LINEAR_BEGIN\n"
        "BFS_PHASE_V10 iteration=1 phase=COUPLED_LINEAR_END status=2 iterations=52 true_relative=0.0106\n"
        "CFDX nonlinear iteration 1 rejected; retry=1 reason=solve_steady_incompressible: "
        "coupled momentum-continuity solve did not converge "
        "(status=2, iterations=52, residual=1.248193e-02, relative=1.055851e-02)"
    )
    diagnostic = classify_log(message, execution="FAILED", exit_code=2)
    assert diagnostic["first_blocker"] == "KRYLOV_SCHUR"
    assert diagnostic["linear_evidence"]["linear_iterations"] == 52
    assert diagnostic["linear_evidence"]["true_residual_relative"] == pytest.approx(
        0.01055851
    )
    assert diagnostic["proof_status"] == "NOT_QUALIFIED"


@pytest.mark.parametrize(
    ("phase", "expected"),
    [
        ("PRESSURE_LINEAR_BEGIN", "PRESSURE_RUNNING_TIMEOUT"),
        ("COUPLED_LINEAR_BEGIN", "KRYLOV_SCHUR_RUNNING_TIMEOUT"),
        ("SST_TRANSPORT_BEGIN", "SST_RUNNING_TIMEOUT"),
    ],
)
def test_timeout_only_reports_phase_not_divergence(phase, expected):
    record = classify_log(
        f"BFS_PHASE_V10 iteration=1 phase={phase}\n",
        execution="TIMEOUT_INCOMPLETE",
        exit_code=None,
    )
    assert record["first_blocker"] == expected
    assert record["proof_status"] == "NOT_QUALIFIED"


def test_first_failure_before_later_retry():
    log = (
        "solve_steady_incompressible: pressure-correction solve did not converge "
        "(status=1, iterations=100, residual=1.0e-13, relative=4.0e-11)\n"
        "solve_steady_incompressible: pressure-correction solve did not converge "
        "(status=3, iterations=110, residual=0.9e-13, relative=3.0e-11)"
    )
    r = classify_log(log, execution="FAILED", exit_code=2)
    assert r["first_blocker"] == "PRESSURE"
    assert r["linear_evidence"]["linear_iterations"] == 100


def test_no_false_physical_result_when_not_converged(tmp_path):
    source = tmp_path / "bad_case.cfdx.h5"
    source.write_bytes(b"not physical")
    record = make_report(
        output_dir=tmp_path / "none",
        case_path=source,
        algorithm="coupled",
        schur="pcd",
        preconditioner="pcd",
        process_exit_code=2,
        rho=1.2,
        u_ref=44.2,
        step_height=0.0127,
    )
    assert not record["v8_nonlinear_certificate_passed"]
    assert record["comparison"] is None
    assert record["physical_validation"] == "NOT_QUALIFIED"
    assert record["reynolds_shear_stresses"] == "NOT_EXPORTED_FROM_PRODUCTION_SOLVER"


def test_downstream_horizontal_wall_crossing_only(tmp_path):
    faces = tmp_path / "faces.csv"
    with faces.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "patch",
                "x_m",
                "area_y_m2",
                "area_m2",
                "wall_shear_x_Pa",
                "wall_shear_method",
            ],
        )
        writer.writeheader()
        # A vertical face at x=0 must never create spurious reattachment.
        writer.writerow(
            {
                "patch": "lowerWall",
                "x_m": 0,
                "area_y_m2": 0,
                "area_m2": 1,
                "wall_shear_x_Pa": 0,
                "wall_shear_method": "first_order_one_sided_diagnostic",
            }
        )
        for x, tau in [(3, -3), (4, -2), (5, -1), (6, 1), (7, 2)]:
            writer.writerow(
                {
                    "patch": "lowerWall",
                    "x_m": x * 0.0127,
                    "area_y_m2": 1,
                    "area_m2": 1,
                    "wall_shear_x_Pa": tau,
                    "wall_shear_method": "first_order_one_sided_diagnostic",
                }
            )
    r = cf_and_reattachment(faces, rho=1.2, u_ref=44.2, step_height=0.0127)
    assert r["reattachment_x_over_h_diagnostic"] == pytest.approx(5.5)
    assert len(r["cf_series"]) == 5
    assert r["physical_certification"].startswith("NOT_QUALIFIED")


def test_missing_reattachment_is_not_assumed(tmp_path):
    faces = tmp_path / "faces.csv"
    with faces.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "patch",
                "x_m",
                "area_y_m2",
                "area_m2",
                "wall_shear_x_Pa",
                "wall_shear_method",
            ],
        )
        writer.writeheader()
        for x in (1, 2, 3, 4):
            writer.writerow(
                {
                    "patch": "lowerWall",
                    "x_m": x * 0.0127,
                    "area_y_m2": 1,
                    "area_m2": 1,
                    "wall_shear_x_Pa": -2,
                    "wall_shear_method": "first_order_one_sided_diagnostic",
                }
            )
    r = cf_and_reattachment(faces, rho=1.2, u_ref=44.2, step_height=0.0127)
    assert r["reattachment_x_over_h_diagnostic"] is None


def test_resolved_solver_is_observed_independently_of_requested_controls():
    log = (
        "BFS_PHASE_V10 iteration=1 phase=COUPLED_LINEAR_BEGIN\n"
        "BFS_PHASE_V10 iteration=1 phase=COUPLED_LINEAR_END status=2 "
        "iterations=52 true_relative=0.0106 actual_krylov=fgmres "
        "actual_preconditioner=pcd resolved_schur=pcd\n"
    )
    r = classify_log(log, execution="TIMEOUT_INCOMPLETE", exit_code=None)
    assert r["first_blocker"] == "UNRESOLVED_SETUP_OR_ASSEMBLY_TIMEOUT"
    assert r["first_resolved_linear_plan"] == {
        "iteration": 1,
        "actual_krylov": "fgmres",
        "actual_preconditioner": "pcd",
        "resolved_schur": "pcd",
    }

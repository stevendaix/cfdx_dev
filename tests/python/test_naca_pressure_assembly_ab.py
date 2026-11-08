"""Fail-closed NACA physical output contract for compensated pressure assembly."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "naca_pressure_assembly_ab.py"
spec = importlib.util.spec_from_file_location("naca_pressure_assembly_ab", SCRIPT)
assert spec and spec.loader
ab = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ab)


def example(case="N64", cd=0.11, scheme_order=1.0):
    return {"case": case, "h": 1/64, "iterations": 305,
            "force": {"cd": cd, "cdp": 0.02, "cdv": cd - .02, "cl": 0.0,
                      "cm": 0.0, "res_mom": 1e-10, "res_mass": 1e-10},
            "final_high_order_fraction": scheme_order,
            "final_momentum_component_max": 1e-9,
            "final_face_momentum_y_relative": 1e-9}


def test_comparison_requires_same_grid():
    with pytest.raises(ValueError, match="different mesh"):
        ab.compare(example(), example(case="N128"), 1e-8)


def test_comparison_catches_changed_forces():
    report = ab.compare(example(), example(cd=0.12), 5e-7)
    assert report["status"] == "ONE_GRID_REGRESSION_FAIL"
    assert report["complete_three_grid_gci"] is False


def test_comparison_passes_only_same_physics():
    report = ab.compare(example(), example(cd=0.11000000001), 5e-7)
    assert report["status"] == "ONE_GRID_REGRESSION_PASS"


def test_comparison_rejects_first_order():
    with pytest.raises(ValueError, match="order-two"):
        ab.compare(example(), example(scheme_order=0.5), 5e-7)


def test_missing_converged_csv_fails(tmp_path):
    (tmp_path / "run.log").write_text("NACA_SOLVE_WALL accepted_iterations=3 converged=0 algorithm=SIMPLE\n")
    with pytest.raises(ValueError, match="physically converge"):
        ab.parse_results(tmp_path)


def test_rejects_nonfinite_tolerance():
    with pytest.raises(ValueError, match="comparison tolerance"):
        ab.compare(example(), example(), float("nan"))


def test_five_argument_single_replay_uses_its_actual_certificate(tmp_path):
    # The production --single mesh alpha_u alpha_p max_iterations branch emits
    # NACA_CONVERGED_FORCE_CERTIFICATE, not the unrelated legacy SINGLE: PASS.
    log = ("NACA_SOLVE_WALL elapsed_seconds=0.1 accepted_iterations=1 converged=1 "
           "algorithm=SIMPLE pseudo_momentum=0\n"
           "NACA_FARFIELD_FLUX_AUDIT faces=1 kinematic_match=1 physically_converged=1\n"
           "NACA_CONVERGED_FORCE_CERTIFICATE path=force.csv qualification=ONE_GRID_ONLY\n"
           "NACA0012_RE1000_A0 level=64 Cd=0.1 iterations=1\n")
    (tmp_path / "run.log").write_text(log)
    (tmp_path / "force.csv").write_text(
        "case,h,converged,res_mom,res_mass,cd,cdp,cdv,cl,cm\n"
        "N64,0.015625,yes,1e-11,1e-11,0.1,0.02,0.08,0,0\n")
    (tmp_path / "history.csv").write_text(
        "iteration,algorithm,convection_scheme,high_order_fraction,"
        "momentum_component_max,face_momentum_y_relative\n"
        "1,SIMPLE,MUSCL_CELL_LIMITED,1.0,1e-11,1e-11\n")
    result = ab.parse_results(tmp_path)
    assert result["case"] == "N64" and result["iterations"] == 1
    assert result["force"]["cd"] == pytest.approx(0.1)
    (tmp_path / "run.log").write_text(log.replace("iterations=1\n", "iterations=2\n"))
    with pytest.raises(ValueError, match="iteration mismatch"):
        ab.parse_results(tmp_path)

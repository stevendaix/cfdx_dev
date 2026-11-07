"""Contract tests for fail-closed NACA GCI reporting."""

import importlib.util
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "naca_three_grid_gci", ROOT / "scripts" / "naca_three_grid_gci.py"
)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def data(cds=(0.14, 0.13, 0.1275)):
    return {
        "three_grid_solver_gate": "PASS",
        "results": [
            {
                "N": n,
                "status": "PASS",
                "metrics": {
                    "Cd": cd,
                    "Cl": 0.0,
                    "Cd_pressure": 0.05,
                    "Cd_viscous": cd - 0.05,
                    "continuity": 1e-9,
                },
            }
            for n, cd in zip((64, 128, 256), cds)
        ],
    }


def test_convergent_three_grid_example():
    report = MODULE.verify_three_grid(data())
    assert math.isclose(report["observed_order"], 2.0)
    assert math.isclose(report["richardson_extrapolated_cd"], 0.12666666666666668)
    assert report["fine_grid_relative_gci"] > 0
    assert report["status"] == "DIAGNOSTIC_ONLY"
    assert report["reference_validation"] == "NOT_PERFORMED"


def test_missing_n128_cannot_promote_success():
    d = data()
    d["results"][1]["status"] = "FAIL"
    with pytest.raises(ValueError):
        MODULE.verify_three_grid(d)


def test_nonfinite_forces_fail():
    d = data()
    d["results"][2]["metrics"]["Cd"] = float("nan")
    with pytest.raises(ValueError):
        MODULE.verify_three_grid(d)


def test_oscillatory_refinement_fails():
    with pytest.raises(ValueError):
        MODULE.verify_three_grid(data((0.14, 0.13, 0.135)))


def test_no_solver_success_no_gci():
    d = data()
    d["three_grid_solver_gate"] = "FAIL"
    with pytest.raises(ValueError):
        MODULE.verify_three_grid(d)


def test_negative_continuity_is_not_physical_convergence():
    d = data()
    d["results"][1]["metrics"]["continuity"] = -1e-9
    with pytest.raises(ValueError):
        MODULE.verify_three_grid(d)


def test_nonconservative_drag_partition_rejected():
    d = data()
    d["results"][2]["metrics"]["Cd_viscous"] = 0.0
    with pytest.raises(ValueError):
        MODULE.verify_three_grid(d)

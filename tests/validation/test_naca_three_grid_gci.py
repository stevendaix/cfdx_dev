"""Contract tests for fail-closed NACA GCI reporting."""

import importlib.util
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
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
                    "Cd_pressure": 0.05 + (cd - 0.12)*0.4,
                    "Cd_viscous": cd - 0.05 - (cd - 0.12)*0.4,
                    "Cm": 0.0,
                    "continuity": 1e-9,
                },
            }
            for n, cd in zip((64, 128, 256), cds)
        ],
        "mesh_quality": [
            {"n_surface": n, "n_radial": 3*n//8,
             "radius": 20.0, "radial_power": 1.8,
             "minimum_first_layer_spacing": 1/n,
             "inverted_or_degenerate_cells": 0}
            for n in (64, 128, 256)
        ],
    }


def test_convergent_three_grid_example():
    report = MODULE.verify_three_grid(data())
    assert math.isclose(report["observed_order"], 2.0)
    assert math.isclose(report["richardson_extrapolated_cd"], 0.12666666666666668)
    assert report["fine_grid_relative_gci"] > 0
    assert report["reference_validation"] == "NOT_PERFORMED"
    assert report["quantities"]["Cd_pressure"]["status"] == "RICHARDSON_CANDIDATE"
    assert report["quantities"]["Cd_viscous"]["status"] == "RICHARDSON_CANDIDATE"
    assert report["quantities"]["Cl"]["gci"] is None


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


def test_no_mesh_quality_no_gci():
    d = data()
    del d["mesh_quality"]
    with pytest.raises(ValueError, match="mesh-quality"):
        MODULE.verify_three_grid(d)


def test_original_naca_wall_refinement_blocks_gci():
    d = data()
    for row, spacing in zip(d["mesh_quality"], (0.049996, 0.015458, 0.004733)):
        row["minimum_first_layer_spacing"] = spacing
    with pytest.raises(ValueError, match="not similar"):
        MODULE.verify_three_grid(d)


def test_pressure_viscous_partition_must_close():
    d = data()
    d["results"][2]["metrics"]["Cd_viscous"] = -42.0
    with pytest.raises(ValueError, match="do not sum"):
        MODULE.verify_three_grid(d)

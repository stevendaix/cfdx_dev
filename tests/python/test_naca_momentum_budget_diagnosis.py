from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import pytest

SCRIPT = (
    Path(__file__).resolve().parents[2] / "scripts/naca_momentum_budget_diagnosis.py"
)
spec = importlib.util.spec_from_file_location("naca_momentum_budget_diagnosis", SCRIPT)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def make_file(tmp_path: Path, **changes) -> Path:
    p = tmp_path / "budget.csv"
    row = {
        "iteration": "10",
        "component": "1",
        "rank": "0",
        "cell": "55",
        "x": "0.8",
        "y": "0.13",
        "z": "0.5",
        "convection": "0.15",
        "viscous": "-0.03",
        "transpose_viscous": "0",
        "pressure": "-0.10",
        "body_force": "-0.015",
        "bounded_divergence": "-0.004",
        "mass_divergence": "0.002",
        "collocated_difference": "0.003",
        "physical_residual": "0.005",
        "bounded_operator_residual": "0.001",
        "worst_residual": "0.005",
    }
    row.update(changes)
    with p.open("w", encoding="utf8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=module.COLUMNS)
        writer.writeheader()
        writer.writerow(row)
    return p


def test_accepts_signed_physical_budget_but_does_not_certify(tmp_path):
    data = module.diagnose(module.read_budget(make_file(tmp_path)))
    assert data["solver_convergence_proven"] is False
    assert data["three_grid_gci_proven"] is False
    assert data["budget_series"][0]["worst_cell"] == 55


@pytest.mark.parametrize(
    "change",
    [
        {"convection": "0.2"},
        {"bounded_operator_residual": "0.2"},
        {"physical_residual": "nan"},
        {"rank": "2"},
        {"worst_residual": "-0.005"},
        {"iteration": "0"},
        {"component": "3"},
        {"cell": "1.1"},
        {"pressure": "inf"},
    ],
)
def test_rejects_corrupt_evidence(tmp_path, change):
    with pytest.raises(ValueError):
        module.read_budget(make_file(tmp_path, **change))

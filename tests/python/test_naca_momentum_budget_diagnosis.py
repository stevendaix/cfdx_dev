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


def test_accepts_cellwise_algebraic_closure(tmp_path):
    file = make_file(tmp_path)
    rows = module.read_budget(file)
    assert module.diagnose(rows)["matrix_face_closure_checked"] is False
    with file.open(newline="", encoding="utf8") as stream:
        parsed = list(csv.DictReader(stream))
    row = parsed[0]
    row.update(
        temporal="0",
        algebraic_residual="0.001",
        matrix_face_difference="0",
        matrix_face_scale="0.5",
    )
    with file.open("w", newline="", encoding="utf8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=[*module.COLUMNS, *module.EXTENDED_COLUMNS]
        )
        writer.writeheader()
        writer.writerow(row)
    report = module.diagnose(module.read_budget(file))
    assert report["matrix_face_closure_checked"]
    assert report["max_matrix_face_error"] == 0.0
    assert (
        report["budget_series"][0]["dominant_terms_by_absolute_magnitude"][0]
        == "convection"
    )


@pytest.mark.parametrize(
    "field, value",
    [
        ("algebraic_residual", "0.003"),
        ("matrix_face_difference", "0.002"),
        ("matrix_face_scale", "-1"),
        ("temporal", "0.2"),
    ],
)
def test_rejects_bad_algebraic_closure(tmp_path, field, value):
    file = make_file(tmp_path)
    with file.open(newline="", encoding="utf8") as stream:
        row = list(csv.DictReader(stream))[0]
    row.update(
        temporal="0",
        algebraic_residual="0.001",
        matrix_face_difference="0",
        matrix_face_scale="0.5",
    )
    row[field] = value
    with file.open("w", newline="", encoding="utf8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=[*module.COLUMNS, *module.EXTENDED_COLUMNS]
        )
        writer.writeheader()
        writer.writerow(row)
    with pytest.raises(ValueError):
        module.read_budget(file)


def test_period_two_requires_unit_stride_and_fixed_cell():
    rows = [
        {
            "iteration": k + 1,
            "component": 1,
            "cell": 184,
            "physical_residual": (1 if k % 2 == 0 else -1) * 4.3e-5,
        }
        for k in range(10)
    ]
    result = module.sampled_sign_alternation(rows, 1, 184)
    assert result["status"] == "PERIOD_TWO_SIGNATURE"
    rows = [{**row, "iteration": row["iteration"] * 15} for row in rows]
    result = module.sampled_sign_alternation(rows, 1, 184)
    assert result["status"] == "ALTERNATING_SAMPLES_PERIOD_UNPROVEN"
    rows[4]["cell"] = 201
    assert (
        module.sampled_sign_alternation(rows, 1, 184)["status"] == "IRREGULAR_SAMPLING"
    )


def test_monotone_signed_residual_is_not_period_two():
    rows = [
        {
            "iteration": k + 1,
            "component": 1,
            "cell": 184,
            "physical_residual": 1e-3 / (k + 1),
        }
        for k in range(10)
    ]
    assert (
        module.sampled_sign_alternation(rows, 1, 184)["status"]
        == "NO_STABLE_ALTERNATION"
    )


def test_period_two_signed_term_attribution_and_cancellation():
    rows = []
    for i in range(10):
        even = i % 2 == 0
        convection = 5e-4 if even else 4e-4
        viscous = -2e-5 if even else -1e-5
        pressure = -4e-4 if even else -4.8e-4
        physical = convection + viscous + pressure
        rows.append(
            {
                "iteration": i + 1,
                "component": 1,
                "cell": 184,
                "convection": convection,
                "viscous": viscous,
                "transpose_viscous": 0.0,
                "pressure": pressure,
                "body_force": 0.0,
                "temporal": 0.0,
                "bounded_divergence": -1e-6,
                "collocated_difference": 2e-6,
                "physical_residual": physical,
                "bounded_operator_residual": physical - 1e-6,
            }
        )
    outcome = module.period_two_term_decomposition(rows, 1, 184)
    assert outcome["term_decomposition_available"]
    assert outcome["term_delta_closure_error"] == pytest.approx(0, abs=1e-16)
    assert outcome["terms_by_absolute_cycle_amplitude"][0] == "convection"
    assert outcome["physical_cycle_delta"] == pytest.approx(0.00017)


def test_cycle_attribution_not_inferred_from_sparse_samples():
    rows = [
        {
            "iteration": 15 * (i + 1),
            "component": 1,
            "cell": 184,
            "physical_residual": ((-1) ** i) * 4.3e-5,
        }
        for i in range(8)
    ]
    report = module.period_two_term_decomposition(rows, 1, 184)
    assert report["signature"]["status"] == "ALTERNATING_SAMPLES_PERIOD_UNPROVEN"
    assert report["term_decomposition_available"] is False

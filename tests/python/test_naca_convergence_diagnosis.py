"""Fail-closed convergence history diagnostics, never artificial solver PASS."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = ROOT / "scripts/naca_convergence_diagnosis.py"
SPEC = importlib.util.spec_from_file_location("naca_diag", MODULE_PATH)
assert SPEC and SPEC.loader
DIAG = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DIAG)


def history(tmp_path: Path, **changes) -> Path:
    names = (*DIAG.REQUIRED,)
    rows = []
    for i in range(1, 5):
        values = {key: "1e-11" for key in names}
        values.update(
            {
                "iteration": str(i),
                "convection_scheme": "TVD",
                "algorithm": "SIMPLE",
                "high_order_fraction": "1",
                "momentum_component_max": "1e-11",
            }
        )
        values.update(changes)
        rows.append(values)
    output = tmp_path / "history.csv"
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)
    return output


def test_good_history_stays_diagnostic_not_certified(tmp_path):
    report = DIAG.diagnose(DIAG.read_history(history(tmp_path)))
    assert report["physical_tolerance_checks_pass"]
    assert report["full_order_scheme"]
    assert report["solver_convergence_proven"] is False
    assert report["three_grid_gci_proven"] is False


def test_multidimensional_muscl_is_diagnostic_not_a_certificate(tmp_path):
    rows = DIAG.read_history(history(tmp_path, convection_scheme="MUSCL_CELL_LIMITED"))
    result = DIAG.diagnose(rows)
    assert result["full_order_scheme"]
    assert result["experimental_multidimensional_limiter"]
    assert result["solver_convergence_proven"] is False
    assert result["three_grid_gci_proven"] is False


def test_transverse_momentum_detected_independently(tmp_path):
    path = history(
        tmp_path, momentum_y_rhs_relative="2e-4", momentum_component_max="2e-4"
    )
    report = DIAG.diagnose(DIAG.read_history(path))
    assert not report["physical_tolerance_checks_pass"]
    assert report["unmet_nonlinear_checks"]["momentum_y_rhs_relative"] == 2e-4


@pytest.mark.parametrize(
    "field,value",
    [
        ("momentum_reference_relative", "nan"),
        ("face_momentum_y_relative", "-1"),
        ("continuity_normalized", "inf"),
        ("high_order_fraction", "1.1"),
        ("iteration", "1.5"),
        ("algorithm", "COUPLED"),
    ],
)
def test_bad_history_fails_closed(tmp_path, field, value):
    with pytest.raises(ValueError):
        DIAG.read_history(history(tmp_path, **{field: value}))


def test_first_order_not_full_order_candidate(tmp_path):
    report = DIAG.diagnose(
        DIAG.read_history(history(tmp_path, convection_scheme="UPWIND"))
    )
    assert report["physical_tolerance_checks_pass"]
    assert report["full_order_scheme"] is False


def test_missing_independent_metric_rejected(tmp_path):
    path = history(tmp_path)
    path.write_text("iteration,momentum_rhs_relative\n1,1e-11\n")
    with pytest.raises(ValueError):
        DIAG.read_history(path)

"""Synthetic negative controls for the forensic CFDX convergence gate."""

import copy
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "tools/diagnostics/qualify_convergence_history.py"
spec = importlib.util.spec_from_file_location("cfdx_history_audit", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

LIMITS = {
    "momentum_relative": 1e-8,
    "pressure_relative": 1e-8,
    "continuity_relative": 1e-8,
    "mass_relative": 1e-8,
}


def valid_evidence():
    return {
        "format": "CFDX-CONVERGENCE",
        "schema_version": 1,
        "converged": True,
        "iterations": 6,
        "history": [
            {
                "iteration": iteration,
                "high_order_fraction": fraction,
                "full_order_eligible": iteration >= 6,
                "momentum_equation_residual_rhs_active_component_max": 1e-9,
                "pressure_residual": 2e-9,
                "continuity_normalized": 3e-9,
                "mass_normalized_imbalance": 4e-9,
            }
            for iteration, fraction in enumerate([0.0, 0.5, 1, 1, 1, 1], start=1)
        ],
    }


def test_valid_physical_history_reports_limited_qualification():
    result = module.audit_history(valid_evidence(), LIMITS, 4)
    assert result["qualified"]
    assert result["full_order_tail"] == 4
    assert result["spatial_convergence"] == "not_evaluated"


def test_incomplete_high_order_refused():
    bad = valid_evidence()
    bad["history"][-1]["high_order_fraction"] = 0.75
    bad["history"][-1]["full_order_eligible"] = False
    assert not module.audit_history(bad, LIMITS, 4)["qualified"]


def test_false_eligibility_fails_closed():
    bad = valid_evidence()
    bad["history"][1]["full_order_eligible"] = True
    try:
        module.audit_history(bad, LIMITS, 4)
    except ValueError as exc:
        assert "intermediate" in str(exc)
    else:
        raise AssertionError("accepted eligibility before full-order operator")


def test_momentum_component_failure_refused():
    bad = valid_evidence()
    bad["history"][-1]["momentum_equation_residual_rhs_active_component_max"] = 4e-5
    assert (
        "momentum_relative_above_limit"
        in module.audit_history(bad, LIMITS, 4)["failures"]
    )


def test_mass_balance_failure_refused():
    bad = valid_evidence()
    bad["history"][-1]["mass_normalized_imbalance"] = 1e-3
    assert (
        "mass_relative_above_limit" in module.audit_history(bad, LIMITS, 4)["failures"]
    )


def test_nonfinite_and_missing_fields_refused():
    bad = valid_evidence()
    bad["history"][-1]["mass_normalized_imbalance"] = float("nan")
    try:
        module.audit_history(bad, LIMITS, 4)
    except ValueError:
        pass
    else:
        raise AssertionError("accepted NaN")
    del bad["history"][-1]["mass_normalized_imbalance"]
    try:
        module.audit_history(bad, LIMITS, 4)
    except ValueError:
        pass
    else:
        raise AssertionError("accepted missing mass evidence")


def test_false_solver_verdict_rejected():
    bad = valid_evidence()
    bad["converged"] = False
    assert "solver_not_converged" in module.audit_history(bad, LIMITS, 4)["failures"]


def test_historical_artifact_without_stage_certificate_is_not_promoted():
    bad = valid_evidence()
    del bad["history"][-1]["high_order_fraction"]
    try:
        module.audit_history(bad, LIMITS, 4)
    except ValueError:
        pass
    else:
        raise AssertionError("accepted historical data without ramp evidence")


def test_cli_returns_failure_exit_code(tmp_path):
    path = tmp_path / "history.json"
    path.write_text(json.dumps(valid_evidence()), encoding="utf-8")
    args = [
        sys.executable,
        str(SCRIPT),
        str(path),
        "--momentum-relative",
        "1e-8",
        "--pressure-relative",
        "1e-8",
        "--continuity-relative",
        "1e-8",
        "--mass-relative",
        "1e-8",
    ]
    ok = subprocess.run(args, capture_output=True, text=True, check=False)
    assert ok.returncode == 0, ok.stdout + ok.stderr
    assert json.loads(ok.stdout)["qualified"]
    bad = copy.deepcopy(valid_evidence())
    bad["converged"] = False
    path.write_text(json.dumps(bad), encoding="utf-8")
    rejected = subprocess.run(args, capture_output=True, text=True, check=False)
    assert rejected.returncode == 1
    assert not json.loads(rejected.stdout)["qualified"]

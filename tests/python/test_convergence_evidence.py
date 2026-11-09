"""Strict tests for persisted execution and convergence evidence."""

import json

import pytest

from cfdx.convergence_evidence import (
    EvidenceError,
    baseline_signatures,
    load_run_evidence,
)


def _v8_certificate() -> dict[str, object]:
    return {
        "evaluated": True,
        "algorithm": "coupled",
        "schur": "pcd",
        "momentum_global_relative": 1.0e-6,
        "momentum_max_active_component_relative": 2.0e-6,
        "momentum_component_relative": [1.0e-6, 2.0e-6, None],
        "continuity": 1.0e-6,
        "k": 1.0e-5,
        "omega": 1.0e-5,
        "eddy_viscosity_closure": 1.0e-5,
        "merit": 1.0e-5,
        "sst": True,
        "sst_transport_committed": True,
        "bounded": True,
        "all_finite": True,
        "converged": False,
    }


def _load(tmp_path, certificate):
    baseline = baseline_signatures(tmp_path)
    execution = {
        "format": "CFDX-EXECUTION",
        "schema_version": 1,
        "converged": False,
        "iterations": 1,
        "process_exit_code": 1,
        "convergence_reason": "iteration budget exhausted",
        "convergence_status": 1,
        "artifacts": {"convergence_json": True},
    }
    convergence = {
        "format": "CFDX-CONVERGENCE",
        "schema_version": 1,
        "converged": False,
        "iterations": 1,
        "history": [
            {
                "iteration": 1,
                "momentum_residual": 1.0e-3,
                "global_nonlinear_v8": certificate,
            }
        ],
    }
    (tmp_path / "execution.json").write_text(json.dumps(execution), encoding="utf-8")
    (tmp_path / "convergence.json").write_text(
        json.dumps(convergence), encoding="utf-8"
    )
    return load_run_evidence(tmp_path, process_exit_code=1, baseline=baseline)


def test_load_run_evidence_accepts_strict_v8_certificate(tmp_path):
    evidence = _load(tmp_path, _v8_certificate())
    assert evidence.flow_history_count == 1
    assert evidence.overall_converged is False
    assert dict(evidence.final_metrics)["momentum_residual"] == pytest.approx(1.0e-3)


def test_load_run_evidence_rejects_unknown_v8_certificate_keys(tmp_path):
    certificate = _v8_certificate()
    certificate["unreviewed_claim"] = True
    with pytest.raises(EvidenceError, match="invalid keys"):
        _load(tmp_path, certificate)


def test_load_run_evidence_rejects_nested_non_numeric_v8_metric(tmp_path):
    certificate = _v8_certificate()
    certificate["momentum_global_relative"] = {"value": 1.0e-6}
    with pytest.raises(EvidenceError, match="finite non-negative number"):
        _load(tmp_path, certificate)


def test_load_run_evidence_rejects_inconsistent_v8_sst_flags(tmp_path):
    certificate = _v8_certificate()
    certificate["sst"] = False
    with pytest.raises(EvidenceError, match="requires sst"):
        _load(tmp_path, certificate)


def test_load_run_evidence_rejects_missing_v8_certificate(tmp_path):
    with pytest.raises(EvidenceError, match="global_nonlinear_v8 must be an object"):
        _load(tmp_path, None)

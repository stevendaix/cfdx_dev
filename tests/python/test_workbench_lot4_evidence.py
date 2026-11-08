"""Run-scoped proof contract: both production artifacts, never a clean exit."""

import json
import sys
from pathlib import Path

import pytest

from cfdx import CFDXSession, ExecutionController, SolverRunner
from cfdx.convergence_evidence import (
    EvidenceError, baseline_signatures, load_run_evidence,
)


def _documents(*, code=0, flow=True, overall=True, iterations=3):
    execution = {
        "format": "CFDX-EXECUTION", "schema_version": 1,
        "process_exit_code": code, "converged": overall,
        "iterations": iterations, "convergence_reason": "native solver stage gate",
        "convergence_status": 2 if overall else -5,
        "artifacts": {"restart_dat": True, "convergence_json": True},
    }
    convergence = {
        "format": "CFDX-CONVERGENCE", "schema_version": 1,
        "converged": flow, "iterations": iterations,
        "history": [
            {"iteration": index, "momentum_equation_residual_relative": 10 ** -index,
             "continuity_normalized": 10 ** (-index - 1),
             "velocity_change_inf": 0.001 / index}
            for index in range(1, iterations + 1)
        ],
    }
    return execution, convergence


def _write(root: Path, execution, flow):
    root.mkdir(exist_ok=True)
    (root / "execution.json").write_text(json.dumps(execution), encoding="utf-8")
    (root / "convergence.json").write_text(json.dumps(flow), encoding="utf-8")


def test_fresh_native_artifacts_expose_physical_evidence_without_qualifying_it(tmp_path):
    baseline = baseline_signatures(tmp_path)
    _write(tmp_path, *_documents())
    evidence = load_run_evidence(tmp_path, process_exit_code=0, baseline=baseline)
    assert evidence.overall_converged
    assert evidence.iterations == 3
    assert dict(evidence.final_metrics)["continuity_normalized"] == pytest.approx(1e-4)
    assert len(evidence.execution_sha256) == len(evidence.convergence_sha256) == 64
    assert evidence.to_dict()["independent_physical_qualification"] == "not_evaluated"
    assert evidence.source == "CFDX-EXECUTION v1 + CFDX-CONVERGENCE v1"


def test_stage_fail_cannot_be_masked_by_converged_incompressible_flow(tmp_path):
    baseline = baseline_signatures(tmp_path)
    _write(tmp_path, *_documents(code=1, flow=True, overall=False))
    evidence = load_run_evidence(tmp_path, process_exit_code=1, baseline=baseline)
    assert not evidence.overall_converged
    assert evidence.flow_converged


@pytest.mark.parametrize("mutation, error", [
    (lambda e, f: e.update({"format": "OTHER"}), "format"),
    (lambda e, f: e.update({"converged": "true"}), "boolean"),
    (lambda e, f: e.update({"process_exit_code": 1}), "disagrees"),
    (lambda e, f: e.update({"iterations": 100}), "iteration counts"),
    (lambda e, f: e.update({"artifacts": {"convergence_json": False}}), "artifact"),
    (lambda e, f: f.update({"history": []}), "no final flow metrics"),
    (lambda e, f: f["history"].append(dict(f["history"][0])), "strictly increasing"),
    (lambda e, f: f["history"][-1].update({"continuity_normalized": float('inf')}), "non-finite"),
    (lambda e, f: f.update({"schema_version": 2}), "schema"),
    (lambda e, f: f.update({"converged": False}), "cannot converge"),
    (lambda e, f: f["history"][0].update({"pressure_residual": float("inf")}), "non-finite"),
    (lambda e, f: f["history"].pop(), "final iteration"),
])
def test_artifact_contradictions_fail_closed(tmp_path, mutation, error):
    execution, flow = _documents()
    mutation(execution, flow)
    baseline = baseline_signatures(tmp_path)
    _write(tmp_path, execution, flow)
    with pytest.raises(EvidenceError, match=error):
        load_run_evidence(tmp_path, process_exit_code=0, baseline=baseline)


def test_replay_of_prior_execution_is_not_current_evidence(tmp_path):
    _write(tmp_path, *_documents())
    baseline = baseline_signatures(tmp_path)
    with pytest.raises(EvidenceError, match="previous run"):
        load_run_evidence(tmp_path, process_exit_code=0, baseline=baseline)


def test_missing_baseline_and_symlink_evidence_rejected(tmp_path):
    _write(tmp_path, *_documents())
    with pytest.raises(EvidenceError, match="baseline"):
        load_run_evidence(tmp_path, process_exit_code=0)
    other = tmp_path / "other.json"
    other.write_text((tmp_path / "execution.json").read_text())
    (tmp_path / "execution.json").unlink()
    (tmp_path / "execution.json").symlink_to(other)
    with pytest.raises(EvidenceError, match="symbolic-link"):
        baseline_signatures(tmp_path)


def test_duplicate_json_keys_rejected(tmp_path):
    base = baseline_signatures(tmp_path)
    execution, flow = _documents()
    _write(tmp_path, execution, flow)
    path = tmp_path / "execution.json"
    text = path.read_text().replace('"converged": true', '"converged": true, "converged": true')
    path.write_text(text)
    with pytest.raises(EvidenceError, match="duplicate"):
        load_run_evidence(tmp_path, process_exit_code=0, baseline=base)


def test_unbounded_data_is_rejected_before_read(tmp_path):
    baseline = baseline_signatures(tmp_path)
    execution, flow = _documents()
    execution["padding"] = "x" * 4096
    _write(tmp_path, execution, flow)
    with pytest.raises(EvidenceError, match="bounded evidence read"):
        load_run_evidence(tmp_path, process_exit_code=0, baseline=baseline, max_bytes=1024)


def _run_script(tmp_path, code, *, fresh=True, overall=True, flow=True):
    output = tmp_path / "output"
    output.mkdir(exist_ok=True)
    script = tmp_path / "fake_native.py"
    execution, convergence = _documents(code=code, overall=overall, flow=flow)
    script.write_text(
        "import json, sys\nfrom pathlib import Path\n"
        "p = Path(sys.argv[1])\n"
        + ("(p / 'execution.json').write_text(json.dumps(" + repr(execution) + "))\n"
           "(p / 'convergence.json').write_text(json.dumps(" + repr(convergence) + "))\n"
           if fresh else "")
        + "print('Iteration 3 Time = 0.5 residual p = 1e-5', flush=True)\n"
        + f"raise SystemExit({code})\n", encoding="utf-8",
    )
    session = CFDXSession()
    session.case.execution.output_dir = str(output)
    runner = SolverRunner([sys.executable, str(script), str(output)])
    controller = ExecutionController(session, runner)
    controller.start()
    assert runner._thread is not None
    runner._thread.join(timeout=8)
    return controller


def test_controller_accepts_genuine_fresh_evidence(tmp_path):
    controller = _run_script(tmp_path, 0)
    assert controller.session.state.value == "CONVERGED"
    assert controller.convergence_verdict is True
    assert controller.convergence_evidence is not None
    assert controller.monitor_series.at(3).values["p"] == pytest.approx(1e-5)


def test_controller_reports_physics_nonconvergence_as_stopped_not_crash(tmp_path):
    controller = _run_script(tmp_path, 1, overall=False, flow=True)
    assert controller.session.state.value == "STOPPED"
    assert controller.convergence_verdict is False
    assert controller.error is None


def test_controller_does_not_reuse_previous_run_artifacts(tmp_path):
    output = tmp_path / "output"
    output.mkdir()
    _write(output, *_documents())
    controller = _run_script(tmp_path, 0, fresh=False)
    assert controller.session.state.value == "STOPPED"
    assert controller.convergence_evidence is None
    assert "previous run" in controller.convergence_evidence_error


def test_controller_missing_artifacts_is_never_converged(tmp_path):
    controller = _run_script(tmp_path, 0, fresh=False)
    assert controller.session.state.value == "STOPPED"
    assert controller.convergence_verdict is None
    assert "missing" in controller.convergence_evidence_error


def test_invalid_pre_run_artifact_baseline_does_not_leave_session_running(tmp_path):
    target = tmp_path / "real.json"
    target.write_text("{}", encoding="utf-8")
    (tmp_path / "execution.json").symlink_to(target)
    session = CFDXSession()
    session.case.execution.output_dir = str(tmp_path)
    runner = SolverRunner([sys.executable, "-c", "pass"])
    controller = ExecutionController(session, runner)
    with pytest.raises(EvidenceError, match="symbolic-link"):
        controller.start()
    assert session.state.value == "CREATED"
    assert runner._thread is None


def test_process_spawn_failure_does_not_leave_session_running(tmp_path):
    session = CFDXSession()
    session.case.execution.output_dir = str(tmp_path)
    runner = SolverRunner([str(tmp_path / "missing-executable")])
    controller = ExecutionController(session, runner)
    with pytest.raises(OSError):
        controller.start()
    assert session.state.value == "FAILED"
    assert runner._thread is None

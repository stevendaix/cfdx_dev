"""Headless, deterministic Mission Control and recent-project regressions."""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

from cfdx import CFDXSession
from cfdx.application import Application
from cfdx.convergence_evidence import ConvergenceEvidence
from cfdx.operations_overview import build_operations_overview
from cfdx.recent_projects import normalize_recent_projects, remember_recent_project
from cfdx.session import SimulationState


def test_new_case_does_not_claim_preflight_or_solver_convergence():
    app = Application(CFDXSession())
    state = app.state
    assert not state.validation_performed
    assert not state.validation_current
    board = build_operations_overview(state)
    assert board.check_state == "not_checked"
    assert board.solver_evidence == "missing_or_unverified"
    assert board.configured_sections < board.section_count
    assert board.to_dict()["solver_evidence"]["independent_physical_qualification"] is False
    assert any(item.key == "checks" for item in board.findings)
    assert "not_a_readiness_certificate" in board.to_dict()["configuration"]
    assert json.loads(json.dumps(board.to_dict(), allow_nan=False))["schema_version"] == 1


def test_config_check_is_distinct_from_physical_qualification_and_stales_on_edits():
    app = Application(CFDXSession())
    checked = app.validate_case()
    assert checked.validation_performed and checked.validation_current
    assert not checked.validation_with_mesh
    board = build_operations_overview(checked)
    # The default empty case is not ready; a validation call records real
    # diagnostics and can never be silently promoted to solver convergence.
    assert board.check_state in ("errors", "warnings")
    assert board.solver_evidence == "missing_or_unverified"
    app.set_project_path(None, dirty=True)
    stale = app.state
    assert stale.validation_performed and not stale.validation_current
    assert build_operations_overview(stale).check_state == "stale"
    refreshed = app.validate_case()
    assert refreshed.validation_current
    app.session.case_revision += 1  # independent direct CFD revision
    assert not app.state.validation_current


def test_diagnostic_override_does_not_count_as_successful_check():
    app = Application(CFDXSession())
    app.set_diagnostics(())  # An empty table is not validation evidence.
    assert build_operations_overview(app.state).check_state == "not_checked"


def test_solver_certification_requires_verified_run_evidence():
    app = Application(CFDXSession())
    app.session.iteration = 20
    app.session.state = SimulationState.STOPPED
    missing = SimpleNamespace(convergence_evidence=None, convergence_evidence_error="missing native artifact")
    board = build_operations_overview(app.state, run=missing)
    assert board.solver_evidence == "missing_or_unverified"
    assert any(f.key == "solver.evidence" for f in board.findings)
    evidence = ConvergenceEvidence(
        overall_converged=True, flow_converged=True, iterations=20,
        reason="criteria met", convergence_status=0, final_metrics=(("mass", 0.0),),
        flow_history_count=20, execution_sha256="a" * 64,
        convergence_sha256="b" * 64, output_directory="/test/run",
    )
    verified = SimpleNamespace(convergence_evidence=evidence, convergence_evidence_error=None)
    board = build_operations_overview(app.state, run=verified)
    assert board.solver_evidence == "converged"
    assert board.solver_evidence_sha256 == ("a" * 64, "b" * 64)
    assert not board.to_dict()["solver_evidence"]["independent_physical_qualification"]
    failed = SimpleNamespace(convergence_evidence=SimpleNamespace(
        overall_converged=False, reason="momentum stalled",
        execution_sha256="x", convergence_sha256="y"), convergence_evidence_error=None)
    board = build_operations_overview(app.state, run=failed)
    assert board.solver_evidence == "not_converged"
    assert any(row.key == "solver.not_converged" for row in board.findings)


def test_mesh_and_result_counts_are_actual_snapshot_values():
    app = Application(CFDXSession())
    absent = build_operations_overview(app.state)
    assert not absent.mesh_catalog_loaded
    assert absent.results_frames == 0
    assert any(f.key == "mesh.not_inspected" for f in absent.findings)
    loaded = build_operations_overview(app.state, mesh_catalog_loaded=True)
    assert loaded.mesh_catalog_loaded
    assert all(f.key != "mesh.not_inspected" for f in loaded.findings)


def test_recent_cases_bounded_deduplicated_and_only_existing(tmp_path: Path):
    case_a = tmp_path / "a.cfdx.h5"
    case_b = tmp_path / "b.h5"
    case_a.write_bytes(b"x")
    case_b.write_bytes(b"y")
    missing = tmp_path / "missing.cfdx.h5"
    alias = tmp_path / "linked.cfdx.h5"
    alias.symlink_to(case_a)
    rows = normalize_recent_projects([
        str(case_a), str(case_a), str(case_b), str(missing), str(alias),
        "relative.cfdx.h5", "", 17, None,
    ])
    assert rows == (str(case_a), str(case_b))
    assert remember_recent_project(case_b, rows) == (str(case_b), str(case_a))
    assert remember_recent_project(missing, rows) == rows
    assert normalize_recent_projects("not a list") == ()
    assert normalize_recent_projects([str(case_a)] * 100) == (str(case_a),)


def test_manual_diagnostic_override_invalidates_preflight():
    app = Application(CFDXSession())
    assert app.validate_case().validation_current
    app.set_diagnostics(())
    assert app.state.validation_performed
    assert not app.state.validation_current


def test_recent_project_limit_and_suffix_filter(tmp_path: Path):
    files = [tmp_path / f"project-{n}.cfdx.h5" for n in range(14)]
    for file in files:
        file.write_bytes(b"test")
    seen = normalize_recent_projects([str(file) for file in files])
    assert len(seen) == 8
    assert seen[0] == str(files[0])
    assert seen[-1] == str(files[7])
    other = tmp_path / "notes.txt"
    other.write_text("not an HDF5 project")
    assert normalize_recent_projects([str(other)]) == ()


def test_new_project_emits_result_and_selection_reset_events():
    from cfdx.application import ResultsChanged, SelectionChanged
    app = Application(CFDXSession())
    observed = []
    app.events.subscribe(ResultsChanged, lambda event: observed.append(("results", event.state)))
    app.events.subscribe(SelectionChanged, lambda event: observed.append(("selection", event.state)))
    app.new_project()
    assert [kind for kind, _state in observed] == ["results", "selection"]
    assert observed[0][1].results.frames == ()
    assert observed[1][1].selection.stable_id is None


def test_contradictory_run_verdict_is_fail_closed():
    app = Application(CFDXSession())
    evidence = ConvergenceEvidence(
        overall_converged=True, flow_converged=True, iterations=10,
        reason="converged", convergence_status=0, final_metrics=(),
        flow_history_count=10, execution_sha256="a" * 64,
        convergence_sha256="b" * 64, output_directory="/run",
    )
    run = SimpleNamespace(convergence_evidence=evidence,
                          convergence_verdict=False,
                          convergence_evidence_error=None)
    board = build_operations_overview(app.state, run=run)
    assert board.solver_evidence == "missing_or_unverified"
    assert board.solver_evidence_sha256 is None
    assert any(row.key == "solver.contradiction" and row.level == "error"
               for row in board.findings)

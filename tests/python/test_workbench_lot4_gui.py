"""Actual Qt widgets: execute in the GUI CI with PySide6 and Xvfb."""

import importlib.util
import json
import sys

import pytest


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="PySide6 optional")
def test_boundary_field_matrix_is_lossless_and_can_remove_one_field(monkeypatch):
    from PySide6.QtWidgets import QApplication, QMessageBox
    from cfdx.case import Case
    from cfdx.setup_panel import CaseSetupPanel

    app = QApplication.instance() or QApplication(["cfdx-lot4-boundaries"])
    case = Case()
    panel = CaseSetupPanel(case)
    panel.boundary_name.setText("outflow")
    panel.boundary_type.setCurrentText("outlet")
    panel.boundary_field.setCurrentText("pressure")
    panel.boundary_scalar_type.setCurrentText("FIXED_VALUE")
    panel.boundary_value.setValue(101325.0)
    panel.boundary_button.click()
    panel.boundary_field.setCurrentText("temperature")
    panel.boundary_scalar_type.setCurrentText("FIXED_VALUE")
    panel.boundary_value.setValue(307.5)
    panel.boundary_button.click()
    assert case.boundaries["outflow"]["field_conditions"]["pressure"]["value"] == 101325.0
    assert case.boundaries["outflow"]["field_conditions"]["temperature"]["value"] == 307.5
    assert panel.boundary_matrix.rowCount() == 2
    assert "BLOCKED" in panel.boundary_matrix.item(0, 3).text()
    monkeypatch.setattr(QMessageBox, "question", lambda *args, **kwargs: QMessageBox.StandardButton.Yes)
    panel.boundary_matrix.selectRow(1)
    panel.boundary_field_remove_button.click()
    assert case.boundaries["outflow"]["fields"] == ["pressure"]
    assert case.boundaries["outflow"]["value"] == 101325.0
    assert panel.boundary_matrix.rowCount() == 1
    panel.close()
    app.processEvents()


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="PySide6 optional")
def test_cockpit_does_not_display_console_as_verified_evidence(tmp_path):
    from PySide6.QtWidgets import QApplication
    from cfdx import CFDXSession, ExecutionController, SolverRunner
    from cfdx.application import Application
    from cfdx.convergence_evidence import baseline_signatures, load_run_evidence
    from cfdx.run_center_panel import RunCenterPanel

    app = QApplication.instance() or QApplication(["cfdx-lot4-evidence"])
    session = CFDXSession()
    controller = ExecutionController(session, SolverRunner([sys.executable, "-c", "pass"]))
    controller._output("Iteration 1 Time = 0.1 residual p = 1e-12", False)
    application = Application(session, controller=controller)
    panel = RunCenterPanel(application)
    assert not panel.export_evidence.isEnabled()
    assert "NOT VERIFIED" in panel.verdict.text()

    baseline = baseline_signatures(tmp_path)
    (tmp_path / "execution.json").write_text(json.dumps({
        "format": "CFDX-EXECUTION", "schema_version": 1, "process_exit_code": 0,
        "converged": True, "iterations": 1, "convergence_reason": "native solver",
        "convergence_status": 2, "artifacts": {"convergence_json": True},
    }), encoding="utf8")
    (tmp_path / "convergence.json").write_text(json.dumps({
        "format": "CFDX-CONVERGENCE", "schema_version": 1,
        "converged": True, "iterations": 1,
        "history": [{"iteration": 1, "continuity_normalized": 1e-9}],
    }), encoding="utf8")
    controller.convergence_evidence = load_run_evidence(
        tmp_path, process_exit_code=0, baseline=baseline,
    )
    controller.convergence_verdict = True
    panel._refresh(application.state)
    assert panel.export_evidence.isEnabled()
    assert panel.evidence_table.rowCount() == 1
    assert "Independent physical qualification remains required" in panel.verdict.text()
    panel.close()
    application.run_center.close()
    app.processEvents()


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="PySide6 optional")
def test_mesh_health_dashboard_actions_present_without_loading_mesh():
    from PySide6.QtWidgets import QApplication
    from cfdx.case_file_panel import CaseFilePanel

    app = QApplication.instance() or QApplication(["cfdx-lot4-mesh"])
    panel = CaseFilePanel()
    assert panel.volume_health_button.objectName() == "workbench.case_file.volume_health"
    assert panel.mesh_health_export.objectName() == "workbench.case_file.export_mesh_health"
    assert not panel.volume_health_button.isEnabled()
    assert not panel.mesh_health_export.isEnabled()
    panel.close()
    app.processEvents()

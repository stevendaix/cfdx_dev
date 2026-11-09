"""Qt/OpenGL-exercisable contracts for the native cell hotspot cockpit."""

from __future__ import annotations

import importlib.util
import sys

import pytest
from test_cell_convergence import fixture_trace


@pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 optional"
)
def test_hotspot_panel_source_selection_and_history(tmp_path):
    from cfdx.cell_hotspots_panel import CellHotspotsPanel
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication(["cfdx-hotspots-qt"])
    source = fixture_trace(tmp_path)
    panel = CellHotspotsPanel()
    panel.set_source(source)
    # The reader runs on the Qt worker pool. Exercise actual signal handoff.
    from time import monotonic, sleep

    deadline = monotonic() + 4
    while panel.trace is None and monotonic() < deadline:
        app.processEvents()
        sleep(0.01)
    assert panel.trace is not None
    assert panel.trace.complete
    assert panel.frame_choice.currentData() == 4
    assert panel.table.rowCount() == 1
    panel.table.selectRow(0)
    app.processEvents()
    assert panel.history_table.rowCount() == 1
    assert panel.history_table.item(0, 0).text() == "4"
    assert "Cell 2" in panel.history_note.text()
    panel.close()
    app.processEvents()


@pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 optional"
)
def test_run_checkbox_gates_native_opt_in_without_spawning(monkeypatch):
    from cfdx.run_center_panel import RunCenterPanel
    from PySide6.QtWidgets import QApplication

    from cfdx import CFDXSession, ExecutionController, SolverRunner
    from cfdx.application import Application

    app = QApplication.instance() or QApplication(["cfdx-hotspots-run"])
    session = CFDXSession()
    runner = SolverRunner(["cfdx_production_solver", "--mesh", "native.h5"])
    application = Application(session, controller=ExecutionController(session, runner))
    panel = RunCenterPanel(application)
    requests = []
    monkeypatch.setattr(panel, "_invoke", requests.append)
    panel.spatial_diagnostics.setChecked(True)
    panel.run_button.click()
    assert "--cell-hotspots" in runner.command
    assert requests == ["run"]
    panel.spatial_diagnostics.setChecked(False)
    panel.run_button.click()
    assert "--cell-hotspots" not in runner.command
    assert requests == ["run", "run"]
    runner.command = (sys.executable, "-c", "pass")
    panel.spatial_diagnostics.setChecked(True)
    panel.run_button.click()
    assert requests == ["run", "run"]
    assert "requires the native" in panel.status.text()
    panel.close()
    app.processEvents()

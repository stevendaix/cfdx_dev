"""Optional CI Qt tests for the actual connected Mission Control widgets."""
from __future__ import annotations

import importlib.util

import pytest


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="PySide6 optional")
def test_overview_live_checks_palette_and_window_lifecycle(tmp_path, monkeypatch):
    from PySide6.QtCore import QSettings
    from PySide6.QtWidgets import QApplication, QDockWidget
    from cfdx.session import CFDXSession
    from cfdx.workbench import CFDXWorkbenchWindow

    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    app = QApplication.instance() or QApplication(["cfdx-mission-control-qt"])
    settings = QSettings("CFDX", "Workbench")
    settings.clear()
    window = CFDXWorkbenchWindow(CFDXSession())
    assert window._theme == "black"
    assert window.workspace_actions["overview"].shortcut().toString() == "Ctrl+0"
    assert window.findChild(QDockWidget, "workbench.dock.overview") is not None
    window.workspace_actions["overview"].trigger()
    assert window._workspace == "overview"
    panel = window.overview_panel
    assert panel.snapshot is not None
    assert panel.snapshot.check_state == "not_checked"
    panel.check_button.click()
    app.processEvents()
    assert panel.snapshot.check_state in ("errors", "warnings")
    window.application.set_project_path(None, dirty=True)
    app.processEvents()
    assert panel.snapshot.check_state == "stale"
    window._show_command_palette()
    window.command_palette_query.setText("Overview")
    assert window.command_palette_results.count() == 1
    window.command_palette_dialog.close()
    bus = window.application.events
    assert sum(len(rows) for rows in bus._subscribers.values()) >= 3
    window.close()
    app.processEvents()
    assert all(callback not in bus._subscribers[event]
               for event, callback in window._event_subscriptions)
    settings.clear()

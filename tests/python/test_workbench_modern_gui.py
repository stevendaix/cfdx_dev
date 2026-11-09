"""CI graphical smoke; only executes on real Qt installations."""

import importlib.util

import pytest


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="Qt optional")
def test_methods_and_boundaries_are_navigable_without_solver_evidence():
    from cfdx.case import Case
    from cfdx.gui import create_application
    from cfdx.methods_panel import MethodsPanel
    from cfdx.setup_panel import CaseSetupPanel

    app = create_application(["cfdx-methods-smoke"])
    case = Case()
    setup = CaseSetupPanel(case)
    assert setup.boundary_search.objectName() == "workbench.boundaries.search"
    setup.boundary_preset.setCurrentText("Velocity inlet")
    setup.boundary_preset_button.click()
    assert setup.boundary_type.currentText() == "inlet"
    assert setup.boundary_field.currentText() == "velocity"
    assert case.boundaries == {}  # choosing a preset NEVER silently commits
    setup.boundary_name.setText("inflow")
    setup.boundary_x.setValue(10.0)
    setup.boundary_button.click()
    assert case.boundaries["inflow"]["velocity_value"][0] == 10.0
    assert setup.boundary_list.currentItem().text() == "inflow"
    setup.boundary_field.setCurrentText("temperature")
    assert setup.boundary_scalar_type.currentText() == "ZERO_GRADIENT"
    assert setup.boundary_value.value() == 0.0
    assert setup.boundary_gradient.value() == 0.0
    setup.boundary_field.setCurrentText("velocity")
    assert setup.boundary_scalar_type.currentText() == "FIXED_VALUE"
    assert setup.boundary_x.value() == 10.0
    setup.boundary_search.setText("other")
    assert setup.boundary_list.item(0).isHidden()
    setup.boundary_search.setText("")
    assert not setup.boundary_list.item(0).isHidden()
    methods = MethodsPanel()
    methods.family.setCurrentText("Nonlinear acceleration")
    methods.method.setCurrentIndex(0)
    assert "no native" in methods.details.toPlainText().lower()
    methods.close()
    setup.close()
    app.processEvents()


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="Qt optional")
def test_black_default_appearance_switch_and_layout_reset(tmp_path, monkeypatch):
    from cfdx.session import CFDXSession
    from cfdx.theme import THEMES
    from cfdx.workbench import CFDXWorkbenchWindow
    from PySide6.QtCore import QSettings
    from PySide6.QtGui import QPalette
    from PySide6.QtWidgets import QApplication, QPushButton

    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    app = QApplication.instance() or QApplication(["cfdx-appearance-smoke"])
    # Isolate preferences from any other test or developer's machine.
    settings = QSettings("CFDX", "Workbench")
    settings.remove("appearance/theme")
    settings.remove("workspace")
    window = CFDXWorkbenchWindow(CFDXSession())
    assert window._theme == "black"
    assert window.theme_actions["black"].isChecked()
    assert window.palette().color(QPalette.ColorRole.Window).name() == "#000000"
    assert window.findChild(QPushButton, "workbench.viewport.go.validate")
    assert window.workspace_actions["validate"].shortcut().toString() == "Ctrl+5"
    window._show_command_palette()
    window.command_palette_query.setText("Appearance: Light")
    assert window.command_palette_results.count() == 1
    window.command_palette_dialog.close()
    window.theme_actions["light"].trigger()
    assert window._theme == "light"
    assert settings.value("appearance/theme") == "light"
    assert (
        window.palette().color(QPalette.ColorRole.Window).name()
        == THEMES["light"]["surface.canvas"].lower()
    )
    window.theme_actions["graphite"].trigger()
    assert window._theme == "graphite"
    window._set_workspace("validate")
    window.reset_layout_action.trigger()
    assert window._workspace == "setup"
    # A reset changes the layout, not the explicitly chosen appearance.
    assert window._theme == "graphite"
    window.close()
    app.processEvents()
    assert settings.value("appearance/theme") == "graphite"
    settings.remove("appearance/theme")
    settings.remove("workspace")


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="Qt optional")
def test_methods_search_and_comparison_never_write_case():
    from cfdx.case import Case
    from cfdx.gui import create_application
    from cfdx.methods_panel import MethodsPanel

    app = create_application(["cfdx-methods-filter-test"])
    case = Case()
    before = repr(case)
    panel = MethodsPanel()
    panel.family.setCurrentText("Nonlinear acceleration")
    panel.integration.setCurrentIndex(panel.integration.findData("research"))
    panel.search.setText("Anderson")
    assert panel.method.count() == 1
    assert panel.method.currentData() == "ANDERSON"
    index = panel.compare_with.findData("NGMRES")
    panel.compare_with.setCurrentIndex(index)
    assert panel.compare_button.isEnabled()
    panel.compare_button.click()
    assert "no native integration certified" in panel.comparison.toPlainText()
    assert repr(case) == before
    panel.close()
    app.processEvents()

import importlib.util

import pytest
from cfdx.workbench import CFDXWorkbenchWindow

from cfdx import CFDXSession


def test_workbench_is_importable_without_qt() -> None:
    assert CFDXWorkbenchWindow is not None
    if importlib.util.find_spec("PySide6") is not None:
        pytest.skip("PySide6 is installed")
    with pytest.raises(RuntimeError, match="PySide6 is required"):
        CFDXWorkbenchWindow(CFDXSession())


@pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 optional"
)
def test_workbench_has_stable_docks() -> None:
    from cfdx.gui import create_application
    from cfdx.session import SimulationState
    from PySide6.QtWidgets import QDockWidget, QLabel, QPushButton

    app = create_application(["cfdx-workbench-test"])
    session = CFDXSession()
    session.state = SimulationState.READY
    window = CFDXWorkbenchWindow(session)
    from cfdx.application import WorkflowStatus

    assert window._status_marker(WorkflowStatus.CONDITIONAL) == "[~]"
    assert window._status_marker(WorkflowStatus.COMPLETE) == "[OK]"
    assert window._status_marker(WorkflowStatus.NOT_CONFIGURED) == "[ ]"
    assert window.objectName() == ""
    assert window.workflow_tree.objectName() == "workbench.workflow_tree"
    assert window.workflow_tree.topLevelItem(1).child(2).text(0).endswith("[OK]")
    assert window.workflow_tree.topLevelItem(0).child(4).text(0).endswith("[ ]")
    assert window.findChild(QDockWidget, "workbench.dock.workflow") is not None
    assert window.findChild(QDockWidget, "workbench.dock.setup") is not None
    assert window.findChild(QDockWidget, "workbench.dock.probes") is not None
    assert window.findChild(QDockWidget, "workbench.dock.properties") is not None
    assert window.findChild(QDockWidget, "workbench.dock.monitor") is not None
    assert window.centralWidget().objectName() == "workbench.viewport"
    assert window.viewport_stack.currentWidget() is window.viewport_placeholder
    assert (
        window.viewport_placeholder.findChild(
            QLabel, "workbench.viewport.empty_title"
        ).text()
        == "Your simulation workspace"
    )
    assert window.workflow_tree.topLevelItem(0).font(0).bold()
    assert window._workspace == "setup"
    assert not window.findChild(QDockWidget, "workbench.dock.mesh").isHidden()
    assert window.findChild(QDockWidget, "workbench.dock.results").isHidden()
    assert (
        window.findChild(QDockWidget, "workbench.dock.visualization_pipeline")
        is not None
    )
    assert (
        window.visualization_operations.objectName()
        == "workbench.visualization.operations"
    )
    assert window.visualization_operations.count() == 0
    assert not window.remove_visualization_operation_button.isEnabled()
    assert not window.reset_visualization_button.isEnabled()
    window.workspace_actions["analyze"].trigger()
    assert window._workspace == "analyze"
    assert not window.findChild(QDockWidget, "workbench.dock.results").isHidden()
    assert not window.findChild(QDockWidget, "workbench.dock.probes").isHidden()
    assert window.findChild(QDockWidget, "workbench.dock.mesh").isHidden()
    window.workspace_actions["setup"].trigger()
    assert (
        window.viewport_placeholder.findChild(
            QPushButton, "workbench.empty_state.open_mesh"
        ).text()
        == "Open Mesh"
    )
    assert (
        window.viewport_placeholder.findChild(
            QPushButton, "workbench.empty_state.open_dat_result"
        ).text()
        == "Load DAT Result"
    )
    assert (
        window.viewport_placeholder.findChild(
            QPushButton, "workbench.empty_state.open_results"
        ).text()
        == "Open Results Folder"
    )
    window.workflow_tree.setCurrentItem(window.workflow_tree.topLevelItem(0).child(4))
    assert window.search_action.isEnabled()
    assert window.search_action.shortcut().toString() == "Ctrl+K"
    window._show_command_palette()
    assert window.command_palette_dialog.objectName() == "workbench.command_palette"
    window.command_palette_query.setText("Analyze")
    assert window.command_palette_results.count() == 1
    assert window.command_palette_results.item(0).text() == "Analyze workspace"
    window.command_palette_dialog.close()
    assert window.application.state.selection.stable_id == "boundaries"
    window.close()
    app.quit()


@pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 optional"
)
def test_properties_inspector_applies_typed_cfl_value() -> None:
    from cfdx.gui import create_application
    from cfdx.session import SimulationState
    from PySide6.QtWidgets import QDoubleSpinBox, QPushButton

    app = create_application(["cfdx-workbench-test"])
    session = CFDXSession()
    session.state = SimulationState.READY
    window = CFDXWorkbenchWindow(session)
    window.application.select("numerics", "workflow", "Numerics")
    editor = window.findChild(QDoubleSpinBox, "workbench.property.numerics.cfl")
    assert editor is not None
    editor.setValue(2.5)
    window.findChild(QPushButton, "workbench.property.apply.numerics.cfl").click()
    assert session.case.numerics["cfl"] == 2.5
    assert isinstance(session.case.numerics["cfl"], float)
    window.close()
    app.quit()


@pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 optional"
)
def test_properties_inspector_applies_typed_physics_toggle() -> None:
    from cfdx.gui import create_application
    from cfdx.session import SimulationState
    from PySide6.QtWidgets import QCheckBox, QPushButton

    app = create_application(["cfdx-workbench-test"])
    session = CFDXSession()
    session.state = SimulationState.READY
    session.case.physics["energy"] = {"enabled": True}
    window = CFDXWorkbenchWindow(session)
    window.application.select("physics", "workflow", "Physics")
    toggle = window.findChild(QCheckBox, "workbench.property.physics.energy.enabled")
    assert toggle is not None
    assert toggle.isChecked()
    toggle.setChecked(False)
    window.findChild(
        QPushButton, "workbench.property.apply.physics.energy.enabled"
    ).click()
    assert session.case.physics["energy"]["enabled"] is False
    window.close()
    app.quit()


@pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 optional"
)
def test_mesh_selection_updates_typed_inspector_context() -> None:
    from cfdx.gui import create_application
    from cfdx.mesh_model import MeshCatalog, MeshPatch, MeshRegion
    from cfdx.session import SimulationState
    from PySide6.QtWidgets import QLineEdit, QSpinBox

    app = create_application(["cfdx-workbench-test"])
    session = CFDXSession()
    session.state = SimulationState.READY
    window = CFDXWorkbenchWindow(session)
    patch = MeshPatch(0, "inlet", "inlet", 2, "fluid")
    region = MeshRegion("fluid", 16, 19, 5, (patch,), ("c0", "c1"))
    catalog = MeshCatalog(16, 19, 5, (patch,), (region,))
    window.application.mesh_catalog = catalog
    window.mesh_browser.set_catalog(catalog)
    browser = window.mesh_browser.patches
    browser.setCurrentItem(browser.topLevelItem(0).child(0))
    assert window.application.selection.stable_id == "region:fluid/patch:inlet"
    patch_type = window.findChild(QLineEdit, "workbench.property.mesh.patch_type")
    assert patch_type is not None
    assert patch_type.text() == "inlet"
    assert not patch_type.isEnabled()
    faces = window.findChild(QSpinBox, "workbench.property.mesh.faces")
    assert faces is not None
    assert faces.value() == 2
    browser.setCurrentItem(browser.topLevelItem(0))
    assert window.application.selection.stable_id == "region:fluid"
    cells = window.findChild(QSpinBox, "workbench.property.mesh.cells")
    assert cells is not None
    assert cells.value() == 5
    window.close()
    app.quit()


@pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 optional"
)
def test_undo_and_redo_restore_previous_property_value() -> None:
    from cfdx.gui import create_application
    from cfdx.session import SimulationState
    from PySide6.QtCore import QCoreApplication, QEvent
    from PySide6.QtWidgets import QDoubleSpinBox, QPushButton

    app = create_application(["cfdx-workbench-test"])
    session = CFDXSession()
    session.state = SimulationState.READY
    window = CFDXWorkbenchWindow(session)
    window.application.select("numerics", "workflow", "Numerics")
    initial = session.case.numerics.get("cfl", 1.0)
    assert not window.undo_action.isEnabled()
    assert not window.redo_action.isEnabled()
    assert window.undo_action.statusTip() == "Requires: an edit to undo"
    editor = window.findChild(QDoubleSpinBox, "workbench.property.numerics.cfl")
    assert editor is not None
    editor.setValue(2.5)
    window.findChild(QPushButton, "workbench.property.apply.numerics.cfl").click()
    assert session.case.numerics["cfl"] == 2.5
    assert window.undo_action.isEnabled()
    assert not window.redo_action.isEnabled()
    window.undo_action.trigger()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    assert session.case.numerics.get("cfl", 1.0) == initial
    editor = window.findChild(QDoubleSpinBox, "workbench.property.numerics.cfl")
    assert editor is not None
    assert editor.value() == initial
    assert window.statusBar().currentMessage() == "Undid: set numerics.cfl to 2.5"
    assert not window.undo_action.isEnabled()
    assert window.redo_action.isEnabled()
    window.redo_action.trigger()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    assert session.case.numerics["cfl"] == 2.5
    editor = window.findChild(QDoubleSpinBox, "workbench.property.numerics.cfl")
    assert editor is not None
    assert editor.value() == 2.5
    assert window.statusBar().currentMessage() == "Redid: set numerics.cfl to 2.5"
    assert window.undo_action.isEnabled()
    assert not window.redo_action.isEnabled()
    window.close()
    app.quit()


@pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 optional"
)
def test_command_palette_lists_undo_with_history_prerequisite() -> None:
    from cfdx.gui import create_application

    app = create_application(["cfdx-workbench-test"])
    window = CFDXWorkbenchWindow(CFDXSession())
    window._show_command_palette()
    window.command_palette_query.setText("Undo")
    assert window.command_palette_results.count() == 1
    assert (
        window.command_palette_results.item(0).text()
        == "Undo — Requires: an edit to undo"
    )
    window.command_palette_dialog.close()
    window.close()
    app.quit()


@pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 optional"
)
def test_undo_and_redo_actions_have_stable_ids_and_shortcuts() -> None:
    from cfdx.gui import create_application

    app = create_application(["cfdx-workbench-test"])
    window = CFDXWorkbenchWindow(CFDXSession())
    assert window.undo_action.objectName() == "workbench.action.undo"
    assert window.redo_action.objectName() == "workbench.action.redo"
    assert window.undo_action.text() == "Undo"
    assert window.redo_action.text() == "Redo"
    assert [shortcut.toString() for shortcut in window.undo_action.shortcuts()] == [
        "Ctrl+Z"
    ]
    assert [shortcut.toString() for shortcut in window.redo_action.shortcuts()] == [
        "Ctrl+Shift+Z",
        "Ctrl+Y",
    ]
    window.close()
    app.quit()

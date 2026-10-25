"""Behavior tests for the capability registry and solver state views.

Covers gap G3 of docs/gui/PHASE0_AUDIT.md (no silently disabled action) and
the P0 solver-state review requirement (RUNNING is never convergence).
"""

import importlib.util

import pytest
from cfdx.capabilities import (
    CAPABILITIES,
    CapabilityContext,
    evaluate_capability,
    get_capability,
)
from cfdx.gui import create_application
from cfdx.solver_states import (
    convergence_status_view,
    session_state_view,
    unknown_state_view,
)
from cfdx.workbench import CFDXWorkbenchWindow

from cfdx import CFDXSession


def test_capability_ids_are_stable_and_unique() -> None:
    ids = list(CAPABILITIES)
    assert len(ids) == len(set(ids))
    assert "workbench.action.contour" in CAPABILITIES
    for capability in CAPABILITIES.values():
        assert capability.id
        assert capability.label
        assert capability.lot >= 1


def test_contour_reports_missing_prerequisites_in_clear_text() -> None:
    status = evaluate_capability("workbench.action.contour", CapabilityContext())
    assert not status.available
    expected = "requires a 3D renderer and a mesh or result loaded in the 3D viewport"
    assert status.reason == expected
    assert status.tip.startswith("Requires:")


def test_contour_available_once_renderer_and_dataset_are_ready() -> None:
    context = CapabilityContext(renderer_ready=True, dataset_loaded=True)
    status = evaluate_capability("workbench.action.contour", context)
    assert status.available
    assert status.missing == ()
    assert status.tip == ""


def test_run_unavailable_while_simulation_is_active() -> None:
    context = CapabilityContext(simulation_active=True)
    status = evaluate_capability("workbench.action.run", context)
    assert not status.available
    assert "no simulation currently running" in status.reason
    assert evaluate_capability("workbench.action.run", CapabilityContext()).available


def test_stop_requires_an_active_simulation() -> None:
    idle = CapabilityContext()
    assert not evaluate_capability("workbench.action.stop", idle).available
    active = CapabilityContext(simulation_active=True)
    assert evaluate_capability("workbench.action.stop", active).available
    assert not evaluate_capability("workbench.action.run", active).available


def test_undo_and_redo_report_their_history_prerequisites() -> None:
    empty = CapabilityContext()
    undo = evaluate_capability("workbench.action.undo", empty)
    redo = evaluate_capability("workbench.action.redo", empty)
    assert not undo.available
    assert not redo.available
    assert undo.tip == "Requires: an edit to undo"
    assert redo.tip == "Requires: an edit to redo"
    ready = CapabilityContext(undo_available=True, redo_available=True)
    assert evaluate_capability("workbench.action.undo", ready).available
    assert evaluate_capability("workbench.action.redo", ready).available


def test_unknown_capability_is_rejected_explicitly() -> None:
    import pytest

    assert get_capability("workbench.action.nope") is None
    with pytest.raises(KeyError, match="unknown capability"):
        evaluate_capability("workbench.action.nope", CapabilityContext())


def test_solver_states_are_never_conflated() -> None:
    keys = ("RUNNING", "CONVERGED", "STAGNATED", "DIVERGED", "STOPPED")
    markers = [convergence_status_view(key).marker for key in keys]
    assert len(set(markers)) == len(markers)

    running = session_state_view("RUNNING")
    assert not running.converged
    assert not running.terminal
    assert "convergence not established" in running.label
    assert session_state_view("CONVERGED").converged
    for key in ("STAGNATED", "DIVERGED", "MAX_ITERATIONS"):
        view = convergence_status_view(key)
        assert view.terminal
        assert not view.converged
        assert view.token in {"status.warning", "status.error"}


def test_unknown_solver_values_are_surfaced_verbatim() -> None:
    view = convergence_status_view("SOMETHING_NEW")
    assert view.key == "SOMETHING_NEW"
    assert view.label == "SOMETHING_NEW"
    assert view.marker == "?"
    assert not view.converged
    assert unknown_state_view(None).key == "None"


@pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 optional"
)
def test_disabled_actions_explain_their_prerequisites() -> None:
    app = create_application(["cfdx-capabilities-test"])
    session = CFDXSession()
    window = CFDXWorkbenchWindow(session)
    for name in ("contour_action", "slice_action", "glyph_action"):
        action = getattr(window, name)
        assert not action.isEnabled()
        assert action.statusTip().startswith("Requires:")
        assert "mesh or result" in action.statusTip()

    assert window.run_action.isEnabled()
    assert not window.stop_action.isEnabled()
    assert "Requires:" in window.stop_action.statusTip()

    window._show_command_palette()
    window.command_palette_query.setText("Contour")
    assert window.command_palette_results.count() == 1
    item_text = window.command_palette_results.item(0).text()
    assert item_text.startswith("Contour —")
    assert "Requires:" in item_text
    window.command_palette_dialog.close()

    window.workspace_actions["analyze"].trigger()
    assert window._workspace == "analyze"
    # Restore the default workspace before closing: closeEvent persists the
    # workspace and dock layout through QSettings, and later test modules
    # assert the fresh-launch layout of a brand new window.
    window.workspace_actions["setup"].trigger()
    window.close()
    app.quit()


def test_registry_covers_every_workbench_command() -> None:
    assert set(CAPABILITIES) == {
        "workbench.action.contour",
        "workbench.action.slice",
        "workbench.action.glyph",
        "workbench.action.run",
        "workbench.action.check",
        "workbench.action.save",
        "workbench.action.stop",
        "workbench.action.undo",
        "workbench.action.redo",
        "workbench.action.new_project",
        "workbench.action.open",
        "workbench.action.open_mesh",
        "workbench.action.open_dat_result",
        "workbench.action.open_dat_checkpoint",
        "workbench.workspace.setup",
        "workbench.workspace.mesh",
        "workbench.workspace.run",
        "workbench.workspace.analyze",
    }


def test_file_and_workspace_commands_are_always_available() -> None:
    context = CapabilityContext()
    for capability_id in (
        "workbench.action.new_project",
        "workbench.action.open",
        "workbench.action.open_mesh",
        "workbench.action.open_dat_result",
        "workbench.action.open_dat_checkpoint",
        "workbench.workspace.setup",
        "workbench.workspace.mesh",
        "workbench.workspace.run",
        "workbench.workspace.analyze",
    ):
        status = evaluate_capability(capability_id, context)
        assert status.available
        assert status.tip == ""


def test_capability_labels_are_unique() -> None:
    labels = [capability.label for capability in CAPABILITIES.values()]
    assert len(labels) == len(set(labels))


@pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 optional"
)
def test_palette_lists_the_registry_and_only_the_registry() -> None:
    from PySide6.QtGui import QAction

    app = create_application(["cfdx-palette-test"])
    window = CFDXWorkbenchWindow(CFDXSession())
    assert set(window._capability_actions) == set(CAPABILITIES)
    # Unnamed QActions belong to Qt itself (QDockWidget view toggles,
    # QToolBar builtin actions) and the command palette never lists them:
    # every named action must be either a registered capability or the
    # palette search action.
    named_actions = {
        action.objectName()
        for action in window.findChildren(QAction)
        if action.objectName()
    }
    assert named_actions - set(CAPABILITIES) == {"workbench.action.search"}

    window._show_command_palette()
    results = window.command_palette_results
    assert results.count() == len(CAPABILITIES)
    labels = [results.item(i).text().split(" — ")[0] for i in range(results.count())]
    expected_labels = {capability.label for capability in CAPABILITIES.values()}
    assert set(labels) == expected_labels
    assert labels == sorted(labels, key=str.casefold)

    contour_text = results.item(labels.index("Contour")).text()
    assert "Requires:" in contour_text

    window.command_palette_query.setText("Mesh workspace")
    assert results.count() == 1
    window.command_palette_query.returnPressed.emit()
    assert window._workspace == "mesh"
    # Restore the default workspace before closing (see
    # test_disabled_actions_explain_their_prerequisites).
    window.workspace_actions["setup"].trigger()
    window.command_palette_dialog.close()
    window.close()
    app.quit()

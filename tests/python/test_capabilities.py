"""Behavior tests for the capability registry and solver state views.

Covers gap G3 of docs/gui/PHASE0_AUDIT.md (no silently disabled action) and
the P0 solver-state review requirement (RUNNING is never convergence).
"""

from cfdx.capabilities import (
    CAPABILITIES,
    CapabilityContext,
    evaluate_capability,
    get_capability,
)
from cfdx.solver_states import (
    convergence_status_view,
    session_state_view,
    unknown_state_view,
)


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


def test_disabled_actions_explain_their_prerequisites() -> None:
    import importlib.util

    import pytest
    from cfdx import CFDXSession
    from cfdx.gui import create_application
    from cfdx.workbench import CFDXWorkbenchWindow

    if importlib.util.find_spec("PySide6") is None:
        pytest.skip("PySide6 optional")
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
    window.close()
    app.quit()

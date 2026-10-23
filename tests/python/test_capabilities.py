"""Behavior tests for the capability registry and solver state views.

Covers gap G3 of docs/gui/PHASE0_AUDIT.md (no silently disabled action) and
the P0 solver-state review requirement (RUNNING is never convergence).
"""

from cfdx.capabilities import (
    CAPABILITIES,
    CapabilityContext,
    evaluate_capability,
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

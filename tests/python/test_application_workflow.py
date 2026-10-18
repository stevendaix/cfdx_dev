"""Workflow state reflects the post-processing registry's real availability.

The Workbench renders the post-processing capability registry as workflow
nodes. A node is "complete" only when every capability it represents can
actually be produced from a real production dataset through the results
contract; "conditional" means the capability is declared and owned by a known
kernel but has no export path yet. This test pins that mapping so a future
refactor cannot render a conditional capability as complete.
"""

from __future__ import annotations

from cfdx.application.workflow import (
    WorkflowStatus,
    _post_processing_availability,
    build_workflow_state,
    workflow_children,
)
from cfdx.session import CFDXSession


def test_post_processing_availability_maps_registry_to_workflow_status() -> None:
    availability = _post_processing_availability()
    # probes is complete: the probe CSV contract is real and round-trips
    # against artifacts actually emitted by the production solver.
    assert availability["probes"] is WorkflowStatus.COMPLETE
    # fields is conditional: the production solver's per-iteration VTU frames
    # carry pressure only, so velocity, the velocity magnitude and the
    # mesh-quality metrics are conditional on a richer artifact or an export
    # path that has not been built.
    assert availability["fields"] is WorkflowStatus.CONDITIONAL


def test_build_workflow_state_renders_fields_as_conditional() -> None:
    state = build_workflow_state(CFDXSession())
    by_id = {step.id: step for step in state}

    assert by_id["workflow.results"].status is WorkflowStatus.CONDITIONAL
    assert by_id["fields"].status is WorkflowStatus.CONDITIONAL
    assert by_id["probes"].status is WorkflowStatus.COMPLETE


def test_build_workflow_state_propagates_diagnostics_over_availability() -> None:
    """A validation diagnostic overrides the registry-derived status.

    The registry says a capability is conditional; a validation diagnostic
    that targets the same node is authoritative for the current case and must
    override the registry-derived status, because the Workbench renders the
    case-specific problem, not the capability declaration.
    """

    class Diagnostic:
        node_id = "fields"
        path = "fields.U"
        message = "velocity field not present in result series"
        severity = "error"

    state = build_workflow_state(CFDXSession(), diagnostics=[Diagnostic()])
    by_id = {step.id: step for step in state}
    assert by_id["fields"].status is WorkflowStatus.ERROR
    assert by_id["fields"].message == "velocity field not present in result series"
    assert by_id["workflow.results"].status is WorkflowStatus.ERROR


def test_workflow_children_filters_by_parent() -> None:
    state = build_workflow_state(CFDXSession())
    setup_children = workflow_children(state, "workflow.setup")
    assert [step.id for step in setup_children] == [
        "geometry",
        "mesh",
        "physics",
        "materials",
        "boundaries",
        "numerics",
        "solver",
    ]
    assert all(step.parent == "workflow.setup" for step in setup_children)

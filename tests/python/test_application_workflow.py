from dataclasses import dataclass

from cfdx import CFDXSession, SimulationState
from cfdx.application import WorkflowStatus, build_workflow_state, workflow_children


@dataclass
class Diagnostic:
    node_id: str
    severity: str
    message: str


def test_workflow_state_has_stable_groups_and_children() -> None:
    session = CFDXSession()
    state = build_workflow_state(session)

    roots = workflow_children(state, None)
    assert [step.id for step in roots] == ["workflow.setup", "workflow.run", "workflow.results"]
    assert [step.id for step in workflow_children(state, "workflow.setup")] == [
        "geometry", "mesh", "physics", "materials", "boundaries", "numerics", "solver"
    ]
    assert all(step.status is WorkflowStatus.NOT_CONFIGURED for step in roots)


def test_workflow_statuses_come_from_session_and_diagnostics() -> None:
    session = CFDXSession()
    session.state = SimulationState.READY
    state = build_workflow_state(
        session,
        diagnostics=[Diagnostic("boundaries", "error", "inlet is incomplete")],
    )
    by_id = {step.id: step for step in state}

    assert by_id["run"].status is WorkflowStatus.COMPLETE
    assert by_id["boundaries"].status is WorkflowStatus.ERROR
    assert by_id["boundaries"].message == "inlet is incomplete"

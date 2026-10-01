"""Workflow navigation state shared by Workbench and non-Qt clients."""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Iterable

from ..session import CFDXSession, SimulationState


class WorkflowStatus(str, Enum):
    COMPLETE = "complete"
    WARNING = "warning"
    ERROR = "error"
    NOT_CONFIGURED = "not_configured"


@dataclass(frozen=True)
class WorkflowStepState:
    id: str
    label: str
    kind: str
    status: WorkflowStatus
    parent: str | None = None
    message: str | None = None


_WORKFLOW_GROUPS = {
    "setup": ("geometry", "mesh", "physics", "materials", "boundaries", "numerics", "solver"),
    "run": ("check", "initialize", "run", "monitor", "checkpoint"),
    "results": ("fields", "contours", "vectors", "slices", "probes", "reports"),
}
_LABELS = {
    "setup": "SETUP",
    "run": "RUN",
    "results": "RESULTS",
    "geometry": "Geometry",
    "mesh": "Mesh",
    "physics": "Physics",
    "materials": "Materials",
    "boundaries": "Boundaries",
    "numerics": "Numerics",
    "solver": "Solver",
    "check": "Check",
    "initialize": "Initialize",
    "run": "Run",
    "monitor": "Monitor",
    "checkpoint": "Checkpoint",
    "fields": "Fields",
    "contours": "Contours",
    "vectors": "Vectors",
    "slices": "Slices",
    "probes": "Probes",
    "reports": "Reports",
}
_GROUP_IDS = {"setup": "workflow.setup", "run": "workflow.run", "results": "workflow.results"}


def build_workflow_state(
    session: CFDXSession,
    *,
    diagnostics: Iterable[object] = (),
) -> tuple[WorkflowStepState, ...]:
    """Build navigation state from the session and validation diagnostics.

    Diagnostics may optionally expose ``node_id`` and ``message`` attributes;
    unknown diagnostic shapes are ignored rather than copied into widgets.
    """
    diagnostic_by_node: dict[str, tuple[WorkflowStatus, str]] = {}
    for diagnostic in diagnostics:
        node_id = getattr(diagnostic, "node_id", None)
        message = getattr(diagnostic, "message", None)
        severity = str(getattr(diagnostic, "severity", "warning")).lower()
        if node_id is None:
            continue
        status = WorkflowStatus.ERROR if severity in {"error", "critical"} else WorkflowStatus.WARNING
        diagnostic_by_node[str(node_id)] = (status, str(message or ""))

    run_configured = session.state is not SimulationState.CREATED
    steps: list[WorkflowStepState] = []
    for group, node_ids in _WORKFLOW_GROUPS.items():
        group_id = _GROUP_IDS[group]
        steps.append(WorkflowStepState(group_id, _LABELS[group], group, WorkflowStatus.NOT_CONFIGURED))
        for node_id in node_ids:
            if node_id == "run":
                status = WorkflowStatus.COMPLETE if run_configured else WorkflowStatus.NOT_CONFIGURED
            else:
                status = WorkflowStatus.NOT_CONFIGURED
            message = None
            if node_id in diagnostic_by_node:
                status, message = diagnostic_by_node[node_id]
            steps.append(WorkflowStepState(node_id, _LABELS[node_id], node_id, status, group_id, message))
    return tuple(steps)


def workflow_children(state: Iterable[WorkflowStepState], parent: str | None) -> tuple[WorkflowStepState, ...]:
    """Return stable children for a navigation tree without GUI dependencies."""
    return tuple(step for step in state if step.parent == parent)

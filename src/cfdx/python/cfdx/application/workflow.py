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
    "setup": "SETUP", "run": "RUN", "results": "RESULTS", "geometry": "Geometry",
    "mesh": "Mesh", "physics": "Physics", "materials": "Materials", "boundaries": "Boundaries",
    "numerics": "Numerics", "solver": "Solver", "check": "Check", "initialize": "Initialize",
    "run": "Run", "monitor": "Monitor", "checkpoint": "Checkpoint", "fields": "Fields",
    "contours": "Contours", "vectors": "Vectors", "slices": "Slices", "probes": "Probes",
    "reports": "Reports",
}
_GROUP_IDS = {"setup": "workflow.setup", "run": "workflow.run", "results": "workflow.results"}


def diagnostic_node_id(path: str) -> str | None:
    """Map shared validation paths to stable workflow nodes."""
    root = path.split(".", 1)[0] if path else ""
    return {
        "name": "geometry",
        "mesh": "mesh",
        "physics": "physics",
        "materials": "materials",
        "boundaries": "boundaries",
        "numerics": "numerics",
        "execution": "solver",
    }.get(root)


def build_workflow_state(
    session: CFDXSession,
    *,
    diagnostics: Iterable[object] = (),
) -> tuple[WorkflowStepState, ...]:
    diagnostic_by_node: dict[str, tuple[WorkflowStatus, str]] = {}
    for diagnostic in diagnostics:
        node_id = getattr(diagnostic, "node_id", None) or diagnostic_node_id(getattr(diagnostic, "path", ""))
        message = getattr(diagnostic, "message", None)
        severity = str(getattr(diagnostic, "severity", "warning")).lower()
        if node_id is None:
            continue
        status = WorkflowStatus.ERROR if severity in {"error", "critical"} else WorkflowStatus.WARNING
        diagnostic_by_node[str(node_id)] = (status, str(message or ""))

    configured = {
        "geometry": False,
        "mesh": False,
        "physics": bool(session.case.physics),
        "materials": bool(session.case.materials),
        "boundaries": bool(session.case.boundaries),
        "numerics": bool(session.case.numerics),
        "solver": bool(session.case.execution.solver),
    }
    steps: list[WorkflowStepState] = []
    for group, node_ids in _WORKFLOW_GROUPS.items():
        group_id = _GROUP_IDS[group]
        steps.append(WorkflowStepState(group_id, _LABELS[group], group, WorkflowStatus.NOT_CONFIGURED))
        for node_id in node_ids:
            if node_id in {"check", "initialize", "monitor", "checkpoint"}:
                status = WorkflowStatus.NOT_CONFIGURED
            elif node_id == "run":
                status = WorkflowStatus.COMPLETE if session.state is not SimulationState.CREATED else WorkflowStatus.NOT_CONFIGURED
            else:
                status = WorkflowStatus.COMPLETE if configured.get(node_id, False) else WorkflowStatus.NOT_CONFIGURED
            message = None
            if node_id in diagnostic_by_node:
                status, message = diagnostic_by_node[node_id]
            steps.append(WorkflowStepState(node_id, _LABELS[node_id], node_id, status, group_id, message))
        children = steps[-len(node_ids):]
        if any(child.status is WorkflowStatus.ERROR for child in children):
            group_status = WorkflowStatus.ERROR
        elif any(child.status is WorkflowStatus.WARNING for child in children):
            group_status = WorkflowStatus.WARNING
        elif children and all(child.status is WorkflowStatus.COMPLETE for child in children):
            group_status = WorkflowStatus.COMPLETE
        else:
            group_status = WorkflowStatus.NOT_CONFIGURED
        steps[-len(node_ids) - 1] = WorkflowStepState(
            group_id, _LABELS[group], group, group_status
        )
    return tuple(steps)


def workflow_children(state: Iterable[WorkflowStepState], parent: str | None) -> tuple[WorkflowStepState, ...]:
    return tuple(step for step in state if step.parent == parent)

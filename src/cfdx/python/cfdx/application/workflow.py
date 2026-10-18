"""Workflow navigation state shared by Workbench and non-Qt clients."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum

from ..session import CFDXSession, SimulationState
from .postprocessing import build_post_processing_registry


class WorkflowStatus(str, Enum):
    COMPLETE = "complete"
    WARNING = "warning"
    ERROR = "error"
    NOT_CONFIGURED = "not_configured"
    # A capability that the registry declares but for which no export path
    # exists yet. It is not "not_configured" (the quantity is real and owned by
    # a known kernel) and not "complete" (nothing can actually be produced
    # from the production dataset through the results contract). Surfacing it
    # separately keeps the Workbench from showing a green check for a path
    # that has not been built.
    CONDITIONAL = "conditional"


@dataclass(frozen=True)
class WorkflowStepState:
    id: str
    label: str
    kind: str
    status: WorkflowStatus
    parent: str | None = None
    message: str | None = None


_WORKFLOW_GROUPS = {
    "setup": (
        "geometry",
        "mesh",
        "physics",
        "materials",
        "boundaries",
        "numerics",
        "solver",
    ),
    "run": ("check", "initialize", "run", "monitor", "checkpoint"),
    "results": ("fields", "contours", "vectors", "slices", "probes", "reports"),
}
_GROUP_IDS = {
    "setup": "workflow.setup",
    "run": "workflow.run",
    "results": "workflow.results",
}
# Group labels and node labels are kept separate because the "run" key is
# shared by both the RUN group and the Run node, and they need different
# casing. A single dict would silently drop one of the two.
_GROUP_LABELS = {
    "setup": "SETUP",
    "run": "RUN",
    "results": "RESULTS",
}
_NODE_LABELS = {
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


def _post_processing_availability() -> dict[str, WorkflowStatus]:
    """Map the post-processing registry's availability to workflow status.

    The registry is capability declaration, not executable evidence: a
    spec's ``availability`` says whether the quantity can actually be produced
    from a real production dataset. The Workbench renders that as a workflow
    node, so the mapping is the single place a frontend learns what is
    genuinely available versus what is only declared.

    The registry's stable_ids are mapped to the Workbench's result-group node
    ids (``fields``, ``probes``). A node is "conditional" when at least one of
    the capabilities it groups is conditional, because the node cannot be
    reported complete until every capability it represents is actually
    producible from the production dataset through the results contract.
    """
    registry = build_post_processing_registry()
    node_for_spec: dict[str, str] = {
        "field:p": "fields",
        "field:U": "fields",
        "derived:u_magnitude": "fields",
        "mesh_quality:skewness": "fields",
        "mesh_quality:non_orthogonality": "fields",
        "probe:p": "probes",
        "probe:u_x": "probes",
    }
    per_node: dict[str, set[WorkflowStatus]] = {}
    for spec in registry.all():
        node = node_for_spec.get(spec.stable_id)
        if node is None:
            continue
        status = {
            "available": WorkflowStatus.COMPLETE,
            "conditional": WorkflowStatus.CONDITIONAL,
            "unavailable": WorkflowStatus.NOT_CONFIGURED,
        }[spec.availability]
        per_node.setdefault(node, set()).add(status)
    mapping: dict[str, WorkflowStatus] = {}
    for node, statuses in per_node.items():
        if WorkflowStatus.NOT_CONFIGURED in statuses:
            mapping[node] = WorkflowStatus.NOT_CONFIGURED
        elif WorkflowStatus.CONDITIONAL in statuses:
            mapping[node] = WorkflowStatus.CONDITIONAL
        else:
            mapping[node] = WorkflowStatus.COMPLETE
    return mapping


def build_workflow_state(
    session: CFDXSession,
    *,
    diagnostics: Iterable[object] = (),
) -> tuple[WorkflowStepState, ...]:
    diagnostic_by_node: dict[str, tuple[WorkflowStatus, str]] = {}
    for diagnostic in diagnostics:
        node_id = getattr(diagnostic, "node_id", None) or diagnostic_node_id(
            getattr(diagnostic, "path", "")
        )
        message = getattr(diagnostic, "message", None)
        severity = str(getattr(diagnostic, "severity", "warning")).lower()
        if node_id is None:
            continue
        status = (
            WorkflowStatus.ERROR
            if severity in {"error", "critical"}
            else WorkflowStatus.WARNING
        )
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
    availability = _post_processing_availability()
    steps: list[WorkflowStepState] = []
    for group, node_ids in _WORKFLOW_GROUPS.items():
        group_id = _GROUP_IDS[group]
        steps.append(
            WorkflowStepState(
                group_id, _GROUP_LABELS[group], group, WorkflowStatus.NOT_CONFIGURED
            )
        )
        for node_id in node_ids:
            if node_id in {"check", "initialize", "monitor", "checkpoint"}:
                status = WorkflowStatus.NOT_CONFIGURED
            elif node_id == "run":
                status = (
                    WorkflowStatus.COMPLETE
                    if session.state is not SimulationState.CREATED
                    else WorkflowStatus.NOT_CONFIGURED
                )
            elif node_id in availability:
                status = availability[node_id]
            else:
                status = (
                    WorkflowStatus.COMPLETE
                    if configured.get(node_id, False)
                    else WorkflowStatus.NOT_CONFIGURED
                )
            message = None
            if node_id in diagnostic_by_node:
                status, message = diagnostic_by_node[node_id]
            steps.append(
                WorkflowStepState(
                    node_id, _NODE_LABELS[node_id], node_id, status, group_id, message
                )
            )
        children = steps[-len(node_ids) :]
        if any(child.status is WorkflowStatus.ERROR for child in children):
            group_status = WorkflowStatus.ERROR
        elif any(child.status is WorkflowStatus.WARNING for child in children):
            group_status = WorkflowStatus.WARNING
        elif any(child.status is WorkflowStatus.CONDITIONAL for child in children):
            # A conditional capability is real but not yet producible from the
            # production dataset through the results contract, so the group
            # cannot be reported complete even if every declared node is.
            group_status = WorkflowStatus.CONDITIONAL
        elif children and all(
            child.status is WorkflowStatus.COMPLETE for child in children
        ):
            group_status = WorkflowStatus.COMPLETE
        else:
            group_status = WorkflowStatus.NOT_CONFIGURED
        steps[-len(node_ids) - 1] = WorkflowStepState(
            group_id, _GROUP_LABELS[group], group, group_status
        )
    return tuple(steps)


def workflow_children(
    state: Iterable[WorkflowStepState], parent: str | None
) -> tuple[WorkflowStepState, ...]:
    return tuple(step for step in state if step.parent == parent)

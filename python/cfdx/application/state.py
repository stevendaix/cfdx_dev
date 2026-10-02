"""Read-only application snapshots for GUI, TUI, and CLI clients."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, TYPE_CHECKING

from ..session import CFDXSession, SimulationState
from .results import ResultsState
from .workflow import WorkflowStepState, build_workflow_state

if TYPE_CHECKING:
    from ..execution import ExecutionController


@dataclass(frozen=True)
class ProjectState:
    name: str
    path: str | None
    dirty: bool


@dataclass(frozen=True)
class SelectionState:
    stable_id: str | None = None
    kind: str | None = None
    label: str | None = None


@dataclass(frozen=True)
class ExecutionState:
    iteration: int
    time: float
    latest_metrics: Any | None
    monitor_sample_count: int
    error: Any | None


@dataclass(frozen=True)
class ApplicationState:
    """Immutable rendering snapshot; CFDXSession remains the source of truth."""

    project: ProjectState
    simulation_state: SimulationState
    case_revision: int
    mesh_revision: int
    physics_revision: int
    numerics_revision: int
    requires_restart: bool
    execution: ExecutionState
    selection: SelectionState
    capabilities: tuple[str, ...]
    diagnostics: tuple[Any, ...]
    case_tree: tuple[Any, ...]
    workflow: tuple[WorkflowStepState, ...]
    results: ResultsState


def build_application_state(
    session: CFDXSession,
    *,
    project_path: str | Path | None = None,
    dirty: bool = False,
    selection: SelectionState | None = None,
    controller: ExecutionController | None = None,
    capabilities: Iterable[str] = (),
    diagnostics: Iterable[Any] = (),
    results: ResultsState | None = None,
) -> ApplicationState:
    """Build a snapshot without copying or mutating CFD domain state."""
    diagnostics_tuple = tuple(diagnostics)
    monitor_count = 0
    latest_metrics = None
    execution_error = None
    if controller is not None:
        latest_metrics = controller.latest_metrics
        execution_error = controller.error
        monitor_count = len(controller.monitor_series.samples)

    return ApplicationState(
        project=ProjectState(
            name=session.case.name,
            path=str(Path(project_path)) if project_path is not None else None,
            dirty=dirty,
        ),
        simulation_state=session.state,
        case_revision=session.case_revision,
        mesh_revision=session.mesh_revision,
        physics_revision=session.physics_revision,
        numerics_revision=session.numerics_revision,
        requires_restart=session.requires_restart,
        execution=ExecutionState(
            iteration=session.iteration,
            time=session.time,
            latest_metrics=latest_metrics,
            monitor_sample_count=monitor_count,
            error=execution_error,
        ),
        selection=selection or SelectionState(),
        capabilities=tuple(sorted(set(capabilities))),
        diagnostics=diagnostics_tuple,
        case_tree=tuple(session.case_tree()),
        workflow=build_workflow_state(session, diagnostics=diagnostics_tuple),
        results=results or ResultsState(),
    )

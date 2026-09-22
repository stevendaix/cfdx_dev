"""Headless application session shared by future CLI/TUI/GUI frontends.

The session owns orchestration state only. Numerical kernels remain in C++.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

from .case import Case


class SimulationState(str, Enum):
    CREATED = "CREATED"
    VALIDATING = "VALIDATING"
    READY = "READY"
    RUNNING = "RUNNING"
    PAUSED = "PAUSED"
    STOPPED = "STOPPED"
    CONVERGED = "CONVERGED"
    FAILED = "FAILED"


class ChangeImpact(str, Enum):
    HOT = "HOT"
    RESTART = "RESTART"
    REBUILD = "REBUILD"


@dataclass(frozen=True)
class CaseNode:
    """Stable logical node in the application case tree."""

    id: str
    label: str
    kind: str
    parent: str | None = None


@dataclass(frozen=True)
class Checkpoint:
    case_revision: int
    mesh_revision: int
    physics_revision: int
    numerics_revision: int
    iteration: int
    time: float


@dataclass
class CFDXSession:
    """Single source of truth for headless/TUI/GUI application state."""

    case: Case = field(default_factory=Case)
    state: SimulationState = SimulationState.CREATED
    case_revision: int = 0
    mesh_revision: int = 0
    physics_revision: int = 0
    numerics_revision: int = 0
    iteration: int = 0
    time: float = 0.0
    _requires_restart: bool = False

    def validate(self) -> None:
        self.state = SimulationState.VALIDATING
        if self.case.execution.mpi_ranks < 1:
            self.state = SimulationState.FAILED
            raise ValueError("mpi_ranks must be positive")
        self.state = SimulationState.READY

    def run(self) -> None:
        if self.state is SimulationState.CREATED:
            self.validate()
        if self.state not in {
            SimulationState.READY,
            SimulationState.PAUSED,
            SimulationState.STOPPED,
        }:
            raise RuntimeError(f"cannot run from state {self.state.value}")
        self.state = SimulationState.RUNNING

    def pause(self) -> None:
        if self.state is SimulationState.RUNNING:
            self.state = SimulationState.PAUSED

    def stop(self) -> None:
        if self.state in {SimulationState.RUNNING, SimulationState.PAUSED}:
            self.state = SimulationState.STOPPED

    def converge(self) -> None:
        if self.state is not SimulationState.RUNNING:
            raise RuntimeError("convergence requires a running session")
        self.state = SimulationState.CONVERGED

    def edit(self, key: str, value: Any, impact: ChangeImpact = ChangeImpact.HOT) -> ChangeImpact:
        if self.state in {SimulationState.RUNNING, SimulationState.VALIDATING}:
            raise RuntimeError("case edits are not allowed while the session is active")
        self.case.numerics[key] = value
        self.case_revision += 1
        if impact is not ChangeImpact.HOT:
            self._requires_restart = True
        return impact

    @property
    def requires_restart(self) -> bool:
        return self._requires_restart

    def acknowledge_restart(self) -> None:
        self._requires_restart = False
        self.numerics_revision += 1

    def advance(self, iterations: int = 1, dt: float = 0.0) -> None:
        if self.state is not SimulationState.RUNNING:
            raise RuntimeError("advance requires a running session")
        if iterations < 0:
            raise ValueError("iterations must be non-negative")
        self.iteration += iterations
        self.time += dt * iterations

    def checkpoint(self) -> Checkpoint:
        return Checkpoint(
            case_revision=self.case_revision,
            mesh_revision=self.mesh_revision,
            physics_revision=self.physics_revision,
            numerics_revision=self.numerics_revision,
            iteration=self.iteration,
            time=self.time,
        )

    def case_tree(self) -> tuple[CaseNode, ...]:
        return (
            CaseNode("geometry", "Geometry", "geometry"),
            CaseNode("mesh", "Mesh", "mesh"),
            CaseNode("physics", "Physics", "physics"),
            CaseNode("materials", "Materials", "materials"),
            CaseNode("boundaries", "Boundaries", "boundaries"),
            CaseNode("numerics", "Numerics", "numerics"),
            CaseNode("solver", "Solver", "solver"),
            CaseNode("run", "Run", "run"),
            CaseNode("monitors", "Monitors", "monitors"),
            CaseNode("results", "Results", "results"),
            CaseNode("reports", "Reports", "reports"),
        )

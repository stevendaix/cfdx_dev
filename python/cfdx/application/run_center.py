"""Headless Run Center state fed by the existing execution controller."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from ..execution import ExecutionController
from ..metrics import SolverMetrics
from ..runner import ProcessResult
from ..session import SimulationState


@dataclass(frozen=True)
class RunCenterState:
    simulation_state: SimulationState
    iteration: int
    time: float
    latest_metrics: SolverMetrics | None
    output: tuple[str, ...]
    error: object | None


class RunCenterModel:
    """Presentation state for output, metrics and run controls.

    The controller remains the sole owner of process execution. This model only
    subscribes to its existing callbacks and keeps a bounded output transcript.
    """

    def __init__(self, controller: ExecutionController, *, max_output_lines: int = 2000) -> None:
        self.controller = controller
        self.max_output_lines = max_output_lines
        self.output: list[str] = []
        self._latest_metrics = getattr(controller, "latest_metrics", None)
        self.on_change: Callable[[RunCenterState], None] | None = None
        self._previous_output = getattr(controller, "on_output", None)
        self._previous_metrics = getattr(controller, "on_metrics", None)
        self._previous_complete = getattr(controller, "on_complete", None)
        controller.on_output = self._on_output
        controller.on_metrics = self._on_metrics
        controller.on_complete = self._on_complete

    @property
    def state(self) -> RunCenterState:
        session = self.controller.session
        return RunCenterState(
            session.state,
            session.iteration,
            session.time,
            self._latest_metrics,
            tuple(self.output),
            getattr(self.controller, "error", None),
        )

    def clear(self) -> RunCenterState:
        self.output.clear()
        self._changed()
        return self.state

    def close(self) -> None:
        """Detach without taking ownership of the controller callbacks."""
        self.controller.on_output = self._previous_output
        self.controller.on_metrics = self._previous_metrics
        self.controller.on_complete = self._previous_complete

    def _on_output(self, line: str, is_stderr: bool) -> None:
        self.output.append(f"STDERR | {line}" if is_stderr else line)
        if len(self.output) > self.max_output_lines:
            del self.output[: len(self.output) - self.max_output_lines]
        if self._previous_output:
            self._previous_output(line, is_stderr)
        self._changed()

    def _on_metrics(self, metrics: SolverMetrics) -> None:
        self._latest_metrics = metrics
        if self._previous_metrics:
            self._previous_metrics(metrics)
        self._changed()

    def _on_complete(self, result: ProcessResult) -> None:
        if self._previous_complete:
            self._previous_complete(result)
        self._changed()

    def _changed(self) -> None:
        if self.on_change:
            self.on_change(self.state)

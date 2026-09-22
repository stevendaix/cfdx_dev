"""Execution controller connecting the session, runner and progress parser."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from .metrics import SolverMetrics, SolverMetricsParser
from .runner import ProcessResult, SolverRunner
from .session import CFDXSession, SimulationState


@dataclass(frozen=True)
class ExecutionError:
    returncode: int
    command: tuple[str, ...]


class ExecutionController:
    """Drive solver execution without blocking the caller."""

    def __init__(self, session: CFDXSession, runner: SolverRunner) -> None:
        self.session = session
        self.runner = runner
        self.parser = SolverMetricsParser()
        self.latest_metrics: SolverMetrics | None = None
        self.error: ExecutionError | None = None
        self.on_output: Callable[[str, bool], None] | None = None
        self.on_metrics: Callable[[SolverMetrics], None] | None = None
        self.on_complete: Callable[[ProcessResult], None] | None = None

    def start(self) -> None:
        self.session.run()
        self.error = None
        self.latest_metrics = None
        self.runner.start(self._output, self._complete)

    def _output(self, line: str, is_stderr: bool) -> None:
        metrics = self.parser.parse(line)
        if metrics is not None:
            self.latest_metrics = metrics
            if metrics.iteration is not None:
                self.session.iteration = metrics.iteration
            if metrics.time is not None:
                self.session.time = metrics.time
            if self.on_metrics:
                self.on_metrics(metrics)
        if self.on_output:
            self.on_output(line, is_stderr)

    def _complete(self, result: ProcessResult) -> None:
        if result.returncode == 0:
            self.session.state = SimulationState.CONVERGED
        else:
            self.error = ExecutionError(result.returncode, result.command)
            self.session.state = SimulationState.FAILED
        if self.on_complete:
            self.on_complete(result)

    def stop(self, timeout: float = 5.0) -> None:
        self.runner.stop(timeout)
        if self.session.state is SimulationState.RUNNING:
            self.session.stop()

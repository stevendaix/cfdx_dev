""""Execution controller connecting the session, runner and progress parser."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from threading import Lock
from typing import Callable

from .metrics import SolverMetrics, SolverMetricsParser
from .monitors import MonitorSample, MonitorSeries
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
        self.monitor_series = MonitorSeries("solver", [])
        self.error: ExecutionError | None = None
        self.on_output: Callable[[str, bool], None] | None = None
        self.on_metrics: Callable[[SolverMetrics], None] | None = None
        self.on_complete: Callable[[ProcessResult], None] | None = None
        self._stop_requested = False
        self._monitor_lock = Lock()

    def restart(self, dat_path: str | Path) -> None:
        """Restart through the runner using the case-configured checkpoint option."""
        path = Path(dat_path)
        if not path.is_file():
            raise FileNotFoundError(path)
        if self.session.state in {SimulationState.RUNNING, SimulationState.PAUSED}:
            raise RuntimeError("cannot restart an active session")
        restart_option = self.session.case.execution.restart_option
        if restart_option is None:
            raise ValueError("A DAT restart is loaded but no solver restart option is configured")
        self.session.run()
        self.error = None
        self.latest_metrics = None
        with self._monitor_lock:
            self.monitor_series = MonitorSeries("solver", [])
        self._stop_requested = False
        self.runner.start(self._output, self._complete, restart_path=path, restart_option=restart_option)

    def start(self) -> None:
        self.session.run()
        self.error = None
        self.latest_metrics = None
        with self._monitor_lock:
            self.monitor_series = MonitorSeries("solver", [])
        self._stop_requested = False
        self.runner.start(self._output, self._complete)

    def _output(self, line: str, is_stderr: bool) -> None:
        metrics = self.parser.parse(line)
        if metrics is not None:
            self.latest_metrics = metrics
            if metrics.iteration is not None:
                self.session.iteration = metrics.iteration
            if metrics.time is not None:
                self.session.time = metrics.time
            if metrics.iteration is not None or metrics.time is not None or metrics.residuals:
                with self._monitor_lock:
                    iteration = metrics.iteration if metrics.iteration is not None else (
                        self.monitor_series.samples[-1].iteration
                        if self.monitor_series.samples else 0
                    )
                    time_value = metrics.time if metrics.time is not None else (
                        self.monitor_series.samples[-1].time
                        if self.monitor_series.samples else 0.0
                    )
                    values = {name: value for name, value in metrics.residuals}
                    if metrics.cfl is not None:
                        values["CFL"] = metrics.cfl
                    self.monitor_series.upsert(MonitorSample(iteration, time_value, values))
            if self.on_metrics:
                self.on_metrics(metrics)
        if self.on_output:
            self.on_output(line, is_stderr)

    def _complete(self, result: ProcessResult) -> None:
        if self._stop_requested:
            self.session.state = SimulationState.STOPPED
        elif result.returncode == 0:
            self.session.state = SimulationState.CONVERGED
        else:
            self.error = ExecutionError(result.returncode, result.command)
            self.session.state = SimulationState.FAILED
        if self.on_complete:
            self.on_complete(result)

    def pause(self) -> None:
        if self.session.state is not SimulationState.RUNNING:
            raise RuntimeError("pause requires a running session")
        self.runner.pause()
        self.session.pause()

    def resume(self) -> None:
        if self.session.state is not SimulationState.PAUSED:
            raise RuntimeError("resume requires a paused session")
        self.runner.resume()
        self.session.run()

    def stop(self, timeout: float = 5.0) -> None:
        if self.session.state not in {SimulationState.RUNNING, SimulationState.PAUSED}:
            return
        self._stop_requested = True
        self.runner.stop(timeout)
        self.session.stop()

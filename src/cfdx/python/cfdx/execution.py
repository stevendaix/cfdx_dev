"""Execution controller connecting the session, runner and progress parser."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from threading import Lock

from .boundary_setup import validate_boundary_definition
from .convergence_evidence import (
    ConvergenceEvidence,
    EvidenceError,
    baseline_signatures,
    load_run_evidence,
)
from .metrics import SolverMetrics, SolverMetricsParser
from .monitors import MonitorSample, MonitorSeries
from .probe_stream import ProbeStreamSample, parse_probe_stream_line
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
        self.convergence_verdict: bool | None = None
        self.convergence_evidence: ConvergenceEvidence | None = None
        self.convergence_evidence_error: str | None = None
        self._run_output_dir: Path | None = None
        self._run_evidence_baseline: dict = {}
        self.on_output: Callable[[str, bool], None] | None = None
        self.on_metrics: Callable[[SolverMetrics], None] | None = None
        self.on_probe_sample: Callable[[ProbeStreamSample], None] | None = None
        self.on_complete: Callable[[ProcessResult], None] | None = None
        self.latest_probe_samples: dict[str, ProbeStreamSample] = {}
        self.probe_stream_error: str | None = None
        self._stop_requested = False
        self._monitor_lock = Lock()

    def _prepare_run_evidence(self) -> Path | None:
        """Resolve the actual configured output folder before subprocess launch.

        The CLI may specify --output-dir even when the case model does not.
        Do not guess a directory from the current working directory.
        """
        candidate = self.session.case.execution.output_dir
        if not candidate:
            command = tuple(self.runner.command)
            for index, arg in enumerate(command):
                if arg == "--output-dir" and index + 1 < len(command):
                    candidate = command[index + 1]
                elif arg.startswith("--output-dir="):
                    candidate = arg.split("=", 1)[1]
        directory = Path(candidate) if candidate else None
        if directory is not None and not directory.is_absolute() and self.runner.cwd:
            directory = Path(self.runner.cwd) / directory
        self._run_output_dir = directory
        self._run_evidence_baseline = baseline_signatures(directory)
        self.convergence_evidence = None
        self.convergence_evidence_error = None
        return directory

    def _assert_native_boundary_ready(self) -> None:
        """Block unrepresentable staged fields even when Run bypasses Validate.

        The editable Python case and C++ CaseSetup JSON are distinct formats.
        This check prevents known loss of multiple fields on direct Run/Restart;
        it does not claim that a native bridge exists for legacy boundaries.
        """
        for name, boundary in self.session.case.boundaries.items():
            try:
                validate_boundary_definition(boundary)
            except (ValueError, KeyError, TypeError) as exc:
                raise ValueError(
                    f"boundary {name!r} cannot be sent to the native solver: {exc}"
                ) from exc

    def restart(self, dat_path: str | Path) -> None:
        """Restart through the runner using the case-configured checkpoint option."""
        path = Path(dat_path)
        if not path.is_file():
            raise FileNotFoundError(path)
        if self.session.state in {SimulationState.RUNNING, SimulationState.PAUSED}:
            raise RuntimeError("cannot restart an active session")
        restart_option = self.session.case.execution.restart_option
        if restart_option is None:
            raise ValueError(
                "A DAT restart is loaded but no solver restart option is configured"
            )
        self._assert_native_boundary_ready()
        output_dir = self._prepare_run_evidence()
        self.session.run()
        self.error = None
        self.convergence_verdict = None
        self.latest_metrics = None
        self.latest_probe_samples.clear()
        self.probe_stream_error = None
        with self._monitor_lock:
            self.monitor_series = MonitorSeries("solver", [])
        self._stop_requested = False
        try:
            self.runner.start(
                self._output,
                self._complete,
                restart_path=path,
                restart_option=restart_option,
                output_dir=output_dir,
            )
        except (OSError, RuntimeError, ValueError):
            self.session.state = SimulationState.FAILED
            raise

    def start(self) -> None:
        self._assert_native_boundary_ready()
        output_dir = self._prepare_run_evidence()
        self.session.run()
        self.error = None
        self.convergence_verdict = None
        self.latest_metrics = None
        self.latest_probe_samples.clear()
        self.probe_stream_error = None
        with self._monitor_lock:
            self.monitor_series = MonitorSeries("solver", [])
        self._stop_requested = False
        try:
            self.runner.start(
                self._output,
                self._complete,
                output_dir=output_dir,
            )
        except (OSError, RuntimeError, ValueError):
            self.session.state = SimulationState.FAILED
            raise

    def _output(self, line: str, is_stderr: bool) -> None:
        is_probe_stream_line = not is_stderr and (
            line == "CFDX_PROBE_SAMPLE_V1" or line.startswith("CFDX_PROBE_SAMPLE_V1\t")
        )
        stream_warning: str | None = None
        if is_probe_stream_line:
            try:
                sample = parse_probe_stream_line(line)
                assert sample is not None
                previous = self.latest_probe_samples.get(sample.name)
                if previous is not None and (
                    sample.iteration <= previous.iteration
                    or sample.time < previous.time
                ):
                    raise ValueError(f"probe {sample.name!r} samples are not monotonic")
            except ValueError as exc:
                self.probe_stream_error = str(exc)
                stream_warning = f"Probe stream warning: {exc}"
            else:
                self.latest_probe_samples[sample.name] = sample
                self.probe_stream_error = None
                if self.on_probe_sample:
                    self.on_probe_sample(sample)

        metrics = self.parser.parse(line)
        if metrics is not None:
            self.latest_metrics = metrics
            if metrics.iteration is not None:
                self.session.iteration = metrics.iteration
            if metrics.time is not None:
                self.session.time = metrics.time
            if (
                metrics.iteration is not None
                or metrics.time is not None
                or metrics.cfl is not None
                or metrics.residuals
            ):
                with self._monitor_lock:
                    iteration = (
                        metrics.iteration
                        if metrics.iteration is not None
                        else (
                            self.monitor_series.samples[-1].iteration
                            if self.monitor_series.samples
                            else 0
                        )
                    )
                    time_value = (
                        metrics.time
                        if metrics.time is not None
                        else (
                            self.monitor_series.samples[-1].time
                            if self.monitor_series.samples
                            else 0.0
                        )
                    )
                    values = {name: value for name, value in metrics.residuals}
                    if metrics.cfl is not None:
                        values["CFL"] = metrics.cfl
                    self.monitor_series.upsert(
                        MonitorSample(iteration, time_value, values)
                    )
            if self.on_metrics:
                self.on_metrics(metrics)
        if self.on_output and stream_warning:
            self.on_output(stream_warning, True)
        if self.on_output and not is_probe_stream_line:
            self.on_output(line, is_stderr)

    def _complete(self, result: ProcessResult) -> None:
        """Never equate clean process exit with scientific convergence.

        Exit code 1 means nonconvergence only when the native execution JSON
        proves it. Otherwise unknown nonzero returns are execution errors.
        """
        self.convergence_evidence = None
        self.convergence_verdict = None
        self.convergence_evidence_error = None
        if self._stop_requested:
            self.session.state = SimulationState.STOPPED
        elif result.returncode in (0, 1) and result.output_dir is not None:
            try:
                certificate = load_run_evidence(
                    result.output_dir,
                    process_exit_code=result.returncode,
                    baseline=self._run_evidence_baseline,
                )
            except (EvidenceError, OSError, ValueError) as exc:
                self.convergence_evidence_error = str(exc)
            else:
                self.convergence_evidence = certificate
                self.convergence_verdict = certificate.overall_converged

            if self.convergence_verdict is True:
                self.session.state = SimulationState.CONVERGED
            elif self.convergence_verdict is False:
                self.session.state = SimulationState.STOPPED
            elif result.returncode == 0:
                # Legacy or invalid artifact: no certified convergence.
                self.session.state = SimulationState.STOPPED
            else:
                self.error = ExecutionError(result.returncode, result.command)
                self.session.state = SimulationState.FAILED
        elif result.returncode == 0:
            self.convergence_evidence_error = "no output directory attached to this run"
            self.session.state = SimulationState.STOPPED
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

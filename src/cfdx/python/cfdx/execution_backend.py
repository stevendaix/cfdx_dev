"""Execution backend contracts for local and scheduler-based runs."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, Sequence

from .runner import SolverRunner


class ExecutionBackend(Protocol):
    def start(self) -> None: ...
    def stop(self, timeout: float = 5.0) -> None: ...


class LocalExecutionBackend:
    def __init__(self, command: Sequence[str], cwd: Path | None = None) -> None:
        self.runner = SolverRunner(command, cwd)

    def start(self) -> None:
        self.runner.start()

    def stop(self, timeout: float = 5.0) -> None:
        self.runner.stop(timeout)


@dataclass(frozen=True)
class SlurmConfig:
    partition: str | None = None
    nodes: int = 1
    tasks: int = 1
    time_limit: str | None = None


def slurm_command(command: Sequence[str], config: SlurmConfig) -> tuple[str, ...]:
    """Build an sbatch command without submitting it."""
    if not command:
        raise ValueError("command must not be empty")
    result = ["sbatch", "--wait", "--nodes", str(config.nodes), "--ntasks", str(config.tasks)]
    if config.partition:
        result += ["--partition", config.partition]
    if config.time_limit:
        result += ["--time", config.time_limit]
    result += ["--wrap", " ".join(str(item) for item in command)]
    return tuple(result)

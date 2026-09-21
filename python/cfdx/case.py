from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import subprocess
from typing import Any


@dataclass
class ExecutionConfig:
    policy: str = "AUTO"
    solver: str | None = None
    mpi_ranks: int = 1
    deterministic: bool = True


@dataclass
class Case:
    """Python-side case orchestration; numerical kernels remain in C++."""

    name: str = "untitled"
    physics: dict[str, Any] = field(default_factory=dict)
    numerics: dict[str, Any] = field(default_factory=dict)
    boundaries: dict[str, Any] = field(default_factory=dict)
    execution: ExecutionConfig = field(default_factory=ExecutionConfig)

    def enable(self, model: str) -> "Case":
        self.physics[model] = {"enabled": True}
        return self

    def set_numerics(self, **values: Any) -> "Case":
        self.numerics.update(values)
        return self

    def set_boundary(self, name: str, **values: Any) -> "Case":
        self.boundaries[name] = dict(values)
        return self

    def as_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "physics": self.physics,
            "numerics": self.numerics,
            "boundaries": self.boundaries,
            "execution": {
                "policy": self.execution.policy,
                "solver": self.execution.solver,
                "mpi_ranks": self.execution.mpi_ranks,
                "deterministic": self.execution.deterministic,
            },
        }

    def run(self, case_file: Path, extra_args: list[str] | None = None) -> int:
        if not self.execution.solver:
            raise ValueError("execution.solver must be set before run()")
        command = [self.execution.solver, str(case_file)]
        if extra_args:
            command.extend(extra_args)
        if self.execution.mpi_ranks > 1:
            command = [
                "mpiexec", "-n", str(self.execution.mpi_ranks), *command
            ]
        return subprocess.call(command)

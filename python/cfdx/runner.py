"""Non-blocking local solver runner for CFDX orchestration."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import subprocess
import threading
from typing import Callable, Sequence


@dataclass(frozen=True)
class ProcessResult:
    returncode: int
    command: tuple[str, ...]


class SolverRunner:
    """Own a solver subprocess without blocking the application thread."""

    def __init__(self, command: Sequence[str], cwd: Path | None = None) -> None:
        if not command:
            raise ValueError("command must not be empty")
        self.command = tuple(str(item) for item in command)
        self.cwd = cwd
        self._process: subprocess.Popen[str] | None = None
        self._thread: threading.Thread | None = None

    @property
    def running(self) -> bool:
        return self._process is not None and self._process.poll() is None

    def start(
        self,
        on_output: Callable[[str, bool], None] | None = None,
        on_complete: Callable[[ProcessResult], None] | None = None,
    ) -> None:
        if self.running:
            raise RuntimeError("solver is already running")

        self._process = subprocess.Popen(
            self.command,
            cwd=self.cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        def consume() -> None:
            assert self._process is not None
            assert self._process.stdout is not None
            for line in self._process.stdout:
                if on_output:
                    on_output(line.rstrip("\n"), False)
            returncode = self._process.wait()
            if on_complete:
                on_complete(ProcessResult(returncode, self.command))

        self._thread = threading.Thread(target=consume, name="cfdx-solver-output", daemon=True)
        self._thread.start()

    def stop(self, timeout: float = 5.0) -> None:
        if not self._process or self._process.poll() is not None:
            return
        self._process.terminate()
        try:
            self._process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            self._process.kill()
            self._process.wait()

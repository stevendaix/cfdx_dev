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


OutputCallback = Callable[[str, bool], None]
CompletionCallback = Callable[[ProcessResult], None]


class SolverRunner:
    """Own a local solver subprocess without blocking the application thread.

    stdout and stderr are consumed independently. The second callback argument
    identifies stderr, allowing a TUI/GUI to preserve stream semantics.
    """

    def __init__(self, command: Sequence[str], cwd: Path | None = None) -> None:
        if not command:
            raise ValueError("command must not be empty")
        self.command = tuple(str(item) for item in command)
        self.cwd = cwd
        self._process: subprocess.Popen[str] | None = None
        self._thread: threading.Thread | None = None
        self._lock = threading.Lock()

    @property
    def running(self) -> bool:
        process = self._process
        return process is not None and process.poll() is None

    def start(
        self,
        on_output: OutputCallback | None = None,
        on_complete: CompletionCallback | None = None,
    ) -> None:
        with self._lock:
            if self.running:
                raise RuntimeError("solver is already running")
            self._process = subprocess.Popen(
                self.command,
                cwd=self.cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
            )
            process = self._process

        def consume(stream, is_stderr: bool) -> None:
            for line in iter(stream.readline, ""):
                if on_output:
                    on_output(line.rstrip("\r\n"), is_stderr)
            stream.close()

        def monitor() -> None:
            assert process.stdout is not None
            assert process.stderr is not None
            readers = [
                threading.Thread(
                    target=consume,
                    args=(process.stdout, False),
                    name="cfdx-solver-stdout",
                    daemon=True,
                ),
                threading.Thread(
                    target=consume,
                    args=(process.stderr, True),
                    name="cfdx-solver-stderr",
                    daemon=True,
                ),
            ]
            for reader in readers:
                reader.start()
            returncode = process.wait()
            for reader in readers:
                reader.join()
            if on_complete:
                on_complete(ProcessResult(returncode, self.command))

        self._thread = threading.Thread(
            target=monitor,
            name="cfdx-solver-monitor",
            daemon=True,
        )
        self._thread.start()

    def stop(self, timeout: float = 5.0) -> None:
        process = self._process
        if process is None or process.poll() is not None:
            return
        process.terminate()
        try:
            process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()
        thread = self._thread
        if thread is not None and thread is not threading.current_thread():
            thread.join(timeout=timeout)

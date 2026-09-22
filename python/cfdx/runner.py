"""Non-blocking local solver runner for CFDX orchestration."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import os
import signal
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

    stdout and stderr are consumed independently. On POSIX systems the solver
    is placed in its own process group so pause/resume/stop also apply to MPI
    children. Windows exposes the lifecycle contract but reports pause/resume
    as unsupported until a native process-suspension backend is added.
    """

    def __init__(self, command: Sequence[str], cwd: Path | None = None) -> None:
        if not command:
            raise ValueError("command must not be empty")
        self.command = tuple(str(item) for item in command)
        self.cwd = cwd
        self._process: subprocess.Popen[str] | None = None
        self._thread: threading.Thread | None = None
        self._lock = threading.Lock()
        self._paused = False

    @property
    def running(self) -> bool:
        process = self._process
        return process is not None and process.poll() is None

    @property
    def paused(self) -> bool:
        return self._paused and self.running

    @property
    def supports_pause(self) -> bool:
        return os.name == "posix"

    def start(
        self,
        on_output: OutputCallback | None = None,
        on_complete: CompletionCallback | None = None,
        restart_path: Path | None = None,
        restart_option: str | None = "--restart",
    ) -> None:
        if restart_path is not None and restart_option is None:
            raise ValueError("restart_option must be configured for a restart")
        command = (
            self.command
            if restart_path is None
            else self.command + (str(restart_option), str(restart_path))
        )
        with self._lock:
            if self.running:
                raise RuntimeError("solver is already running")
            kwargs = {}
            if os.name == "posix":
                kwargs["start_new_session"] = True
            self._process = subprocess.Popen(
                command,
                cwd=self.cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
                **kwargs,
            )
            self._paused = False
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
                threading.Thread(target=consume, args=(process.stdout, False), name="cfdx-solver-stdout", daemon=True),
                threading.Thread(target=consume, args=(process.stderr, True), name="cfdx-solver-stderr", daemon=True),
            ]
            for reader in readers:
                reader.start()
            returncode = process.wait()
            for reader in readers:
                reader.join()
            if on_complete:
                on_complete(ProcessResult(returncode, command))

        self._thread = threading.Thread(target=monitor, name="cfdx-solver-monitor", daemon=True)
        self._thread.start()

    def _send_group_signal(self, sig: int) -> None:
        process = self._process
        if process is None or process.poll() is not None:
            raise RuntimeError("solver is not running")
        if os.name != "posix":
            raise NotImplementedError("pause/resume is not supported on Windows")
        os.killpg(process.pid, sig)

    def pause(self) -> None:
        if not self.supports_pause:
            raise NotImplementedError("pause/resume is not supported on Windows")
        if not self.paused:
            self._send_group_signal(signal.SIGSTOP)
            self._paused = True

    def resume(self) -> None:
        if not self.supports_pause:
            raise NotImplementedError("pause/resume is not supported on Windows")
        if not self.paused:
            return
        self._send_group_signal(signal.SIGCONT)
        self._paused = False

    def stop(self, timeout: float = 5.0) -> None:
        process = self._process
        if process is None or process.poll() is not None:
            return
        if os.name == "posix":
            os.killpg(process.pid, signal.SIGTERM)
        else:
            process.terminate()
        try:
            process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            if os.name == "posix":
                os.killpg(process.pid, signal.SIGKILL)
            else:
                process.kill()
            process.wait()
        self._paused = False
        thread = self._thread
        if thread is not None and thread is not threading.current_thread():
            thread.join(timeout=timeout)

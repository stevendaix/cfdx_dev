from __future__ import annotations

import sys
import time

import pytest

from cfdx.runner import SolverRunner


def wait_for_completion(runner: SolverRunner, timeout: float = 5.0) -> None:
    deadline = time.monotonic() + timeout
    while runner.running and time.monotonic() < deadline:
        time.sleep(0.01)
    assert not runner.running


def test_runner_streams_stdout_and_stderr() -> None:
    command = [
        sys.executable,
        "-u",
        "-c",
        "import sys; print('out', flush=True); print('err', file=sys.stderr, flush=True)",
    ]
    output: list[tuple[str, bool]] = []
    completed = []

    runner = SolverRunner(command)
    runner.start(on_output=lambda line, is_stderr: output.append((line, is_stderr)),
                 on_complete=completed.append)
    wait_for_completion(runner)

    assert sorted(output) == [("err", True), ("out", False)]
    assert completed[0].returncode == 0


def test_runner_rejects_empty_command() -> None:
    with pytest.raises(ValueError):
        SolverRunner([])


def test_runner_rejects_concurrent_start() -> None:
    command = [sys.executable, "-u", "-c", "import time; time.sleep(30)"]
    runner = SolverRunner(command)
    runner.start()
    try:
        with pytest.raises(RuntimeError):
            runner.start()
    finally:
        runner.stop()


def test_runner_stop_terminates_process() -> None:
    command = [sys.executable, "-u", "-c", "import time; time.sleep(30)"]
    runner = SolverRunner(command)
    runner.start()
    deadline = time.monotonic() + 5
    while not runner.running and time.monotonic() < deadline:
        time.sleep(0.01)
    assert runner.running
    runner.stop()
    assert not runner.running

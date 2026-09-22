from __future__ import annotations

import sys
import time

from cfdx.runner import SolverRunner


def test_runner_streams_output_and_completes() -> None:
    command = [sys.executable, "-u", "-c", "print('step 1', flush=True); print('step 2', flush=True)"]
    output: list[str] = []
    completed = []

    runner = SolverRunner(command)
    runner.start(on_output=lambda line, _: output.append(line), on_complete=completed.append)

    deadline = time.monotonic() + 5
    while runner.running and time.monotonic() < deadline:
        time.sleep(0.01)

    assert completed
    assert completed[0].returncode == 0
    assert output == ["step 1", "step 2"]


def test_runner_rejects_empty_command() -> None:
    try:
        SolverRunner([])
    except ValueError:
        pass
    else:
        raise AssertionError("empty command must be rejected")


def test_runner_stop_terminates_process() -> None:
    command = [sys.executable, "-u", "-c", "import time; print('ready', flush=True); time.sleep(30)"]
    runner = SolverRunner(command)
    runner.start()
    deadline = time.monotonic() + 5
    while not runner.running and time.monotonic() < deadline:
        time.sleep(0.01)
    runner.stop()
    assert not runner.running

from pathlib import Path
import os
import sys
import time

import pytest

from cfdx import CFDXSession, ExecutionController, SolverRunner


def test_controller_propagates_metrics_and_success(tmp_path: Path) -> None:
    script = tmp_path / "solver.py"
    script.write_text("print('Iteration 4 Time = 0.25 CFL: 0.5')\n", encoding="utf-8")
    session = CFDXSession()
    runner = SolverRunner([sys.executable, str(script)])
    controller = ExecutionController(session, runner)
    controller.start()
    assert runner._thread is not None
    runner._thread.join(timeout=5)
    assert session.state.value == "CONVERGED"
    assert session.iteration == 4
    assert session.time == 0.25


def test_controller_propagates_failure(tmp_path: Path) -> None:
    script = tmp_path / "solver.py"
    script.write_text("raise SystemExit(3)\n", encoding="utf-8")
    session = CFDXSession()
    runner = SolverRunner(["python", str(script)])
    controller = ExecutionController(session, runner)
    controller.start()
    assert runner._thread is not None
    runner._thread.join(timeout=5)
    assert session.state.value == "FAILED"
    assert controller.error is not None
    assert controller.error.returncode == 3



def test_controller_stop_reports_stopped(tmp_path: Path) -> None:
    script = tmp_path / "solver.py"
    script.write_text("import time; time.sleep(30)\n", encoding="utf-8")
    session = CFDXSession()
    runner = SolverRunner(["python", str(script)])
    controller = ExecutionController(session, runner)
    controller.start()
    deadline = time.monotonic() + 5
    while not runner.running and time.monotonic() < deadline:
        time.sleep(0.01)
    assert runner.running
    controller.stop()
    assert session.state.value == "STOPPED"


@pytest.mark.skipif(os.name != "posix", reason="POSIX process groups provide pause semantics")
def test_controller_pause_resume(tmp_path: Path) -> None:
    script = tmp_path / "solver.py"
    script.write_text("import time; [print(i, flush=True) or time.sleep(0.05) for i in range(1000)]\n", encoding="utf-8")
    session = CFDXSession()
    runner = SolverRunner([sys.executable, "-u", str(script)])
    controller = ExecutionController(session, runner)
    controller.start()
    deadline = time.monotonic() + 5
    while not runner.running and time.monotonic() < deadline:
        time.sleep(0.01)
    assert runner.running
    controller.pause()
    assert session.state.value == "PAUSED"
    assert runner.paused
    controller.resume()
    assert session.state.value == "RUNNING"
    assert not runner.paused
    controller.stop()
    assert session.state.value == "STOPPED"

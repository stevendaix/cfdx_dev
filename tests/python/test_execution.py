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
    runner = SolverRunner([sys.executable, str(script)])
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
    runner = SolverRunner([sys.executable, str(script)])
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


def test_controller_restart_uses_configured_option(tmp_path: Path) -> None:
    script = tmp_path / "solver.py"
    script.write_text(
        "import sys\n"
        "assert sys.argv[1:] == ['--from-dat', 'checkpoint.dat']\n"
        "print('Iteration 12 Time = 3.5')\n",
        encoding="utf-8",
    )
    dat = tmp_path / "checkpoint.dat"
    dat.write_text("CFDX DAT restart\n", encoding="utf-8")
    session = CFDXSession()
    session.case.execution.restart_option = "--from-dat"
    runner = SolverRunner([sys.executable, str(script)], cwd=tmp_path)
    controller = ExecutionController(session, runner)
    controller.restart(dat)
    assert runner._thread is not None
    runner._thread.join(timeout=5)
    assert session.state.value == "CONVERGED"
    assert session.iteration == 12
    assert session.time == 3.5


def test_controller_builds_monitor_series_from_solver_metrics(tmp_path: Path) -> None:
    script = tmp_path / "solver.py"
    script.write_text(
        "print('Iteration 1 Time = 0.1 CFL: 0.5 residual p = 1.0e-2', flush=True)\\n"
        "print('Iteration 2 Time = 0.2 CFL: 0.3 residual p = 2.0e-3', flush=True)\\n",
        encoding="utf-8",
    )
    session = CFDXSession()
    runner = SolverRunner([sys.executable, str(script)])
    controller = ExecutionController(session, runner)
    controller.start()
    assert runner._thread is not None
    runner._thread.join(timeout=5)
    assert len(controller.monitor_series.samples) == 2
    assert controller.monitor_series.samples[-1].iteration == 2
    assert controller.monitor_series.samples[-1].time == pytest.approx(0.2)
    assert controller.monitor_series.at(2).values["p"] == pytest.approx(2.0e-3)
    assert controller.monitor_series.at(2).values["CFL"] == pytest.approx(0.3)

def test_controller_merges_metrics_emitted_on_separate_lines(tmp_path: Path) -> None:
    script = tmp_path / "solver.py"
    script.write_text(
        "print('Iteration 7 Time = 0.7', flush=True)\n"
        "print('CFL: 0.4', flush=True)\n"
        "print('residual p = 3.0e-4', flush=True)\n",
        encoding="utf-8",
    )
    session = CFDXSession()
    runner = SolverRunner([sys.executable, str(script)])
    controller = ExecutionController(session, runner)
    controller.start()
    assert runner._thread is not None
    runner._thread.join(timeout=5)
    assert len(controller.monitor_series.samples) == 1
    sample = controller.monitor_series.at(7)
    assert sample.time == pytest.approx(0.7)
    assert sample.values["CFL"] == pytest.approx(0.4)
    assert sample.values["p"] == pytest.approx(3.0e-4)

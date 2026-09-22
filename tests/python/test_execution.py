from pathlib import Path

from cfdx import CFDXSession, ExecutionController, SolverRunner


def test_controller_propagates_metrics_and_success(tmp_path: Path) -> None:
    script = tmp_path / "solver.py"
    script.write_text("print('Iteration 4 Time = 0.25 CFL: 0.5')\n", encoding="utf-8")
    session = CFDXSession()
    runner = SolverRunner(["python", str(script)])
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

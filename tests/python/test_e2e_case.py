from pathlib import Path
import sys

from cfdx import CFDXSession, ExecutionController, SolverRunner


def test_short_headless_case(tmp_path: Path) -> None:
    solver = tmp_path / "solver.py"
    solver.write_text(
        "print('Iteration 1 Time = 0.1 CFL: 0.2')\n"
        "print('Iteration 2 Time = 0.2 CFL: 0.2')\n",
        encoding="utf-8",
    )
    session = CFDXSession()
    controller = ExecutionController(session, SolverRunner([sys.executable, str(solver)]))
    controller.start()
    assert controller.runner._thread is not None
    controller.runner._thread.join(timeout=10)
    assert session.state.value == "CONVERGED"
    assert session.iteration == 2
    assert session.time == 0.2

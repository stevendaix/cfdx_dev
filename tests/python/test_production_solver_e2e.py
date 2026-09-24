from __future__ import annotations

import os
import tempfile
from pathlib import Path

from cfdx import CFDXSession, ExecutionController, SolverRunner
from cfdx.dat_io import read_dat_restart


def _run(controller: ExecutionController) -> None:
    output: list[str] = []
    controller.on_output = lambda line, is_stderr: output.append(
        ("stderr: " if is_stderr else "stdout: ") + line
    )
    controller.start()
    thread = controller.runner._thread
    assert thread is not None
    thread.join(timeout=30)
    assert not thread.is_alive()
    assert controller.session.state.value == "CONVERGED", (
        f"production solver failed: error={controller.error!r}; "
        f"output={output[-40:]!r}"
    )


def test_production_solver_full_application_e2e(tmp_path: Path) -> None:
    solver = os.environ.get("CFDX_PRODUCTION_SOLVER")
    mesh = os.environ.get("CFDX_PRODUCTION_MESH")
    assert solver and Path(solver).is_file()
    assert mesh and Path(mesh).is_file()

    first_dir = tmp_path / "first"
    session = CFDXSession()
    controller = ExecutionController(
        session,
        SolverRunner([
            solver, "--mesh", mesh,
            "--output-dir", str(first_dir),
            "--iterations", "20",
        ]),
    )
    _run(controller)

    assert controller.latest_metrics is not None
    # --iterations is the solver's maximum iteration count, not an exact
    # iteration target. The authoritative completion line reports the actual
    # nonlinear iteration reached by the production solver.
    assert controller.latest_metrics.iteration is not None
    assert 1 <= controller.latest_metrics.iteration <= 20
    assert controller.session.iteration == controller.latest_metrics.iteration

    checkpoint = first_dir / "restart.dat"
    assert checkpoint.is_file()

    restart = read_dat_restart(checkpoint)
    assert restart.cells > 0
    assert restart.iteration == controller.session.iteration
    assert 1 <= restart.iteration <= 20
    assert "U" in restart.fields
    assert "p" in restart.fields

    second_dir = tmp_path / "restart"
    restart_session = CFDXSession()
    restart_session.case.execution.restart_option = "--restart"
    restart_controller = ExecutionController(
        restart_session,
        SolverRunner([
            solver, "--mesh", mesh,
            "--output-dir", str(second_dir),
            "--iterations", "5",
        ]),
    )
    restart_controller.restart(checkpoint)
    thread = restart_controller.runner._thread
    assert thread is not None
    thread.join(timeout=30)
    assert not thread.is_alive()
    assert restart_session.state.value == "CONVERGED"
    assert restart_controller.latest_metrics is not None
    assert restart_controller.latest_metrics.iteration is not None
    assert 1 <= restart_controller.latest_metrics.iteration <= 5
    assert restart_session.iteration == restart_controller.latest_metrics.iteration

    outputs = sorted(second_dir.glob("result_*.vtu"))
    assert outputs
    xml = outputs[-1].read_text(encoding="utf-8")
    assert 'Name="physical_time"' in xml
    assert 'Name="iteration"' in xml


if __name__ == "__main__":
    with tempfile.TemporaryDirectory(prefix="cfdx-production-e2e-") as directory:
        test_production_solver_full_application_e2e(Path(directory))

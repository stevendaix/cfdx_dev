from __future__ import annotations

import os
from pathlib import Path

from cfdx import CFDXSession, ExecutionController, SolverRunner
from cfdx.dat_io import read_dat_restart


def _run(controller: ExecutionController) -> None:
    controller.start()
    thread = controller.runner._thread
    assert thread is not None
    thread.join(timeout=30)
    assert not thread.is_alive()
    assert controller.session.state.value == "CONVERGED"


def test_production_solver_full_application_e2e(tmp_path: Path) -> None:
    solver = os.environ.get("CFDX_PRODUCTION_SOLVER")
    mesh = os.environ.get("CFDX_PRODUCTION_MESH")
    assert solver and Path(solver).is_file()
    assert mesh and Path(mesh).is_file()

    first_dir = tmp_path / "first"
    session = CFDXSession()
    session.run()
    controller = ExecutionController(
        session,
        SolverRunner([
            solver, "--mesh", mesh,
            "--output-dir", str(first_dir),
            "--iterations", "2",
        ]),
    )
    _run(controller)

    assert controller.latest_metrics is not None
    assert controller.latest_metrics.iteration == 2
    checkpoint = first_dir / "restart.dat"
    assert checkpoint.is_file()

    restart = read_dat_restart(checkpoint)
    assert restart.cells > 0
    assert restart.iteration == 2
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
            "--iterations", "1",
        ]),
    )
    # ExecutionController appends the configured restart option and DAT path.
    # The solver therefore consumes the real production DAT, not a fixture.
    restart_controller.restart(checkpoint)
    thread = restart_controller.runner._thread
    assert thread is not None
    thread.join(timeout=30)
    assert not thread.is_alive()
    assert restart_session.state.value == "CONVERGED"

    outputs = sorted(second_dir.glob("result_*.vtu"))
    assert outputs
    xml = outputs[-1].read_text(encoding="utf-8")
    assert 'Name="physical_time"' in xml
    assert 'Name="iteration"' in xml

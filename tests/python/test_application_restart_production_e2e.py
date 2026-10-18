"""Application-routed DAT/restart qualification against a real production run.

The Application owns the DAT/restart contract through ``open_project`` and
``ensure_controller``: a loaded DAT checkpoint must become the solver's
restart state, and the restarted command must carry the checkpoint path. This
test qualifies that contract against an artifact actually emitted by the
production solver, not against a synthetic DAT.

The production solver and mesh are provided by ctest via environment variables,
so the test skips cleanly when the fixtures are unavailable.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest
from cfdx.session import CFDXSession

from cfdx.application import Application


def test_application_restart_command_carries_production_checkpoint(
    tmp_path: Path,
) -> None:
    solver = os.environ.get("CFDX_PRODUCTION_SOLVER")
    mesh = os.environ.get("CFDX_PRODUCTION_MESH")
    if not solver or not mesh:
        pytest.skip("CFDX_PRODUCTION_SOLVER/CFDX_PRODUCTION_MESH are provided by ctest")
    assert Path(solver).is_file(), f"production solver missing: {solver}"
    assert Path(mesh).is_file(), f"production mesh missing: {mesh}"

    run_dir = tmp_path / "restart_run"
    run_dir.mkdir()
    first = subprocess.run(
        [
            solver,
            "--mesh",
            mesh,
            "--output-dir",
            str(run_dir / "out1"),
            "--iterations",
            "5",
        ],
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
    )
    assert first.returncode == 0, first.stderr[-4000:]
    checkpoint = run_dir / "out1" / "restart.dat"
    assert checkpoint.is_file(), "production solver did not emit a restart checkpoint"

    # The Application owns the DAT/restart contract. A bare DAT checkpoint
    # written by the production solver must become the solver's restart state
    # through the Application, without requiring the caller to first save a
    # paired .cfdx.h5 project. The production solver requires --mesh and
    # --output-dir even when restarting from a checkpoint, so the bare DAT
    # restart carries those options alongside the checkpoint path.
    application = Application(CFDXSession())
    restart = application.load_restart(
        checkpoint, solver=solver, mesh=mesh, output_dir=str(run_dir / "out2")
    )
    assert restart.iteration >= 0

    controller = application.ensure_controller()
    command = controller.runner.command
    assert command[0] == solver
    assert "--mesh" in command
    assert "--output-dir" in command
    assert "--restart" in command
    assert str(checkpoint) in command

    # The restarted run must actually consume the checkpoint: the solver's
    # initialization strategy must be Restart, not the case-derived uniform
    # state, otherwise the checkpoint was loaded but ignored.
    restarted = subprocess.run(
        list(command),
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
        cwd=run_dir,
    )
    assert restarted.returncode == 0, restarted.stderr[-4000:]
    assert "Initialization strategy restart" in restarted.stdout

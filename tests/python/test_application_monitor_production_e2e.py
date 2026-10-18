"""Application-routed monitor history qualification against a real production run.

The Application owns the monitor contract through the ExecutionController: the
monitor series is the single place a frontend reads solver convergence evidence,
and it must carry the quantities the production solver actually reports rather
than only iteration progress. This test qualifies that contract against an
artifact actually emitted by the production solver.

The production solver and mesh are provided by ctest via environment variables,
so the test skips cleanly when the fixtures are unavailable.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest
from cfdx.execution import ExecutionController
from cfdx.runner import SolverRunner
from cfdx.session import CFDXSession


def test_application_monitor_carries_production_solver_convergence_evidence(
    tmp_path: Path,
) -> None:
    solver = os.environ.get("CFDX_PRODUCTION_SOLVER")
    mesh = os.environ.get("CFDX_PRODUCTION_MESH")
    if not solver or not mesh:
        pytest.skip("CFDX_PRODUCTION_SOLVER/CFDX_PRODUCTION_MESH are provided by ctest")
    assert Path(solver).is_file(), f"production solver missing: {solver}"
    assert Path(mesh).is_file(), f"production mesh missing: {mesh}"

    output_dir = tmp_path / "monitor_run"
    output_dir.mkdir()
    completed = subprocess.run(
        [
            solver,
            "--mesh",
            mesh,
            "--output-dir",
            str(output_dir),
            "--iterations",
            "5",
        ],
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr[-4000:]

    # The Application owns the monitor contract through the ExecutionController:
    # the monitor series is the single place a frontend reads solver
    # convergence evidence, so the contract is exercised here through the same
    # controller the Application would construct.
    session = CFDXSession()
    controller = ExecutionController(
        session,
        SolverRunner(
            [
                solver,
                "--mesh",
                mesh,
                "--output-dir",
                str(output_dir),
                "--iterations",
                "5",
            ]
        ),
    )
    controller.start()
    controller.runner._thread.join(timeout=60.0)

    assert session.state.value == "CONVERGED"
    series = controller.monitor_series
    assert series.samples, "production run produced no monitor samples"

    # The monitor series must carry the production solver's convergence
    # vocabulary, not just iteration progress. The solver reports
    # momentum_residual_raw and momentum_residual_reference on every trace
    # line; the Application's monitor contract must surface them.
    last = series.samples[-1]
    assert last.iteration >= 1
    assert "momentum_residual_raw" in last.values
    assert "momentum_residual_reference" in last.values

    # The persisted monitor JSON must round-trip through the same contract
    # the Application exposes, preserving iteration/time provenance.
    json_path = tmp_path / "monitor.json"
    series.write_json(json_path)
    reloaded = type(series).read_json(json_path)
    assert reloaded.name == series.name
    assert len(reloaded.samples) == len(series.samples)
    assert reloaded.samples[-1].iteration == series.samples[-1].iteration
    assert reloaded.samples[-1].values == series.samples[-1].values

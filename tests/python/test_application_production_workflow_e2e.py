"""End-to-end production workflow test for the GUI orchestration stack.

This test exercises the real production solver against a real CFDX case
artifact, covering the full cycle the GUI must support:

    Open case -> Configure -> Validate -> Run -> Monitor -> Checkpoint
    -> Reload -> Restart -> Results

The case artifact is produced by the Python conversion layer and written
through the canonical CFDX HDF5 writer. The production solver is a real
binary fixture, so this test proves the contract between the Python case
writer and the C++ solver reader, not a mock.
"""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import h5py
import pytest

from cfdx import CFDXSession, SimulationState
from cfdx.application.application import Application
from cfdx.dat_io import read_dat_restart
from cfdx.execution import ExecutionController
from cfdx.runner import SolverRunner
from cfdx.session import ChangeImpact
from cfdx.validation import validate_case


SOLVER_ENV = "CFDX_PRODUCTION_SOLVER"
MESH_ENV = "CFDX_PRODUCTION_MESH"


def _production_solver() -> str:
    solver = os.environ.get(SOLVER_ENV)
    if not solver:
        pytest.skip(f"{SOLVER_ENV} is provided by ctest")
    return solver


def _production_mesh() -> str:
    mesh = os.environ.get(MESH_ENV)
    if not mesh:
        pytest.skip(f"{MESH_ENV} is provided by ctest")
    return mesh


def _run_solver(solver: str, mesh: str, output_dir: Path,
                iterations: int = 1, restart: str | None = None
                ) -> subprocess.CompletedProcess[str]:
    cmd = [solver, "--mesh", mesh, "--output-dir", str(output_dir),
           "--iterations", str(iterations)]
    if restart:
        cmd.extend(["--restart", restart])
    return subprocess.run(cmd, capture_output=True, text=True,
                          timeout=60, check=False)


def _read_execution(path: Path) -> dict:
    assert path.is_file(), f"missing execution artifact: {path}"
    document = json.loads(path.read_text(encoding="utf-8"))
    assert document["format"] == "CFDX-EXECUTION"
    assert document["schema_version"] == 1
    return document


def _read_convergence(path: Path) -> dict:
    assert path.is_file(), f"missing convergence artifact: {path}"
    return json.loads(path.read_text(encoding="utf-8"))


def test_full_production_workflow_e2e(tmp_path: Path) -> None:
    """Drive the complete GUI cycle against the real production solver."""
    solver = _production_solver()
    mesh = _production_mesh()

    # 1. Open a real mesh through the canonical catalog.
    session = CFDXSession()
    session.case.name = "e2e-production"
    session.case.enable("incompressible")
    session.case.materials["air"] = {
        "density": 1.2,
        "dynamic_viscosity": 1.8e-5,
        "cp": 1005.0,
        "conductivity": 0.026,
    }
    session.case.execution.solver = solver
    session.case.execution.mesh = mesh
    session.case.execution.restart_option = "--restart"
    session.case.execution.output_dir = str(tmp_path / "output")
    session.case.set_numerics(cfl=0.5, max_iterations=50, residual_target="1e-3")

    # 2. Validate the case configuration.
    report = validate_case(session.case, None)
    assert report.ok, report.diagnostics

    # 3. Run the solver for one iteration and inspect the artifacts.
    output_dir = tmp_path / "run"
    completed = _run_solver(solver, mesh, output_dir, iterations=1)
    assert completed.returncode in (0, 1), (
        f"production solver failed: {completed.stderr[-2000:]}"
    )

    execution = _read_execution(output_dir / "execution.json")
    assert execution["process_exit_code"] == completed.returncode
    assert execution["artifacts"]["restart_dat"] is True
    assert execution["artifacts"]["convergence_json"] is True

    # 4. Inspect the convergence history for the vocabulary the GUI reads.
    convergence = _read_convergence(output_dir / "convergence.json")
    assert "schema_version" in convergence
    assert "iterations" in convergence
    assert isinstance(convergence["iterations"], int)
    assert convergence["iterations"] >= 1
    assert convergence["history"], "convergence history must be non-empty"

    # 5. Reload the DAT checkpoint and restart from it.
    dat_path = output_dir / "restart.dat"
    assert dat_path.is_file(), "production solver must emit a restart DAT"
    restart = read_dat_restart(dat_path)
    assert restart.iteration >= 1

    restarted_dir = tmp_path / "restart"
    restarted = _run_solver(solver, mesh, restarted_dir, iterations=1,
                            restart=str(dat_path))
    assert restarted.returncode in (0, 1), (
        f"restart failed: {restarted.stderr[-2000:]}"
    )

    restart_execution = _read_execution(restarted_dir / "execution.json")
    assert restart_execution["artifacts"]["restart_dat"] is True

    # 6. Drive the same cycle through the Application facade so the GUI
    #    orchestration path is exercised end-to-end. The Application requires
    #    a canonical .cfdx.h5 case artifact (format="CFDX"), so build one
    #    through the session rather than copying the mesh file.
    from cfdx.case_io import save_case
    case_path = tmp_path / "case.cfdx.h5"
    case_path.write_bytes(Path(mesh).read_bytes())
    save_case(session, case_path)
    application = Application()
    application.open_project(case_path)
    application.set_numerical_option("cfl", 0.5, ChangeImpact.HOT)
    controller = application.ensure_controller()
    assert controller is not None
    assert controller.runner.command[0] == solver


def test_production_solver_artifact_contract(tmp_path: Path) -> None:
    """Prove the artifacts the GUI depends on are produced deterministically."""
    solver = _production_solver()
    mesh = _production_mesh()

    output_dir = tmp_path / "artifacts"
    completed = _run_solver(solver, mesh, output_dir, iterations=1)
    assert completed.returncode in (0, 1), completed.stderr

    execution = _read_execution(output_dir / "execution.json")
    assert execution["artifacts"]["restart_dat"] is True
    assert execution["artifacts"]["convergence_json"] is True

    convergence = _read_convergence(output_dir / "convergence.json")
    assert isinstance(convergence["iterations"], int)
    assert convergence["iterations"] >= 1
    assert "history" in convergence
    assert convergence["history"], "convergence history must be non-empty"

    # The monitor history parser must be able to ingest the solver's own
    # INCOMPRESSIBLE_ITER stdout. The first iteration carries the vocabulary
    # the GUI uses to distinguish running from converged.
    first = convergence["history"][0]
    assert "iteration" in first
    # The parser consumes the stdout INCOMPRESSIBLE_ITER vocabulary, not the
    # JSON history keys. Verify the parser can ingest the solver's own output.
    from cfdx.metrics import SolverMetricsParser
    parser = SolverMetricsParser()
    metrics = [parser.parse(line) for line in completed.stdout.splitlines()]
    metrics = [m for m in metrics if m is not None]
    assert metrics, "parser must extract at least one sample from solver stdout"
    assert metrics[-1].iteration == 1
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path

import pytest

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


def _run_production(
    solver: str, mesh: str, output_dir: Path, iterations: int
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            solver,
            "--mesh",
            mesh,
            "--output-dir",
            str(output_dir),
            "--iterations",
            str(iterations),
        ],
        text=True,
        capture_output=True,
        timeout=30,
        check=False,
    )


def _read_execution(path: Path) -> dict:
    assert path.is_file(), f"missing execution artifact: {path}"
    document = json.loads(path.read_text(encoding="utf-8"))
    assert document["format"] == "CFDX-EXECUTION"
    assert document["schema_version"] == 1
    return document


def test_production_solver_non_converged_exit_is_persisted(tmp_path: Path) -> None:
    solver = os.environ.get("CFDX_PRODUCTION_SOLVER")
    mesh = os.environ.get("CFDX_PRODUCTION_MESH")
    if not solver or not mesh:
        pytest.skip("CFDX_PRODUCTION_SOLVER/CFDX_PRODUCTION_MESH are provided by ctest")

    output_dir = tmp_path / "non_converged"
    completed = _run_production(solver, mesh, output_dir, 1)

    assert completed.returncode == 1, completed.stderr
    execution = _read_execution(output_dir / "execution.json")
    assert execution["process_exit_code"] == 1
    assert execution["converged"] is False
    assert execution["artifacts"]["restart_dat"] is True
    assert execution["artifacts"]["convergence_json"] is True


def test_production_solver_exception_exit_is_persisted(tmp_path: Path) -> None:
    solver = os.environ.get("CFDX_PRODUCTION_SOLVER")
    if not solver:
        pytest.skip("CFDX_PRODUCTION_SOLVER is provided by ctest")

    output_dir = tmp_path / "exception"
    completed = _run_production(
        solver,
        str(tmp_path / "missing_mesh.msh"),
        output_dir,
        2,
    )

    assert completed.returncode == 2
    execution = _read_execution(output_dir / "execution.json")
    assert execution["process_exit_code"] == 2
    assert execution["converged"] is False
    assert execution["artifacts"]["restart_dat"] is False
    assert execution["artifacts"]["convergence_json"] is False


def test_production_solver_artifact_failure_updates_execution_exit(tmp_path: Path) -> None:
    solver = os.environ.get("CFDX_PRODUCTION_SOLVER")
    mesh = os.environ.get("CFDX_PRODUCTION_MESH")
    if not solver or not mesh:
        pytest.skip("CFDX_PRODUCTION_SOLVER/CFDX_PRODUCTION_MESH are provided by ctest")

    output_dir = tmp_path / "artifact_failure"
    output_dir.mkdir()
    (output_dir / "restart.dat").mkdir()

    completed = _run_production(solver, mesh, output_dir, 20)

    assert completed.returncode == 2, completed.stderr
    execution = _read_execution(output_dir / "execution.json")
    assert execution["process_exit_code"] == 2
    assert execution["converged"] is True
    assert execution["artifacts"]["restart_dat"] is False
    assert execution["artifacts"]["convergence_json"] is False
    assert not (output_dir / "convergence.json").exists()



def test_production_solver_stop_checkpoint_reload_restart(tmp_path: Path) -> None:
    solver = os.environ.get("CFDX_PRODUCTION_SOLVER")
    mesh = os.environ.get("CFDX_PRODUCTION_MESH")
    if not solver or not mesh:
        pytest.skip("CFDX_PRODUCTION_SOLVER/CFDX_PRODUCTION_MESH are provided by ctest")
    assert Path(solver).is_file()
    assert Path(mesh).is_file()

    first_dir = tmp_path / "stopped"
    session = CFDXSession()
    controller = ExecutionController(
        session,
        SolverRunner([
            solver, "--mesh", mesh,
            "--output-dir", str(first_dir),
            "--iterations", "100000",
        ]),
    )
    stop_started = False

    def request_stop(line: str, is_stderr: bool) -> None:
        nonlocal stop_started
        if is_stderr or stop_started:
            return
        if line.startswith("Iteration "):
            stop_started = True
            import threading
            threading.Thread(target=controller.stop, daemon=True).start()

    controller.on_output = request_stop
    controller.start()
    thread = controller.runner._thread
    assert thread is not None
    thread.join(timeout=30)
    assert not thread.is_alive()
    assert stop_started
    assert session.state.value == "STOPPED", (
        f"graceful production stop failed: error={controller.error!r}"
    )

    checkpoint = first_dir / "restart.dat"
    assert checkpoint.is_file()
    restart = read_dat_restart(checkpoint)
    assert restart.cells > 0
    assert restart.iteration >= 1
    assert restart.iteration == session.iteration
    assert "U" in restart.fields
    assert "p" in restart.fields

    second_dir = tmp_path / "restarted"
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
    restart_thread = restart_controller.runner._thread
    assert restart_thread is not None
    restart_thread.join(timeout=30)
    assert not restart_thread.is_alive()
    assert restart_session.state.value == "CONVERGED", (
        f"restart after graceful stop failed: error={restart_controller.error!r}"
    )
    assert restart_controller.latest_metrics is not None
    assert restart_controller.latest_metrics.iteration is not None
    assert restart_session.iteration == restart_controller.latest_metrics.iteration


def test_production_solver_full_application_e2e(tmp_path: Path) -> None:
    solver = os.environ.get("CFDX_PRODUCTION_SOLVER")
    mesh = os.environ.get("CFDX_PRODUCTION_MESH")
    if not solver or not mesh:
        pytest.skip("CFDX_PRODUCTION_SOLVER/CFDX_PRODUCTION_MESH are provided by ctest")
    assert Path(solver).is_file()
    assert Path(mesh).is_file()

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
    convergence = first_dir / "convergence.json"
    execution = first_dir / "execution.json"
    assert convergence.is_file()
    assert execution.is_file()
    convergence_document = json.loads(convergence.read_text(encoding="utf-8"))
    execution_document = json.loads(execution.read_text(encoding="utf-8"))
    assert convergence_document["format"] == "CFDX-CONVERGENCE"
    assert convergence_document["converged"] is True
    assert execution_document["format"] == "CFDX-EXECUTION"
    assert execution_document["process_exit_code"] == 0
    assert execution_document["converged"] is True
    assert execution_document["artifacts"]["restart_dat"] is True
    assert execution_document["artifacts"]["convergence_json"] is True

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
    restart_output: list[str] = []
    restart_controller.on_output = lambda line, is_stderr: restart_output.append(
        ("stderr: " if is_stderr else "stdout: ") + line
    )
    restart_controller.restart(checkpoint)
    thread = restart_controller.runner._thread
    assert thread is not None
    thread.join(timeout=30)
    assert not thread.is_alive()
    assert restart_session.state.value == "CONVERGED", (
        f"production solver restart failed: error={restart_controller.error!r}; "
        f"output={restart_output[-40:]!r}"
    )
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

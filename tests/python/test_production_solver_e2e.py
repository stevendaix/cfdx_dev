from __future__ import annotations

import os
import tempfile
from pathlib import Path

from cfdx import CFDXSession, ExecutionController, SolverRunner
from cfdx.io.converter import convert
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


def test_production_solver_explicit_cfdx_case_e2e(tmp_path: Path) -> None:
    solver = os.environ.get("CFDX_PRODUCTION_SOLVER")
    assert solver and Path(solver).is_file()

    root = Path(__file__).resolve().parents[2]
    source = root / "tests" / "data" / "su2" / "mesh_NACA0012_inv.su2"
    config = root / "tests" / "data" / "su2" / "inv_NACA0012_basic.cfg"
    assert source.is_file()
    assert config.is_file()

    case_dir = tmp_path / "explicit_case"
    case_dir.mkdir()
    source_case = case_dir / source.name
    config_case = case_dir / "case.cfg"
    source_case.write_bytes(source.read_bytes())
    config_case.write_bytes(config.read_bytes())
    case_path = case_dir / "case.cfdx.h5"

    converted = convert(source_case, output=case_path)
    assert converted.case is not None
    assert converted.case.numerics.selection.entries
    assert converted.case.numerics.selection.required_families

    selected = {
        (entry.family, entry.configuration_key)
        for entry in converted.case.numerics.selection.entries
    }
    assert ("gradient", "numerics.gradient.gauss") in selected
    assert ("convection", "numerics.convection.upwind") in selected
    assert ("pressure_velocity", "pressure_velocity.simple") in selected

    output_dir = tmp_path / "explicit_run"
    controller = ExecutionController(
        CFDXSession(),
        SolverRunner([
            solver, "--mesh", str(case_path),
            "--output-dir", str(output_dir),
            "--iterations", "20",
        ]),
    )
    output: list[str] = []
    controller.on_output = lambda line, is_stderr: output.append(
        ("stderr: " if is_stderr else "stdout: ") + line
    )
    controller.start()
    thread = controller.runner._thread
    assert thread is not None
    thread.join(timeout=30)
    assert not thread.is_alive(), f"solver timed out: output={output[-40:]!r}"
    assert controller.session.state.value == "CONVERGED", (
        f"explicit CFDX case failed: error={controller.error!r}; "
        f"output={output[-40:]!r}"
    )

    diagnostics = "".join(output)
    assert "Resolved numerical selections:" in diagnostics
    assert "scheme[gradient]=numerics.gradient.gauss" in diagnostics
    assert "scheme[convection]=numerics.convection.upwind" in diagnostics
    assert "scheme[pressure_velocity]=pressure_velocity.simple" in diagnostics
    assert "scheme[linear_solver]=linear.fgmres" in diagnostics


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

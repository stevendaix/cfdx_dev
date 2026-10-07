from __future__ import annotations

import os
import stat
import time
from pathlib import Path
from typing import Any

import pytest
from cfdx import CFDXSession


pytestmark = pytest.mark.skipif(
    os.environ.get("CFDX_PRODUCTION_SOLVER") is None
    or os.environ.get("CFDX_PRODUCTION_MESH") is None,
    reason="production solver GUI E2E requires CFDX_PRODUCTION_SOLVER and CFDX_PRODUCTION_MESH",
)


def _wait_for_convergence(window: Any, timeout: float = 45.0) -> None:
    from cfdx.session import SimulationState
    from PySide6.QtWidgets import QApplication

    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        QApplication.processEvents()
        if window.session.state is SimulationState.CONVERGED:
            return
        if window.session.state is SimulationState.FAILED:
            raise AssertionError("production solver entered FAILED state")
        time.sleep(0.05)
    raise AssertionError(
        f"production solver did not converge within {timeout}s; "
        f"state={window.session.state}, iteration={window.session.iteration}"
    )


def _write_solver_adapter(path: Path, output_dir: Path) -> None:
    solver = Path(os.environ["CFDX_PRODUCTION_SOLVER"]).resolve()
    mesh = Path(os.environ["CFDX_PRODUCTION_MESH"]).resolve()
    path.write_text(
        "#!/usr/bin/env python3\n"
        "import subprocess, sys\n"
        f"solver = {str(solver)!r}\n"
        f"mesh = {str(mesh)!r}\n"
        f"output_dir = {str(output_dir)!r}\n"
        "args = [solver, '--mesh', mesh, '--output-dir', output_dir, '--iterations', '5']\n"
        "restart = list(sys.argv[2:])\n"
        "if restart:\n"
        "    args.extend(restart)\n"
        "raise SystemExit(subprocess.call(args))\n",
        encoding="utf-8",
    )
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


@pytest.mark.skipif(
    __import__("importlib.util").util.find_spec("PySide6") is None,
    reason="PySide6 optional",
)
def test_workbench_production_solver_and_3d_result_e2e(tmp_path: Path) -> None:
    from cfdx.results_series import discover_result_series
    from cfdx.session import SimulationState
    from cfdx.workbench import CFDXWorkbenchWindow
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication(["cfdx-p2i-e2e"])
    case_path = tmp_path / "production.cfdx.h5"
    results_dir = tmp_path / "results"
    adapter = tmp_path / "production_solver_adapter.py"
    _write_solver_adapter(adapter, results_dir)

    session = CFDXSession()
    session.case.name = "p2i-production"
    session.case.execution.solver = str(adapter)
    session.case.set_numerics(cfl=0.5)

    window = CFDXWorkbenchWindow(session)
    try:
        window.application.save_project(case_path)
        assert case_path.is_file()

        window._run()
        _wait_for_convergence(window)

        assert window.session.state is SimulationState.CONVERGED
        assert window.session.iteration > 0

        outputs = sorted(results_dir.glob("result_*.vtu"))
        assert outputs, "production solver produced no VTU result"
        series = discover_result_series(results_dir, inspect_fields=True)
        assert series.frames
        frame = series.frames[-1]
        assert frame.complete

        assert window.view3d is not None
        window.view3d.load(str(frame.path))
        assert window.view3d._current_result == frame.path
        dataset = window.view3d._dataset()
        assert dataset.n_cells > 0
        assert dataset.n_points > 0

        window.application.open_results(results_dir)
        assert window.application.state.results.frames
        window.view3d.set_field("p")
        assert window.view3d._selected_field == "p"
    finally:
        window.close()
        app.processEvents()
        app.quit()
